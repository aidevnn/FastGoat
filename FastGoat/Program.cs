using System.CodeDom;
using System.Collections;
using System.Diagnostics;
using System.Numerics;
using FastGoat.Commons;
using FastGoat.Examples;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.Naming;
using FastGoat.Structures.Subgroups;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Characters;
using FastGoat.UserGroup.GModuleN;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Matrix;
using FastGoat.UserGroup.Perms;
using FastGoat.UserGroup.Polynoms;
using FastGoat.UserGroup.Words;
using static FastGoat.Commons.IntExt;
using static FastGoat.Commons.EnumerableExt;
using FastGoat.UserGroup.Padic;
using FastGoat.UserGroup.Words.Tools;
using GroupRegX = System.Text.RegularExpressions;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");

ConcreteGroup<Mat> ProductDiagPerm(ConcreteGroup<Mat> group0, ConcreteGroup<Mat> group1)
{
    var (gl0, gl1) = (group0.Neutral().GL, group1.Neutral().GL);
    var dim = gl0.N + gl1.N;
    var (p0, p1) = (gl0.P, gl1.P);
    var p = Primes10000.First(p => (p - 1) % (p0 - 1) == 0 && (p - 1) % (p1 - 1) == 0);
    var Up0 = FG.UnInt(p0);
    var Up1 = FG.UnInt(p1);
    var Up = FG.UnInt(p);
    var iso0 = Group.AllMorphisms(Up0, Up, Group.MorphismType.Isomorphism)
        .First()
        .HomMap.ToDictionary(e => e.Key.K, e => e.Value.K);
    var iso1 = Group.AllMorphisms(Up1, Up, Group.MorphismType.Isomorphism)
        .First()
        .HomMap.ToDictionary(e => e.Key.K, e => e.Value.K);

    iso0[0] = 0;
    iso1[0] = 0;

    var gl = new GL(dim, p);
    var gens0 = group0.GetGenerators().Select(e => e.Table.Select(k => iso0[k]).ToArray())
        .Select(e => (e, gl0.N)).ToArray();
    var gens1 = group1.GetGenerators().Select(e => e.Table.Select(k => iso1[k]).ToArray())
        .Select(e => (e, gl1.N)).ToArray();

    var Gens0 = gens0.Select(e => MatrixExt.MergeDiagonalBlocks(e, (gl1.Neutral().Table, gl1.N))).Select(e => gl.Create(e)).ToArray();
    var Gens1 = gens1.Select(e => MatrixExt.MergeDiagonalBlocks((gl0.Neutral().Table, gl0.N), e)).Select(e => gl.Create(e)).ToArray();
    return Group.Generate($"{group0.NameParenthesis()} x {group1.NameParenthesis()}", gl, Gens0.Concat(Gens1).ToArray());
}

(ConcreteGroup<Mat> mat, bool isDiagByPerm) MatrixFormFromNames(ANameElt name)
{
    if (name is Leaf leaf)
    {
        var (pr, coefs) = leaf.LeafDetails();
        if (pr == "C")
            return (FG.AbelianMat(coefs), true);
        else if (pr == "Q" || pr == "Dic")
        {
            var m = pr == "Dic" ? coefs[0] : coefs[0] / 4;
            return (FG.DicyclicGL2p(m), true);
        }
        else if (pr == "A" || pr == "S")
        {
            if (coefs[0] == 4)
            {
                var gl = new GL(3, 3);
                var (a, b) = pr == "S"
                    ? (gl[0, 1, 0, 0, 0, 1, 1, 0, 0], gl[1, 0, 0, 0, 0, 1, 0, 2, 0])
                    : (gl[0, 1, 0, 0, 0, 1, 1, 0, 0], gl[0, 0, 1, 2, 0, 0, 0, 2, 0]);
                return (Group.Generate(name.Name, gl, a, b), true);
            }

            if (coefs[0] == 5)
            {
                var so35 = FG.SO3p(5);
                if (pr == "S")
                    return (so35, false);
                else
                    return (Group.IsomorphicSubgroup(so35, FG.Alternate(5)), false);
            }
        }
        else if (pr == "GL" || pr == "SL")
        {
            var gMat = pr == "GL" ? FG.GLnp(coefs[0], coefs[1]) : FG.SLnp(coefs[0], coefs[1]);
            return (gMat, false);
        }
    }
    else if (name is SemiDirectProductOp sdpOp)
    {
        var mtCyc = sdpOp.MetaCyclicDetails();
        if (mtCyc.Length == 3)
            return (FG.MetaCyclicSdpMat(mtCyc[0], mtCyc[1], mtCyc[2]), true);
    }
    else if (name is DirectProductOp dpOp)
    {
        var eltsMat = dpOp.Elts.Select(name0 => MatrixFormFromNames(name0)).ToArray();
        if (eltsMat.All(e => e.isDiagByPerm) || name.Name == "C2 x SL(2,3)")
        {
            var mat0 = ProductDiagPerm(eltsMat[0].mat, eltsMat[1].mat);
            foreach (var (mat, _) in eltsMat.Skip(2))
                mat0 = ProductDiagPerm(mat0, mat);

            return (mat0, true);
        }
    }

    return (Group.Generate(new GL(1, 2)), false);
}

void MatrixFormTinyGroups(int maxOrder)
{
    var total = 0;
    var missing = new List<ANameElt[]>();
    foreach (var g in maxOrder.Range(1).SelectMany(o => FG.AllGroupsOfOrder(o)))
    {
        ++total;
        var gSubgrs = g.AllSubgroups().ToGroupWrapper();
        var names = NamesTree.BuildName(gSubgrs);
        var (mat, check) = MatrixFormFromNames(names[0]);
        
        if (mat.Count() == 1 && !check)
            missing.Add(names);
        else
        {
            FG.DisplayName(gSubgrs.Parent, gSubgrs, names, false, false, 20);
            DisplayGroup.Generators(mat);
            Console.WriteLine();

            if (!mat.IsIsomorphicTo(g))
                throw new();
        }
    }

    missing.Println(e => e[0].ContentGroup!.ShortName, $"Missing:{missing.Count} Found:{total - missing.Count}/{total}");
    missing.Where(e => e[0].ContentType == ANameElt.NodeType.DirectProduct).Println(e => e[0].ContentGroup!.ShortName, "Missing Direct Product");
    missing.SelectMany(e => e.Where(f => f is ExtensionOp f0 && f0.Lhs.ContentGroup!.GetGenerators().Count() == 1 && f0.Rhs.ContentGroup!.GetGenerators().Count() == 1).Take(1))
        .Println(e => $"{e.Name} -> {e.ContentGroup!.ShortName}", "Missing Non Split Metacyclic");
}

IEnumerable<AllSubgroups<Mat>> MtCyclicSubgroups(int order, int factor = 4)
{
    var total = 0;
    foreach (var g in factor.Range(1).SelectMany(k => FG.MetaCyclicSdpMat(k * order)))
    {
        ++total;
        var gSubgrs = g.AllSubgroups();
        foreach (var sg in gSubgrs.Where(cj => cj.Order == order).Select(cj => gSubgrs.Restriction(cj.Representative)))
            yield return sg;
    }
}

void TestProdAbMtCyc()
{
    var c2c4 = FG.AbelianMat(2, 4);
    var d8 = FG.DihedralGL2p(4);
    var c2c4d8 = ProductDiagPerm(c2c4, d8);
    DisplayGroup.HeadNames(c2c4d8);
    DisplayGroup.Generators(c2c4d8);

    var c3c3 = FG.AbelianMat(3, 3);
    var dic3 = FG.DicyclicGL2p(3);
    var c3c3dic3 = ProductDiagPerm(c3c3, dic3);
    DisplayGroup.HeadNames(c3c3dic3);
    DisplayGroup.Generators(c3c3dic3);
}

{
    Group.ActivedStorage(false);
    MatrixFormTinyGroups(32);
    // MtCyclicSubgroups(32, factor: 5).Select(e => e.ToGroupWrapper()).FilterIsomorphic().Naming().DisplayNames();
    // MtCyclicSubgroups(48, factor: 5).Select(e => e.ToGroupWrapper()).FilterIsomorphic().Naming().DisplayNames();
    // MtCyclicSubgroups(64, factor: 5).Select(e => e.ToGroupWrapper()).FilterIsomorphic().Naming().DisplayNames();
}

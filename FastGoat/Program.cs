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

#region Method1

ConcreteGroup<Mat> ProductMatrixBlock(ConcreteGroup<Mat> group0, ConcreteGroup<Mat> group1)
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
            var mat0 = ProductMatrixBlock(eltsMat[0].mat, eltsMat[1].mat);
            foreach (var (mat, _) in eltsMat.Skip(2))
                mat0 = ProductMatrixBlock(mat0, mat);

            return (mat0, true);
        }
    }

    return (Group.Generate(new GL(1, 2)), false);
}

#endregion

#region Method2

IEnumerable<Mat[]> MGenerators(int p, int[] type, int dim)
{
    if ((dim < 4 && BigInteger.Pow(p, dim) > 10000) || BigInteger.Pow(p, dim) > 100000)
    {
        yield break;
    }

    var sz = type.Length;
    var m = Gcd(type);
    var Up = FG.UnInt(p);
    var sn = FG.Symmetric(dim);
    var gl = new GL(dim, p);

    var byOrders = Up.MultiLoop(dim).Select(l => gl.Create(MatrixExt.Diagonal(l.Select(z => z.K).ToArray())))
                    .GroupBy(mat => Group.Cycle(gl, mat).Count)
                    .ToDictionary(e => e.Key, e => e.ToArray());

    if (type.Any(k => !byOrders.ContainsKey(k)))
        yield break;

    foreach(var mat in type.Select(f => byOrders[f]).MultiLoop().Select(ms => ms.ToArray()))
        yield return mat;
}

IEnumerable<(int[] perm, int[][] cycles, Mat mat)> CGenerators(int m, int n, int dim)
{
    var allTypes = IntExt.Partitions32[dim].Select(l => l.Order().ToArray()).OrderBy(l => l.Length).ToArray();
    var nks = allTypes.Select(l => l.Aggregate((a0, a1) => a0 * a1))
        .SelectMany(e => IntExt.Dividors(e).Append(e).Where(j => j != 1)).Append(n).ToHashSet();
    foreach (var p in nks.SelectMany(nk => IntExt.Primes10000.Where(p => (p - 1) % m == 0 && (p - 1) % nk == 0).Take(6)).Distinct().Order().Where(p => p < 62))
    {
        if ((dim < 4 && BigInteger.Pow(p, dim) > 10000) || BigInteger.Pow(p, dim) > 100000)
            continue;

        var Up = FG.UnInt(p);
        var gl = new GL(dim, p);
        var matn = Up.MultiLoop(dim).Select(l => gl.Create(MatrixExt.Diagonal(l.Select(e => e.K).ToArray()))).ToArray();
        // Console.WriteLine($"{gl} press key...");
        // Console.ReadLine();
        var sn = new Sn(dim);
        var m1s = IntExt.Partitions32[dim]//.Where(l => l.Count == l.Distinct().Count())
            .OrderBy(l => l.Count)
            .Select(t => IntExt.PermAndCyclesFromType(t.Order().ToArray()))
            .Select(e =>
            {
                var e0 = gl.Neutral().Table.Chunk(dim).ToArray();
                var perm = sn.CreateElement(e.perm.Select(i => i + 1).ToArray());
                var e1 = perm.Apply(e0);
                var mat0 = gl.Create(e1.SelectMany(v => v).ToArray());
                return matn.Select(mat => gl.Op(mat0, mat))
                    .Where(mat => mat.IsOrder(n))
                    .Select(mat => (e.perm, e.cycles, mat));
            })
            .SelectMany(e => e);

        foreach(var e in m1s)
            yield return e;
    }
}

(Word[] Mgens, Word Cgen, string name)[] ExtractGeneratorsSdp(ANameElt[] names)
{
    foreach(SemiDirectProductOp e0 in names.Where(e => e is SemiDirectProductOp e0 &&
            e0.Lhs.ContentGroup!.GroupType == GroupType.AbelianGroup &&
            e0.Rhs.ContentGroup!.GroupType == GroupType.AbelianGroup))
    {
        var Mgens = e0.Lhs.ContentGroup!.GetGenerators().Select(m0 => (Word)m0.E)
            .Select(e => char.IsUpper(e.Get()[0]) ? new Word(e.WGroup, e.Get().Revert()) : e).ToArray();
        var Cgens = e0.Rhs.ContentGroup!.GetGenerators().Select(m0 => (Word)m0.E)
            .Select(e => char.IsUpper(e.Get()[0]) ? new Word(e.WGroup, e.Get().Revert()) : e).ToArray();

        if (Cgens.Length == 1 && Cgens.Concat(Mgens).All(w => w.Get().Length == 1 && w.Get().All(c => char.IsLower(c))))
            return new (Word[] Mgens, Word Cgen, string name)[1] {(Mgens, Cgens[0], e0.Name)};
    }

    return new (Word[] Mgens, Word Cgen, string name)[0];
}

ConcreteGroup<Mat> MatrixFormFromNamesMeth2(WordGroup g, ANameElt[] names)
{
    var ext = ExtractGeneratorsSdp(names);
    if (ext.Length == 0)
        return Group.Generate(new GL(1,2));

    var (Mgens, Cgen, name) = ext[0];
    // Console.WriteLine($"{g.ShortName} {g.Definition}");
    // Mgens.Append(Cgen).Select(w => (w, g.ElementsOrders[w])).Println(name);

    var type = Mgens.Select(e => g.ElementsOrders[e]).ToArray();
    var m = Gcd(type);
    var n = g.ElementsOrders[Cgen];
    foreach (var dim in 4.Range(1).Where(d => d >= type.Length && (IntExt.Gcd(m, d) != 1 || IntExt.Gcd(m - 1, d) != 1)))
    {
        foreach(var (perm, cycles, m1) in CGenerators(m, n, dim))
        {
            var gl = m1.GL;
            var p = gl.P;
            // Console.WriteLine($"{g} in {gl} type:[{type.Glue(" ")}] cycles:{cycles.Select(c => c.Glue(" ")).Glue("","({0})")} m:{m} n:{n}");
            // Console.WriteLine($"m1:{m1}");
            foreach(var m0s in MGenerators(p, type, dim))
            {
                var map = Mgens.Zip(m0s).ToDictionary(e => e.First.Get()[0], e => e.Second);
                map[Cgen.Get()[0]] = m1;
                if (g.CheckHomomorphism(gl, map))
                {
                    // map.Println("Gens");
                    // Console.ReadLine();
                    var mat = Group.Generate(g.Name, gl, map.Values.ToArray());
                    if (mat.Count() == g.Count())
                        return mat;
                }
            }
        }
    }

    return Group.Generate(new GL(1,2));
}

#endregion

void MatrixFormTinyGroups(int maxOrder)
{
    var total = 0;
    var missing = new List<(WordGroup g, SubGroupsInfos infos, ANameElt[] names)>();
    foreach (var g in maxOrder.Range(1).SelectMany(o => FG.AllGroupsOfOrder(o)))
    {
        ++total;
        var gSubgrs = g.AllSubgroups().ToGroupWrapper();
        var names = NamesTree.BuildName(gSubgrs);
        var (mat, check) = MatrixFormFromNames(names[0]);

        if (mat.Count() == 1 && !check)
        {
            var mat0 = MatrixFormFromNamesMeth2(g, names);
            if (mat0.Count() == 1)
            {
                missing.Add((g, gSubgrs.Infos, names));
                continue;
            }
            else
                mat = mat0;
        }

        FG.DisplayName(mat, mat.AllSubgroups(), names, false, false, 20);

        if (!mat.IsIsomorphicTo(g))
            throw new();
    }

    Console.WriteLine($"Missing:{missing.Count} Found:{total - missing.Count}/{total}");
    foreach (var (g, infos, names) in missing)
    {
        names.Where(e => e is SemiDirectProductOp e0 &&
            e0.Lhs.ContentGroup!.GroupType == GroupType.AbelianGroup &&
            e0.Rhs.ContentGroup!.GroupType == GroupType.AbelianGroup)
        .Println(g.ShortName);
    }
}

{
    Group.ActivedStorage(false);
    GlobalStopWatch.Restart();
    MatrixFormTinyGroups(32);
    GlobalStopWatch.Show("End");
    Console.Beep();
}

/*
    Missing:11 Found:133/144
    |D16 x: C2| = 32
        C8 x: (C2 x C2)
    |QD16 x: C2| = 32

    |C4 x: Q8| = 32

    |Q8 x: C4| = 32

    |D8 x: (C2 x C2)| = 32
        (C4 x C2) x: (C2 x C2)
        (C2 x C2 x C2) x: (C2 x C2)
    |MM16 x: C2| = 32

    |M(4x:4)3 x: C2| = 32

    |(C2 x Q8) x: C2| = 32

    |C8 . C4| = 32

    |C4 . D8| = 32

    |(C4 x C4) . C2| = 32

    # End Time:1m22s
*/

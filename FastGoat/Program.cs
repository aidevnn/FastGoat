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

ConcreteGroup<Mat> SpecialCases(WordGroup g)
{
    // not implemented yet
    return Group.Generate(new GL(1, 2));
}

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
            var mat0 = SpecialCases(g);
            if (mat0.Count() == 1)
            {
                missing.Add((g, gSubgrs.Infos, names));
                continue;
            }
        }

        FG.DisplayName(gSubgrs.Parent, gSubgrs, names, false, false, 20);
        DisplayGroup.Generators(mat);
        Console.WriteLine();

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

    missing.Where(e => e.names[0].ContentType == ANameElt.NodeType.DirectProduct).Println(e => e.names[0].ContentGroup!.ShortName, "Missing Direct Product");
    missing.SelectMany(e => e.names.Where(f => f is ExtensionOp f0 && f0.Lhs.ContentGroup!.GetGenerators().Count() == 1 && f0.Rhs.ContentGroup!.GetGenerators().Count() == 1).Take(1))
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
    var c2c4d8 = ProductMatrixBlock(c2c4, d8);
    DisplayGroup.HeadNames(c2c4d8);
    DisplayGroup.Generators(c2c4d8);

    var c3c3 = FG.AbelianMat(3, 3);
    var dic3 = FG.DicyclicGL2p(3);
    var c3c3dic3 = ProductMatrixBlock(c3c3, dic3);
    DisplayGroup.HeadNames(c3c3dic3);
    DisplayGroup.Generators(c3c3dic3);
}

{
    Group.ActivedStorage(false);
    // MatrixFormTinyGroups(24);
    // MtCyclicSubgroups(16, factor: 8).Select(e => e.ToGroupWrapper()).FilterIsomorphic().Naming().DisplayNames();
    // MtCyclicSubgroups(32, factor: 5).Select(e => e.ToGroupWrapper()).FilterIsomorphic().Naming().DisplayNames();
    // MtCyclicSubgroups(48, factor: 5).Select(e => e.ToGroupWrapper()).FilterIsomorphic().Naming().DisplayNames();
    // MtCyclicSubgroups(64, factor: 5).Select(e => e.ToGroupWrapper()).FilterIsomorphic().Naming().DisplayNames();
}

void ProductGroupSubgroups<T>(params ConcreteGroup<T>[] tp) where T : struct, IElt<T>
{
    var g = Product.GpGenerate(tp);
    var gSubgrs = g.AllSubgroups();
    gSubgrs.Where(cj => cj.GroupType == GroupType.NonAbelianGroup)
        .Select(cj => gSubgrs.Restriction(cj.Representative))
        .FilterIsomorphic()
        .DisplayBoxes();
}

void ManualSearch()
{
    Group.ActivedStorage(false);
    // Ring.MatrixDisplayForm = Ring.MatrixDisplay.OneLineArray;
    
    // |(C2 x C2) x: C4| = 16   Gap SmallGroup(16,3)
    // ProductGroupSubgroups(FG.DihedralGL2p(4), FG.AbelianMat(4));
    var c4d8 = ProductMatrixBlock(FG.DihedralGL2p(4), FG.AbelianMat(4));
    var g_16_3 = FG.WordGroup("(C2 x C2) x: C4", "a4, b2, c2, bcbc, caca-1, abca-1b");
    // DisplayGroup.HeadOrders(c4d8);
    // DisplayGroup.HeadOrders(g_16_3);
    DisplayGroup.Generators(Group.IsomorphicSubgroup(c4d8, g_16_3));

    // |D8 x: C2| = 16          Gap SmallGroup(16,13)
    // ProductGroupSubgroups(FG.GL2p(5));
    var gl25 = FG.GL2p(5);
    var g_16_13 = FG.WordGroup("D8 x: C2", "a4, b2, c2, a2bcbc, baba-1, caca-1");
    // DisplayGroup.HeadOrders(gl25);
    // DisplayGroup.HeadOrders(g_16_13);
    DisplayGroup.Generators(Group.IsomorphicSubgroup(gl25, g_16_13));

    // |(C3 x C3) x: C2| = 18   Gap mallGroup(18,4)
    // ProductGroupSubgroups(FG.DihedralGL2p(3), FG.DihedralGL2p(3));
    var s3s3 = ProductMatrixBlock(FG.DihedralGL2p(3), FG.DihedralGL2p(3));
    var g_18_4 = FG.WordGroup("(C3 x C3) x: C2", "a3, b3, c2, acac, bcbc, bab-1a-1");
    // DisplayGroup.HeadOrders(s3s3);
    // DisplayGroup.HeadOrders(g_18_4);
    DisplayGroup.Generators(Group.IsomorphicSubgroup(s3s3, g_18_4));

    // |C3 x: D8| = 24          Gap SmallGroup(24,8)
    // ProductGroupSubgroups(FG.DihedralGL2p(3), FG.DihedralGL2p(4));
    var s3d8 = ProductMatrixBlock(FG.DihedralGL2p(3), FG.DihedralGL2p(4));
    var g_24_8 = FG.WordGroup("C3 x: D8", "a6, b2, c2, abab, a3bcbc, caca-1");
    // DisplayGroup.HeadOrders(s3d8);
    // DisplayGroup.HeadOrders(g_24_8);
    DisplayGroup.Generators(Group.IsomorphicSubgroup(s3d8, g_24_8));

    // |(C3 x C3) x: C3| = 27   Gap SmallGroup(27,3)
    // Character Table with max dimension 3, in GL(3,Q(ξ3)), 
    // lucky candidat in GL(3,3)
    var gl33 = FG.GLnp(3, 3);
    var g_27_3 = FG.WordGroup("(C3 x C3) x: C3", "a3, b3, c3, bab-1a-1, cbc-1b-1, ab-1ca-1c-1");
    // DisplayGroup.HeadOrders(gl33);
    // DisplayGroup.HeadOrders(g_27_3);
    DisplayGroup.Generators(Group.IsomorphicSubgroup(gl33, g_27_3));
}
/*
    Generators of (C2 x C2) x: C4 in GL(3,5)
    gen1 of order 2
    [0, 3, 0]
    [2, 0, 0]
    [0, 0, 1]
    gen2 of order 4
    [0, 1, 0]
    [1, 0, 0]
    [0, 0, 2]

    Generators of D8 x: C2 in GL(2,5)
    gen1 of order 2
    [0, 3]
    [2, 0]
    gen2 of order 2
    [4, 0]
    [0, 1]
    gen3 of order 4
    [3, 0]
    [0, 3]

    Generators of (C3 x C3) x: C2 in GL(4,7)
    gen1 of order 2
    [0, 1, 0, 0]
    [1, 0, 0, 0]
    [0, 0, 0, 1]
    [0, 0, 1, 0]
    gen2 of order 3
    [1, 0, 0, 0]
    [0, 1, 0, 0]
    [0, 0, 4, 0]
    [0, 0, 0, 2]
    gen3 of order 3
    [4, 0, 0, 0]
    [0, 2, 0, 0]
    [0, 0, 1, 0]
    [0, 0, 0, 1]

    Generators of C3 x: D8 in GL(4,13)
    gen1 of order 2
    [ 0,  1,  0,  0]
    [ 1,  0,  0,  0]
    [ 0,  0,  0,  1]
    [ 0,  0,  1,  0]
    gen2 of order 2
    [ 1,  0,  0,  0]
    [ 0,  1,  0,  0]
    [ 0,  0,  0,  5]
    [ 0,  0,  8,  0]
    gen3 of order 6
    [ 9,  0,  0,  0]
    [ 0,  3,  0,  0]
    [ 0,  0, 12,  0]
    [ 0,  0,  0, 12]

    Generators of (C3 x C3) x: C3 in GL(3,3)
    gen1 of order 3
    [0, 2, 2]
    [1, 1, 0]
    [0, 1, 2]
    gen2 of order 3
    [1, 0, 0]
    [1, 1, 0]
    [2, 0, 1]
*/

IEnumerable<Mat[]> MGenerators(int p, int[] type, int dim)
{
    if (dim < type.Length)
        throw new();
        
    var sz = type.Length;
    var m = Gcd(type);
    var Up = FG.UnInt(p);
    var sn = FG.Symmetric(dim);
    var gl = new GL(dim, p);
                
    var primaries = type.Select(k => Up.Where(e => Up.ElementsOrders[e] == k).ToArray())
                        .MultiLoop()
                        .Select(l => l.Select(z => z.K).ToArray())
                        .ToArray();
                       
    var secondaries = Up.Where(e => m % Up.ElementsOrders[e] == 0)
                .MultiLoop(dim - type.Length)
                .Select(l => l.Select(z => z.K).ToArray())
                .ToArray();
                
    if (secondaries.Length == 0)
        secondaries = new int[1][] { new int[0] };
    
    foreach(var (f, s) in primaries.Grid2D(secondaries))
    {
        foreach(var o in sn)
        {
            foreach(var l in f.Select((k, i) => new[]
                    {
                        o.Apply(sz.Range().Select(j => j == i ? k : 1).Concat(s).ToArray()),
                        o.Apply(Enumerable.Repeat(k, sz).Concat(s).ToArray())
                    }).MultiLoop())
                yield return l.Select(v => gl.Create(MatrixExt.Diagonal(v.ToArray()))).ToArray();
        }
    }
}

IEnumerable<(int[] perm, int[][] cycles, Mat mat)> CGenerators(int m, int n, int dim)
{
    var distinctTypes = IntExt.Partitions32[dim].Where(l => l.Count == l.Distinct().Count()).Select(l => l.Order().ToArray())
            .OrderBy(l => l.Length).ToArray();
    var nks = distinctTypes.Select(l => l.Aggregate((a0, a1) => a0 * a1))
        .SelectMany(e => IntExt.Dividors(e).Append(e).Where(j => j != 1)).Append(n).ToHashSet();
    foreach (var p in nks.SelectMany(nk => IntExt.Primes10000.Where(p => (p - 1) % m == 0 && (p - 1) % nk == 0).Take(6)).Distinct().Order().Take(6))
    {
        var Up = FG.UnInt(p);
        var gl = new GL(dim, p);
        var matn = Up.MultiLoop(dim).Select(l => gl.Create(MatrixExt.Diagonal(l.Select(e => e.K).ToArray()))).ToArray();
        Console.WriteLine($"{gl} press key...");
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

void PrepareWordGroup(string rels, int dim)
{
    var g = FG.WordGroup("H", rels);
    var gSubgrs = g.AllSubgroups().ToGroupWrapper();
    var names = NamesTree.BuildName(gSubgrs);
    g.Name = gSubgrs.Parent.Name = names[0].Name;
    var (Mgens, Cgen, name) = ExtractGeneratorsSdp(names);
    Mgens.Append(Cgen).Select(w => (w, g.ElementsOrders[w])).Println(name);
    AbelianByCyclicSdp(g, Mgens, Cgen, dim);
}

(Word[] Mgens, Word Cgen, string name) ExtractGeneratorsSdp(ANameElt[] names)
{
    foreach(SemiDirectProductOp e0 in names.Where(e => e is SemiDirectProductOp e0 &&
            e0.Lhs.ContentGroup!.GroupType == GroupType.AbelianGroup &&
            e0.Rhs.ContentGroup!.GroupType == GroupType.AbelianGroup))
    {
        var Mgens = e0.Lhs.ContentGroup!.GetGenerators().Select(m0 => (Word)m0.E)
            .Select(e => char.IsUpper(e.Get()[0]) ? new Word(e.WGroup, e.Get().Revert()) : e).ToArray();
        var Cgens = e0.Rhs.ContentGroup!.GetGenerators().Select(m0 => (Word)m0.E)
            .Select(e => char.IsUpper(e.Get()[0]) ? new Word(e.WGroup, e.Get().Revert()) : e).ToArray();
        
        if (Cgens.Length == 1 && Cgens.Concat(Mgens).All(w => w.Get().All(c => char.IsLower(c))))
            return (Mgens, Cgens[0], e0.Name);
    }
    
    throw new();
}

void AbelianByCyclicSdp(WordGroup g, Word[] Mgens, Word Cgen, int dim)
{
    var type = Mgens.Select(e => g.ElementsOrders[e]).ToArray();
    var m = Gcd(type);
    var n = g.ElementsOrders[Cgen];
    // foreach (var dim in 4.Range(1).Where(d => d != 1 && (IntExt.Gcd(m, d) != 1 || IntExt.Gcd(m - 1, d) != 1)))
    {
        foreach(var (perm, cycles, m1) in CGenerators(m, n, dim))
        {
            var gl = m1.GL;
            var p = gl.P;
            Console.WriteLine($"{g} in {gl} type:[{type.Glue(" ")}] cycles:{cycles.Select(c => c.Glue(" ")).Glue("","({0})")} m:{m} n:{n}");
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
                    {
                        Console.WriteLine();
                        Console.WriteLine("############   Found   ############");
                        DisplayGroup.HeadGenerators(mat);
                        Console.WriteLine("############   Found   ############");
                        Console.WriteLine();
                        return;
                    }
                }
            }
        }
    }

    throw new GroupException(GroupExceptionType.GroupDef);
}

{   
    Group.ActivedStorage(false);
    Ring.MatrixDisplayForm = Ring.MatrixDisplay.OneLineArray;
    
    PrepareWordGroup("a4, b2, c2, bcbc, caca-1, abca-1b", 3); // |(C2 x C2) x: C4| = 16 Gap SmallGroup(16,3)
    PrepareWordGroup("a4, b2, c2, a2bcbc, baba-1, caca-1", 2); // |D8 x: C2| = 16 Gap SmallGroup(16,13)
    PrepareWordGroup("a3, b3, c3, bab-1a-1, cbc-1b-1, ab-1ca-1c-1", 3); // |(C3 x C3) x: C3| = 27 Gap SmallGroup(27,3)
    PrepareWordGroup("a4, b4, c2, caca-1, bab-1a-1, ab-1cb-1c", 2); // |(D8 x: C4)| = 32 Gap SmallGroup(32,11)
    PrepareWordGroup("a4, b2, c2, d2, bcbc, bdbd, cdcd, baba-1, dada-1, acda-1c", 4); // |(C2 x C2 x C2) x: C4| = 32 Gap SmallGroup(32,22)
    PrepareWordGroup("a8, b2, c2, bcbc, caca-1, abca-1b", 3); // |(C2 x C2) x: C8| = 32 Gap SmallGroup(32,5)
    PrepareWordGroup("a6, b6, c2, caca-1, bab-1a-1, ab-1cb-1c", 4); // |(D12 x: C6)| = 72 Gap SmallGroup(72,30) super group for C3 x: D8
}


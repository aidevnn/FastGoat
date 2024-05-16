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

int[][] ChangeGL(ConcreteGroup<Mat> m, int p)
{
    var gl0 = m.Neutral().GL;
    var p0 = gl0.P;
    var Up0 = FG.UnInt(p0);
    var Up = FG.UnInt(p);
    var iso0 = Group.AllMorphisms(Up0, Up, Group.MorphismType.Isomorphism)
        .First()
        .HomMap.ToDictionary(e => e.Key.K, e => e.Value.K);
    iso0[0] = 0;
    return m.GetGenerators().Select(e => e.Table.Select(k => iso0[k]).ToArray()).ToArray();
}

ConcreteGroup<Mat> ProductMatrixBlock(ConcreteGroup<Mat> group0, ConcreteGroup<Mat> group1)
{
    var (gl0, gl1) = (group0.Neutral().GL, group1.Neutral().GL);
    var dim = gl0.N + gl1.N;
    var (p0, p1) = (gl0.P, gl1.P);
    var p = Primes10000.First(p => (p - 1) % (p0 - 1) == 0 && (p - 1) % (p1 - 1) == 0);

    var gl = new GL(dim, p);
    var gens0 = ChangeGL(group0, p).Select(e => (e, gl0.N)).ToArray();
    var gens1 = ChangeGL(group1, p).Select(e => (e, gl1.N)).ToArray();

    var Gens0 = gens0.Select(e => MatrixExt.MergeDiagonalBlocks(e, (gl1.Neutral().Table, gl1.N)))
        .Select(e => gl.Create(e)).ToArray();
    var Gens1 = gens1.Select(e => MatrixExt.MergeDiagonalBlocks((gl0.Neutral().Table, gl0.N), e))
        .Select(e => gl.Create(e)).ToArray();
    return Group.Generate($"{group0.NameParenthesis()} x {group1.NameParenthesis()}", gl,
        Gens0.Concat(Gens1).ToArray());
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
    if ((dim < 4 && BigInteger.Pow(p, dim) > 10000) || BigInteger.Pow(p, dim) > 150000)
        yield break;

    var sz = type.Length;
    var m = Gcd(type);
    var Up = FG.UnInt(p);
    var gl = new GL(dim, p);

    var byOrders = Up.MultiLoop(dim).Select(l => gl.Create(MatrixExt.Diagonal(l.Select(z => z.K).ToArray())))
        .GroupBy(mat => Group.Cycle(gl, mat).Count)
        .ToDictionary(e => e.Key, e => e.ToArray());

    if (type.Any(k => !byOrders.ContainsKey(k)))
        yield break;

    foreach (var mat in type.Select(f => byOrders[f]).MultiLoop().Select(ms => ms.ToArray()))
        yield return mat;
}

IEnumerable<(int[] perm, int[][] cycles, Mat mat)> CGenerators(int m, int n, int dim)
{
    var allTypes = IntExt.Partitions32[dim].Select(l => l.Order().ToArray()).OrderBy(l => l.Length).ToArray();
    var nks = allTypes.Select(l => l.Aggregate((a0, a1) => a0 * a1))
        .SelectMany(e => IntExt.Dividors(e).Append(e).Where(j => j != 1)).Append(n).ToHashSet();
    // nks.Println($"nks dim:{dim}");
    foreach (var p in nks
                 .SelectMany(nk => IntExt.Primes10000.Where(p => (p - 1) % m == 0 && (p - 1) % nk == 0).Take(10))
                 .Distinct().Order().Where(p => p < 62))
    {
        foreach (var res in CGeneratorsP(p, n, dim))
            yield return res;
    }
}

IEnumerable<(int[] perm, int[][] cycles, Mat mat)> CGeneratorsP(int p, int n, int dim)
{
    var Up = FG.UnInt(p);
    var gl = new GL(dim, p);
    var matn = Up.Where(e => n % Up.ElementsOrders[e] == 0)
        .OrderBy(e => Up.ElementsOrders[e])
        .Select(e => gl.At(gl.Neutral().Table, 0, e.K))
        .ToArray();

    // Console.WriteLine($"{gl} press key...");
    // Console.ReadLine();
    var sn = new Sn(dim);
    var m1s = IntExt.Partitions32[dim] //.Where(l => l.Count == l.Distinct().Count())
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

    foreach (var e in m1s)
        yield return e;
}

(Word[] Mgens, Word[] Cgen, string name) ExtractGenerators(ANameElt[] names)
{
    foreach (var e0 in names.Where(e => e is SemiDirectProductOp e0 &&
                                        e0.Lhs.ContentGroup!.GroupType == GroupType.AbelianGroup &&
                                        e0.Rhs.ContentGroup!.GroupType == GroupType.AbelianGroup)
                 .Cast<SemiDirectProductOp>())
    {
        var Mgens = e0.Lhs.ContentGroup!.GetGenerators().Select(m0 => (Word)m0.E).ToArray();
        var Cgens = e0.Rhs.ContentGroup!.GetGenerators().Select(m0 => (Word)m0.E).ToArray();

        if (Cgens.Length == 1 && Mgens.Length != 0)
            return (Mgens, Cgens, e0.Name);
    }

    return (new Word[0], new Word[0], "");
}

ConcreteGroup<Mat> MatrixFormFromNamesMeth2(WordGroup g, ANameElt[] names)
{
    var (mgens, cgens, name) = ExtractGenerators(names);
    if (cgens.Length == 0)
        return Group.Generate(new GL(1, 2));

    var lvl = Logger.Level;
    Logger.Level = LogLevel.Off;
    var wg = FG.WordGroup(name, Graph.DefiningRelatorsOfGroup(g, mgens.Concat(cgens).ToArray()));
    Logger.Level = lvl;
    var Mgens = wg.GetGenerators().SkipLast(1).ToArray();
    var Cgen = wg.GetGenerators().Last();

    var mtype = Mgens.Select(e => wg.ElementsOrders[e]).ToArray();
    var m = Gcd(mtype);
    var c = wg.ElementsOrders[Cgen];
    foreach (var dim in 7.Range(1).Where(d => d != 5 && d >= mtype.Length))
    {
        foreach (var (perm, cycles, m1) in CGenerators(m, c, dim))
        {
            var gl = m1.GL;
            var p = gl.P;
            // Console.WriteLine($"{g} in {gl} type:[{type.Glue(" ")}] cycles:{cycles.Select(c => c.Glue(" ")).Glue("","({0})")} m:{m} n:{n}");
            // Console.WriteLine($"m1:{m1}");
            foreach (var m0s in MGenerators(p, mtype, dim))
            {
                var map = Mgens.Zip(m0s).ToDictionary(e => e.First.Get()[0], e => e.Second);
                map[Cgen.Get()[0]] = m1;
                if (wg.CheckHomomorphism(gl, map))
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

    return Group.Generate(new GL(1, 2));
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
        var mat0 = MatrixFormFromNamesMeth2(g, names);
        if (mat0.Count() == 1)
        {
            (mat0, bool check) = MatrixFormFromNames(names[0]);
            if (mat0.Count() == 1 && !check)
            {
                missing.Add((g, gSubgrs.Infos, names));
                continue;
            }
        }

        mat0.Name = g.Name;
        FG.DisplayName(mat0, mat0.AllSubgroups(), names, false, false, 20);

        if (!mat0.IsIsomorphicTo(g))
            throw new();
    }

    Console.WriteLine($"Missing:{missing.Count} Found:{total - missing.Count}/{total}");
    foreach (var (g, infos, names) in missing)
    {
        names.Where(e => e is SemiDirectProductOp e0 &&
                         (e0.Lhs.ContentGroup!.GroupType == GroupType.AbelianGroup ||
                          e0.Rhs.ContentGroup!.GroupType == GroupType.AbelianGroup))
            .Println(g.ShortName);
    }
}

void CheckGroup(WordGroup g)
{
    var gSubgrs = g.AllSubgroups().ToGroupWrapper();
    var names = NamesTree.BuildName(gSubgrs);
    var mat0 = MatrixFormFromNamesMeth2(g, names);
    FG.DisplayName(mat0, mat0.AllSubgroups(), names, false, false, 20);

    if (!mat0.IsIsomorphicTo(g))
        throw new();
}

Mat GLDiagOrdGenerators(GL gl, int ord)
{
    var id = gl.Neutral().Table;
    var e0 = Solve_k_pow_m_equal_one_mod_n_strict(gl.P, ord);
    return gl.At(id, 0, e0);
}

Mat[] GLPermGenerators(GL gl)
{
    var id = gl.Neutral().Table;
    var idc = id.Chunk(gl.N).ToArray();
    var sn = new Sn(gl.N);
    var genSn = sn.GetGenerators().ToArray();
    var M2 = gl.Create(genSn[0].Apply(idc).SelectMany(l => l).ToArray());
    var Mn = gl.Create(genSn[1].Apply(idc).SelectMany(l => l).ToArray());
    return new Mat[2] { M2, Mn };
}

Mat[] GLDiagPermOrdGenerators(int n, int p, int ord)
{
    var gl = new GL(n, p);
    var diag = GLDiagOrdGenerators(gl, ord);
    var perms = GLPermGenerators(gl);
    return perms.Append(diag).ToArray();
}

ConcreteGroup<Mat> GLDiagPermOrd(int n, int p, int ord)
{
    var gens = GLDiagPermOrdGenerators(n, p, ord);
    return Group.Generate($"DPGL({n},{p})", gens[0].GL, gens);
}

ConcreteGroup<Mat> GLDiagPerm(int n, int p) => GLDiagPermOrd(n, p, p - 1);

ConcreteGroup<Mat> SearchDiagPermGL<T>(ConcreteGroup<Mat> dpgl, WordGroup g, ConcreteGroup<T> m1, ConcreteGroup<T> m2, bool show = true)
    where T : struct, IElt<T>
{
    var allIso2 = Group.AllMorphisms(m2, dpgl, Group.MorphismType.Isomorphism).ToArray();
    var k1 = 0;
    foreach (var iso1 in Group.AllMorphisms(m1, dpgl, Group.MorphismType.Isomorphism))
    {
        ++k1;
        var k2 = 0;
        var gens1 = m1.GetGenerators().Select(e => iso1[e]).ToArray();
        foreach (var iso2 in allIso2)
        {
            // Console.WriteLine($"Search Iso1:{k1} Iso2:{++k2}/{allIso2.Length}");
            var gens = m2.GetGenerators().Select(e => iso2[e]).Concat(gens1).Distinct().ToArray();
            var matForm = Group.Generate(g.Name, dpgl, gens);
            if (matForm.IsIsomorphicTo(g))
            {
                if (show)
                {
                    var matFormAllSubgroups = matForm.AllSubgroups();
                    var names = NamesTree.BuildName(matFormAllSubgroups.ToGroupWrapper(), renaming: false);
                    FG.DisplayName(matForm, matFormAllSubgroups, names, false, false, 20);
                }
                
                return matForm;
            }
        }
    }

    return Group.Generate(new GL(1, 2));
}

void MissingGroupsOrder32()
{
    Group.ActivedStorage(false);
    GlobalStopWatch.AddLap();

    var dpgl25 = GLDiagPerm(2, 5);
    DisplayGroup.HeadOrders(dpgl25);
    var dpgl217 = GLDiagPerm(2, 17);
    DisplayGroup.HeadOrders(dpgl217);
    var dpgl35 = GLDiagPerm(3, 5);
    DisplayGroup.HeadOrders(dpgl35);
    var dpgl317 = GLDiagPerm(3, 17);
    DisplayGroup.HeadOrders(dpgl317);
    var dpgl43 = GLDiagPerm(4, 3);
    DisplayGroup.HeadOrders(dpgl43);
    var dpgl45 = GLDiagPerm(4, 5);
    DisplayGroup.HeadOrders(dpgl45);

    var lvl = Logger.Level;
    Logger.Level = LogLevel.Off;
    var g1 = FG.WordGroup("D16 x: C2", "a2, b2, c2, acac, abcbabcb, bcbcbcbc, abababcabc");
    var g2 = FG.WordGroup("QD16 x: C2", "a4, b2, c2, a2ba2b, a2bcbc, caca-1, abababa-1b");
    var g3 = FG.WordGroup("Q8 x: C4", "b4, c2, a4b-2, ababc, caca-1, cbcb-1, ab-1acb-1");
    var g4 = FG.WordGroup("M(4x:4)3 x: C2", "a4, d2, a2b2, a2c2, acac, bcbcd, abab-1, dada-1, dbdb-1, dcdc-1");
    var g5 = FG.WordGroup("C4 x: Q8", "a4, b4, a2c2, acac-1, baba-1, cbc-1b-1");
    var g6 = FG.WordGroup("MM16 x: C2", "a8, b2, a2ba2b, aba-1baba-1b");
    var g7 = FG.WordGroup("C8 . C4", "c2, ab2ac, a2cb-2, cbcb-1, abcab-1");
    var g8 = FG.WordGroup("D8 x: (C2 x C2)", "a4, b2, c2, d2, abab, bcbc, bdbd, a2cdcd, caca-1, dada-1");
    var g9 = FG.WordGroup("(C2 x Q8) x: C2", "a4, c2, d2, a2b2, a2cdcd, abab-1, caca-1, cbcb-1, dada-1, dbdb-1");
    var g10 = FG.WordGroup("(C4 x C4) . C2", "a4, b4, a2c2, baba-1, ab2cac-1, cbc-1b-1");
    var g11 = FG.WordGroup("C4 . D8", "b4, c2, a4b-2, cbcb-1, baba-1c, abca-1b-1");
    
    // (C4 x C4) . C2 is subgroup of (C8 x C4) : C2
    var gbs0 = FG.WordGroup("(C8 x C4) : C2", "a4, c2, a2b2, a2ca2c, abacbc, abababa-1b-1, abcacaca-1cb-1");
    // C4 . D8 is a subgroup of (C2 x QD16) x: C2
    var gbs1 = FG.WordGroup("(C2 x QD16) x: C2", "b2, c2, d2, cdcd, a3dad, a2cbcb, a3ba-1b, acda-1c");
    Logger.Level = lvl;
    
    var (c2, c4, c8) = (FG.AbelianMat(2), FG.AbelianMat(4), FG.AbelianMat(8));
    
    SearchDiagPermGL(dpgl45, g1, FG.DihedralGL2p(8), c2);
    SearchDiagPermGL(dpgl45, g2, FG.SemiDihedralGL2p(4), c2);
    SearchDiagPermGL(dpgl317, g3, FG.Quaternion(8), c4);
    SearchDiagPermGL(dpgl45, g4, FG.MetaCyclicSdpMat(4, 4, 3), c2);
    SearchDiagPermGL(dpgl45, g5, FG.Quaternion(8), c4);
    SearchDiagPermGL(dpgl43, g6, FG.ModularMaxGL2p(4), c2);
    SearchDiagPermGL(dpgl217, g7, c8, c8);
    SearchDiagPermGL(dpgl43, g8, ProductMatrixBlock(FG.DihedralGL2p(4), c2), c2);
    
    var gl = new GL(2, 5);
    var d8byc2 = Group.Generate("D8 x: C2", gl, gl[1, 0, 0, 4], gl[0, 1, 1, 0], gl[2, 0, 0, 2]);
    SearchDiagPermGL(dpgl45, g9, d8byc2, c2); // slow ~ 40s

    var gSubgrs = gbs0.AllSubgroups();
    var names = NamesTree.BuildName(gSubgrs.ToGroupWrapper());
    var matBs = MatrixFormFromNamesMeth2(gbs0, names);
    SearchDiagPermGL(matBs, g10, FG.AbelianMat(4, 4), c4);

    var prod = ProductMatrixBlock(FG.SemiDihedralGL2p(4), FG.AbelianMat(2, 1));
    var gl417 = prod.Neutral().GL;
    var perms = GLPermGenerators(gl417);
    var s = Group.Generate("Sub-GL(4,17)", gl417, prod.GetGenerators().Concat(perms).ToArray());

    var mBs = SearchDiagPermGL(s, gbs1, prod, FG.AbelianMat(2), show: false);
    SearchDiagPermGL(mBs, g11, FG.ModularMaxGL2p(4), FG.AbelianMat(4));
    
    GlobalStopWatch.Show("END missing groups of order 32"); // Time:49.336s
    Console.WriteLine();
}

void MatrixFormUpto63()
{
    Group.ActivedStorage(false);
    GlobalStopWatch.AddLap();
    MatrixFormTinyGroups(63); // Missing:23 Found:296/319
    GlobalStopWatch.Show("END groups upto 63");
    Console.WriteLine();
}

void MissingGroupsOrder48and54()
{
    GlobalStopWatch.AddLap();
    
    var dpgl37 = GLDiagPerm(3, 7);
    var dpgl313 = GLDiagPerm(3, 13);
    var dpgl47 = GLDiagPerm(4, 7);
    
    var (c2, c4, c3) = (FG.AbelianMat(2), FG.AbelianMat(4), FG.AbelianMat(3));
    
    var lvl = Logger.Level;
    Logger.Level = LogLevel.Off;
    var g1 = FG.WordGroup("C3 x: D16", "a8, b2, a2ba2b, ababababa-1ba-1b");
    var g2 = FG.WordGroup("D24 x: C2", "b2, c2, abab, a2ca2c, a2cbcb, a5ca-1c");
    var g3 = FG.WordGroup("(C3 x: QD16)[1]", "c2, a3b-2, b2cb2c, a2ba-1b, ca2ca-2, bcb-1cb2");
    var g4 = FG.WordGroup("(C3 x: QD16)[2]", "c4, bcbc, a3b-2, caca-1, a2ba-1b, bc-1bc-1");
    var g5 = FG.WordGroup("A4 x: C4", "a4, b3, a2ba2b-1, abababab");
    var g6 = FG.WordGroup("C3 x: Q16", "b4, a4b-2, a2ba2b-1, abababab-1a-1b-1a-1b-1");
    var g7 = FG.WordGroup("Dic3 x: C4", "b4, c2, a6b-2, ababc, caca-1, cbcb-1, ab-1acb-1");
    var g8 = FG.WordGroup("Dic6 x: C2", "b4, c2, a2ca2c, ab2cac, abab-1, a2cbcb-1, a5b-1a-1b-1");
    var g9 = FG.WordGroup("SL(2,3) x: C2", "a4, c3, a2b2, abab, acacac, cbc-1b-1");
    var g10 = FG.WordGroup("C3 x: MM16", "c2, a3b-2, b2cbcb, caca-1, a2ba-1b");
    var g11 = FG.WordGroup("(C3 x C3) x: S3", "a3, b3, c3, d2, adad, cdcd, dbdb-1, bab-1a-1, cbc-1b-1, ab-1ca-1c-1");
    Logger.Level = lvl;
    
    var gl = new GL(2, 5);
    var d8byc2 = Group.Generate("D8 x: C2", gl, gl[1, 0, 0, 4], gl[0, 1, 1, 0], gl[2, 0, 0, 2]);
    
    SearchDiagPermGL(dpgl47, g1, FG.DihedralGL2p(12), c2);
    SearchDiagPermGL(ProductMatrixBlock(d8byc2, FG.DihedralGL2p(3)), g2, FG.DihedralGL2p(12), c2);
    SearchDiagPermGL(dpgl47, g3, FG.MetaCyclicSdpMat(3, 8, 2), c2);
    SearchDiagPermGL(ProductMatrixBlock(FG.DihedralGL2p(3), FG.SemiDihedralGL2p(4)), g4, FG.SemiDihedralGL2p(4), c3);
    
    var gl33 = new GL(3, 3);
    var a4 = Group.Generate("A4", gl33, gl33[0, 1, 0, 0, 0, 1, 1, 0, 0], gl33[0, 0, 1, 2, 0, 0, 0, 2, 0]);
    SearchDiagPermGL(dpgl313, g5, a4, c4);
    SearchDiagPermGL(ProductMatrixBlock(FG.DihedralGL2p(3), FG.Quaternion(16)), g6, FG.Quaternion(16), c3);
    SearchDiagPermGL(ProductMatrixBlock(FG.DicyclicGL2p(6), FG.MetaCyclicSdpMat(12, 2, 5)), g7, FG.DicyclicGL2p(3), c4);
    
    SearchDiagPermGL(ProductMatrixBlock(FG.DihedralGL2p(3), d8byc2), g8, FG.DicyclicGL2p(6), c2);
    SearchDiagPermGL(ProductMatrixBlock(FG.GL2p(3), FG.DihedralGL2p(4)), g9, FG.SL2p(3), c2);
    SearchDiagPermGL(ProductMatrixBlock(FG.DihedralGL2p(3), FG.ModularMaxGL2p(4)), g10, FG.ModularMaxGL2p(4), c3);
    
    var c3c3 = FG.AbelianMat(3, 3);
    var s3 = FG.DihedralGL2p(3);
    SearchDiagPermGL(dpgl37, g11, c3c3, s3);

    GlobalStopWatch.Show("END missing groups of order 48 and 54");
}

void Run()
{
    GlobalStopWatch.Restart();
    
    MatrixFormUpto63();
    
    MissingGroupsOrder32();
    MissingGroupsOrder48and54();
    
    Console.Beep();
    GlobalStopWatch.Show("END"); // Time:4m3s, Missing C2 . S4
}

{
    var c2s4 = FG.WordGroup("C2 . S4", "a4, c3, a2b2, abab-1, acacac, cbcb-1");

    DisplayGroup.HeadOrdersNames(c2s4);
    Console.WriteLine($"Generators of ({c2s4}):{c2s4.GetGenerators().ToDictionary(e => e, e => c2s4.ElementsOrders[e]).GlueMap()}");
    Console.WriteLine();
    // C2 . S4 ~ SL(2,3) . C2
    // generators of order 4, 3, and 4
    
    FG.CharacterTable(c2s4).DisplayCells();
    // Max chi dimension 4, and √2 in Q(ξ4)
    // lucky candidat GL(4,5)
    
    var dpgl45 = GLDiagPerm(4, 5);
    SearchDiagPermGL(dpgl45, c2s4, FG.SL2p(3), FG.AbelianMat(4));
}

/*
   |C2 . S4| = 48
   Type        NonAbelianGroup
   BaseGroup   WG[a,b,c]
   
   Elements Orders : [1]:1, [2]:1, [3]:8, [4]:18, [6]:8, [8]:12
   SubGroupsInfos { AllSubGr = 35, AllConjsCl = 13, AllNorms = 5 }
   Group names
       C2 . S4
       SL(2,3) . C2
       Q8 . S3
   
   Generators of (C2 . S4):a->4, b->4, c->3
   
   |C2 . S4| = 48
   Type        NonAbelianGroup
   BaseGroup   WG[a,b,c]
   
   [Class      1   2   3  4a  4b   6   8a   8b]
   [ Size      1   1   8   6  12   8    6    6]
   [                                          ]
   [  Ꭓ.1      1   1   1   1   1   1    1    1]
   [  Ꭓ.2      1   1   1   1  -1   1   -1   -1]
   [  Ꭓ.3      2   2  -1   2   0  -1    0    0]
   [  Ꭓ.4      2  -2  -1   0   0   1   √2  -√2]
   [  Ꭓ.5      2  -2  -1   0   0   1  -√2   √2]
   [  Ꭓ.6      3   3   0  -1   1   0   -1   -1]
   [  Ꭓ.7      3   3   0  -1  -1   0    1    1]
   [  Ꭓ.8      4  -4   1   0   0  -1    0    0]
   All i,                 Sum[g](Xi(g)Xi(g^−1)) = |G|      : True
   All i <> j,            Sum[g](Xi(g)Xj(g^−1)) =  0       : True
   All g, h in Cl(g),     Sum[r](Xr(g)Xr(h^−1)) = |Cl(g)|  : True
   All g, h not in Cl(g), Sum[r](Xr(g)Xr(h^−1)) =  0       : True
   
   ########################################################
   ################        C2 . S4         ################
   ########################################################
   |C2 . S4| = 48
   NotSimple, NonAbelianGroup, NotNilpotent, Solvable
   
   Elements Orders : [1]:1, [2]:1, [3]:8, [4]:18, [6]:8, [8]:12
   SubGroupsInfos { AllSubGr = 35, AllConjsCl = 13, AllNorms = 5 }
   Lower Serie
   SL(2,3) --> C2 . S4
   Upper Serie
   C1      --> C2
   Derived Serie
   C1      --> C2      --> Q8      --> SL(2,3) --> C2 . S4
   Zentrum  Z(G) = C2
   Frattini Φ(G) = C2
   Fitting  F(G) = Q8
   
   Word Group
   < a,b,c | a4, c3, a2b2, abcb-1c >
   
   Generators of C2 . S4 in GL(4,5)
   gen1 of order 3
   [0, 0, 1, 0]
   [0, 1, 0, 0]
   [0, 0, 0, 1]
   [1, 0, 0, 0]
   gen2 of order 4
   [0, 0, 4, 0]
   [0, 0, 0, 1]
   [1, 0, 0, 0]
   [0, 4, 0, 0]
   gen3 of order 4
   [0, 3, 0, 0]
   [3, 0, 0, 0]
   [0, 0, 2, 0]
   [0, 0, 0, 3]
   
   
   Gap SmallGroup(48,28)    Name:C2 . S4 = SL(2,3) . C2
   Group names
       C2 . S4
       SL(2,3) . C2
       Q8 . S3
*/
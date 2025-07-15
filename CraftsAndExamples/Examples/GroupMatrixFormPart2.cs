using System.Numerics;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.Naming;
using FastGoat.Structures.Subgroups;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Characters;
using FastGoat.UserGroup.Matrix;
using FastGoat.UserGroup.Perms;
using FastGoat.UserGroup.Polynoms;
using FastGoat.UserGroup.Words;
using FastGoat.UserGroup.Words.Tools;

namespace CraftsAndExamples.Examples;

public static class GroupMatrixFormPart2
{
    static GroupMatrixFormPart2()
    {
        DPGLs = new();
        DPSLs = new();
    }

    #region GL methods

    static int[][] ChangeGL(ConcreteGroup<Mat> m, int p)
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

    static bool IsDiagPerm(Mat m)
    {
        var n = m.GL.N;
        var rg = n.Range();
        return rg.All(i => rg.Count(j => m.Table[i * n + j] == 0) == n - 1)
               && rg.All(j => rg.Count(i => m.Table[i * n + j] == 0) == n - 1);
    }

    static ConcreteGroup<Mat> ProductMatrixBlock(ConcreteGroup<Mat> group0, ConcreteGroup<Mat> group1)
    {
        var (gl0, gl1) = (group0.Neutral().GL, group1.Neutral().GL);
        var dim = gl0.N + gl1.N;
        var (p0, p1) = (gl0.P, gl1.P);
        var p = IntExt.Primes10000.First(p => (p - 1) % (p0 - 1) == 0 && (p - 1) % (p1 - 1) == 0);

        var gl = new GL(dim, p);
        var gens0 = ChangeGL(group0, p).Select(e => (e, gl0.N)).ToArray();
        var gens1 = ChangeGL(group1, p).Select(e => (e, gl1.N)).ToArray();

        var Gens0 = gens0.Select(e => MatrixExt.MergeDiagonalBlocks(e, (gl1.Neutral().Table, gl1.N)))
            .Select(e => gl.Create(e)).ToArray();
        var Gens1 = gens1.Select(e => MatrixExt.MergeDiagonalBlocks((gl0.Neutral().Table, gl0.N), e))
            .Select(e => gl.Create(e)).ToArray();

        if ((group0.GetGenerators().Any(mat => !IsDiagPerm(mat)) && p0 != p) ||
            (group1.GetGenerators().Any(mat => !IsDiagPerm(mat)) && p1 != p))
            throw new();

        var group2 = Group.Generate($"{group0.NameParenthesis()} x {group1.NameParenthesis()}", gl,
            Gens0.Concat(Gens1).ToArray());

        if (group2.Count() != group0.Count() * group1.Count())
            throw new();

        return group2;
    }

    static Dictionary<(int n, int p), ConcreteGroup<Mat>> DPGLs { get; }

    static Dictionary<(int n, int p), ConcreteGroup<Mat>> DPSLs { get; }

    public static ConcreteGroup<Mat> GetDPGL(int n, int p)
    {
        if (DPGLs.ContainsKey((n, p)))
            return DPGLs[(n, p)];

        var dpgl = DPGLs[(n, p)] = FG.DPGLnp(n, p);
        return dpgl;
    }
    
    public static ConcreteGroup<Mat> GetDPSL(int n, int p)
    {
        if (DPSLs.ContainsKey((n, p)))
            return DPSLs[(n, p)];

        var dpgl = DPSLs[(n, p)] = FG.DPSLnp(n, p);
        return dpgl;
    }

    #endregion

    #region Matrix Form from Group name

    static (ConcreteGroup<Mat> mat, bool isDiagByPerm) MatrixFormFromNames(ANameElt name)
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
                var n = coefs[0];
                var gens = pr == "S" ? FG.SnGensMat(n) : FG.AnGensMat(n);
                var gl = gens[0].GL;
                return (Group.Generate(name.Name, gl, gens), true);
            }
            else if (pr == "GL" || pr == "SL")
            {
                var gMat = pr == "GL" ? FG.GLnp(coefs[0], coefs[1]) : FG.SLnp(coefs[0], coefs[1]);
                if (coefs[0] == 2 && coefs[1] == 3)
                {
                    var lvl = Logger.SetOff();
                    var wg = FG.WordGroup(gMat.Name, Graph.DefiningRelatorsOfGroup(gMat));
                    Logger.Level = lvl;
                    var gMat0 = SearchDiagPermGL(FG.DPSLnp(4, 5), wg, gMat, FG.AbelianMat(1));
                    return (gMat0, true);
                }

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

    static IEnumerable<Mat[]> MGenerators(int p, int[] type, int dim)
    {
        if ((dim < 4 && BigInteger.Pow(p, dim) > 10000) || BigInteger.Pow(p, dim) > 150000)
            yield break;

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

    static IEnumerable<(int[] perm, int[][] cycles, Mat mat)> CGenerators(int m, int n, int dim)
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

    static IEnumerable<(int[] perm, int[][] cycles, Mat mat)> CGeneratorsP(int p, int n, int dim)
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

    static (Word[] Mgens, Word[] Cgen, string name) ExtractGenerators(ANameElt[] names)
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

    static ConcreteGroup<Mat> MatrixFormFromNamesMeth2(WordGroup g, ANameElt[] names)
    {
        var (mgens, cgens, name) = ExtractGenerators(names);
        if (cgens.Length == 0)
            return Group.Generate(new GL(1, 2));

        var lvl = Logger.SetOff();
        var wg = FG.WordGroup(name, Graph.DefiningRelatorsOfGroup(g, mgens.Concat(cgens).ToArray()));
        Logger.Level = lvl;
        var Mgens = wg.GetGenerators().SkipLast(1).ToArray();
        var Cgen = wg.GetGenerators().Last();

        var mtype = Mgens.Select(e => wg.ElementsOrders[e]).ToArray();
        var m = IntExt.Gcd(mtype);
        var c = wg.ElementsOrders[Cgen];
        foreach (var dim in 7.Range(1).Where(d => d != 5 && d >= mtype.Length))
        {
            foreach (var (_, _, m1) in CGenerators(m, c, dim))
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

    #region Manual search for missing Groups

    static ConcreteGroup<Mat> SearchDiagPermGL<T>(ConcreteGroup<Mat> dpgl, WordGroup g, ConcreteGroup<T> m1,
        ConcreteGroup<T> m2)
        where T : struct, IElt<T>
    {
        var allIso2 = Group.AllMorphisms(m2, dpgl, Group.MorphismType.Isomorphism).ToArray();
        foreach (var iso1 in Group.AllMorphisms(m1, dpgl, Group.MorphismType.Isomorphism))
        {
            var gens1 = m1.GetGenerators().Select(e => iso1[e]).ToArray();
            foreach (var iso2 in allIso2)
            {
                var gens = m2.GetGenerators().Select(e => iso2[e]).Concat(gens1).Distinct().ToArray();
                var matForm = Group.Generate(g.Name, dpgl, gens);
                if (matForm.IsIsomorphicTo(g))
                    return matForm;
            }
        }

        return Group.Generate(new GL(1, 2));
    }

    static ConcreteGroup<Mat> MatrixFormMissingOrder32(WordGroup g, AllSubgroups<Word> gSubgrs)
    {
        if (g.Count() != 32)
            return Group.Generate(new GL(1, 2));

        var id = FG.FindIdGroup(g, gSubgrs.Infos)[0];
        if (id.No == 43)
            return SearchDiagPermGL(GetDPSL(4, 5), g, FG.DihedralGL2p(8), FG.AbelianMat(2));
        if (id.No == 44)
            return SearchDiagPermGL(GetDPSL(4, 5), g, FG.SemiDihedralGL2p(4), FG.AbelianMat(2));
        if (id.No == 10)
            return SearchDiagPermGL(GetDPGL(3, 17), g, FG.Quaternion(8), FG.AbelianMat(4));
        if (id.No == 29)
            return SearchDiagPermGL(GetDPGL(4, 5), g, FG.MetaCyclicSdpMat(4, 4, 3), FG.AbelianMat(2));
        if (id.No == 35)
            return SearchDiagPermGL(GetDPSL(4, 5), g, FG.Quaternion(8), FG.AbelianMat(4));
        if (id.No == 7)
            return SearchDiagPermGL(GetDPSL(4, 3), g, FG.ModularMaxGL2p(4), FG.AbelianMat(2));
        if (id.No == 15)
            return SearchDiagPermGL(GetDPGL(2, 17), g, FG.AbelianMat(8), FG.AbelianMat(8));
        if (id.No == 49)
            return SearchDiagPermGL(GetDPGL(4, 3), g, ProductMatrixBlock(FG.DihedralGL2p(4), FG.AbelianMat(2)),
                FG.AbelianMat(2));
        if (id.No == 50)
            return SearchDiagPermGL(GetDPSL(4, 5), g, ProductMatrixBlock(FG.Quaternion(8), FG.AbelianMat(2)),
                FG.AbelianMat(2));

        if (id.No == 32)
        {
            var lvl = Logger.SetOff();
            var gbs0 = FG.WordGroup("(C8 x C4) : C2", "a4, c2, a2b2, a2ca2c, abacbc, abababa-1b-1, abcacaca-1cb-1");
            Logger.Level = lvl;
            var gSubgrs0 = gbs0.AllSubgroups();
            var names0 = NamesTree.BuildName(gSubgrs0.ToGroupWrapper());
            var matBs = MatrixFormFromNamesMeth2(gbs0, names0);
            return SearchDiagPermGL(matBs, g, FG.AbelianMat(4, 4), FG.AbelianMat(4));
        }
        
        if (id.No == 8)
        {
            var prod = ProductMatrixBlock(FG.SemiDihedralGL2p(4), FG.AbelianMat(2, 1));
            var gl417 = prod.Neutral().GL;
            var perms = FG.SnGensMat(gl417.N).Select(e => gl417.Create(e.Table)).ToArray();
            var s = Group.Generate("Sub-GL(4,17)", gl417, prod.GetGenerators().Concat(perms).ToArray());
            var lvl = Logger.SetOff();
            var gbs1 = FG.WordGroup("(C2 x QD16) x: C2", "b2, c2, d2, cdcd, a3dad, a2cbcb, a3ba-1b, acda-1c");
            Logger.Level = lvl;
            var matBs = SearchDiagPermGL(s, gbs1, prod, FG.AbelianMat(2));
            return SearchDiagPermGL(matBs, g, FG.ModularMaxGL2p(4), FG.AbelianMat(4));
        }
        
        return Group.Generate(new GL(1, 2));
    }

    static ConcreteGroup<Mat> MatrixFormMissingOrder48and54(WordGroup g, AllSubgroups<Word> gSubgrs)
    {
        var og = g.Count();
        if (og != 48 && og != 54)
            return Group.Generate(new GL(1, 2));

        var id = FG.FindIdGroup(g, gSubgrs.Infos)[0];
        if (og == 54 && id.No == 8)
            return SearchDiagPermGL(GetDPSL(3, 7), g, FG.AbelianMat(3, 3), FG.DihedralGL2p(3));

        if (og == 54 && id.No != 8)
            return Group.Generate(new GL(1, 2));

        if (id.No == 15)
            return SearchDiagPermGL(GetDPSL(4, 7), g, FG.DihedralGL2p(12), FG.AbelianMat(2));
        if (id.No == 16)
            return SearchDiagPermGL(GetDPSL(4, 7), g, FG.MetaCyclicSdpMat(3, 8, 2), FG.AbelianMat(2));
        if (id.No == 17)
            return SearchDiagPermGL(ProductMatrixBlock(FG.DihedralGL2p(3), FG.SemiDihedralGL2p(4)), g,
                FG.SemiDihedralGL2p(4), FG.AbelianMat(3));
        if (id.No == 41)
        {
            var gl = new GL(2, 5);
            var d8byc2 = Group.Generate("D8 x: C2", gl, gl[1, 0, 0, 4], gl[0, 1, 1, 0], gl[2, 0, 0, 2]);
            return SearchDiagPermGL(ProductMatrixBlock(d8byc2, FG.DihedralGL2p(3)), g, FG.DihedralGL2p(12),
                FG.AbelianMat(2));
        }

        if (id.No == 30)
        {
            var gl33 = new GL(3, 3);
            var a4 = Group.Generate("A4", gl33, gl33[0, 1, 0, 0, 0, 1, 1, 0, 0], gl33[0, 0, 1, 2, 0, 0, 0, 2, 0]);
            return SearchDiagPermGL(GetDPGL(3, 13), g, a4, FG.AbelianMat(4));
        }

        if (id.No == 18)
            return SearchDiagPermGL(ProductMatrixBlock(FG.DihedralGL2p(3), FG.Quaternion(16)), g, FG.Quaternion(16),
                FG.AbelianMat(3));
        if (id.No == 12)
            return SearchDiagPermGL(ProductMatrixBlock(FG.DicyclicGL2p(6), FG.MetaCyclicSdpMat(12, 2, 5)), g,
                FG.DicyclicGL2p(3), FG.AbelianMat(4));
        if (id.No == 33)
            return SearchDiagPermGL(GetDPSL(4, 5), g, FG.SL2p(3), FG.AbelianMat(2));
        if (id.No == 10)
            return SearchDiagPermGL(ProductMatrixBlock(FG.DihedralGL2p(3), FG.ModularMaxGL2p(4)), g,
                FG.ModularMaxGL2p(4),
                FG.AbelianMat(3));
        if (id.No == 28)
            return SearchDiagPermGL(GetDPSL(4, 5), g, FG.SL2p(3), FG.AbelianMat(4));
        if (id.No == 39)
        {
            var gl = new GL(2, 5);
            var d8byc2 = Group.Generate("D8 x: C2", gl, gl[1, 0, 0, 4], gl[0, 1, 1, 0], gl[2, 0, 0, 2]);
            return SearchDiagPermGL(ProductMatrixBlock(FG.DihedralGL2p(3), d8byc2), g, FG.DicyclicGL2p(6),
                FG.AbelianMat(2));
        }

        return Group.Generate(new GL(1, 2));
    }

    #endregion

    public static (WordGroup g, ConcreteGroup<Mat> mat, AllSubgroups<Mat> matSubgrs, ANameElt[] names)
        MatrixFormOfGroup(WordGroup g)
    {
        var gSubgrs = g.AllSubgroups();
        var og = g.Count();
        var names = NamesTree.BuildName(gSubgrs.ToGroupWrapper());
        var mat0 = MatrixFormFromNamesMeth2(g, names);
        if (mat0.Count() == 1)
        {
            if (og == 32)
                mat0 = MatrixFormMissingOrder32(g, gSubgrs);

            if (og == 48 || og == 54)
                mat0 = MatrixFormMissingOrder48and54(g, gSubgrs);

            if (mat0.Count() == 1)
                (mat0, _) = MatrixFormFromNames(names[0]);
        }

        mat0.Name = g.Name;
        if (mat0.GetGenerators().Any(mat => !IsDiagPerm(mat)))
            throw new($"{g} Matrix Form is not Diag Perm");

        if (!mat0.IsIsomorphicTo(g))
            throw new();

        return (g, mat0, mat0.AllSubgroups(), names);
    }

    static void MatrixFormGroupsOfOrder(int minOrd, int maxOrd)
    {
        var missing = new List<(WordGroup g, ANameElt[] names)>();
        var total = 0;
        GlobalStopWatch.Restart();
        foreach (var (g, mat, matSubgrs, names) in
                 FG.AllGroupsOfOrder(minOrd, maxOrd).Select(g => MatrixFormOfGroup(g)))
        {
            ++total;
            if (mat.Count() == 1 && g.Count() != 1)
            {
                missing.Add((g, names));
                continue;
            }

            FG.DisplayName(mat, matSubgrs, names, false, false, true, 20);
        }

        Console.WriteLine($"Missing:{missing.Count} Found:{total - missing.Count}/{total}");
        foreach (var (g, names) in missing)
        {
            names.Where(e => e is SemiDirectProductOp e0 &&
                             (e0.Lhs.ContentGroup!.GroupType == GroupType.AbelianGroup ||
                              e0.Rhs.ContentGroup!.GroupType == GroupType.AbelianGroup))
                .Println(g.ShortName);
        }

        GlobalStopWatch.Show("END");
    }

    static void MatrixFormGroupsOfOrder(int ord) => MatrixFormGroupsOfOrder(ord, ord);

    public static void GetCharacter(ConcreteGroup<Mat> mtGL, AllSubgroups<Mat> mtGLSubgrs)
    {
        var n = mtGL.Neutral().GL.N;
        var p = mtGL.Neutral().GL.P;

        var Up = FG.UnInt(p);
        var e0 = Up.GetGenerators().First();
        var cnf = Cnf.Nth(p - 1);
        var GLnC = FG.GLnK("C", n, cnf);
        var vals = mtGL.SelectMany(mat => mat.Table).Distinct().ToHashSet();
        var iso = (p - 1).Range().Select(k => (k, e0.Pow(k).K)).Where(e => vals.Contains(e.K))
            .ToDictionary(e => e.K, e => cnf.Pow(e.k).Simplify());
        iso[0] = cnf.Zero;
        var isoMt = mtGL.ToDictionary(mat => mat, mat => mat.Table.Select(z => iso[z]).ToKMatrix(n));

        var lvl = Logger.SetOff();
        // Logger.Level = LogLevel.Level1;
        var ct = FG.CharacterTableEmpty(mtGL);
        if (mtGL.GroupType == GroupType.AbelianGroup)
            ct.AbelianTable();
        else
        {
            ct.DerivedSubGroupLift();
            ct.InductionFromStabilizers();
        }

        if (!mtGL.Name.Contains("SL(2,3)"))
            ct.InductionFromSubGroups(mtGLSubgrs);
        else
        {
            var gl = mtGL.Neutral().GL;
            var m0 = new [] { 4, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0 };
            var m1 = gl.N == 4
                ? gl.Create(m0)
                : gl.Create(MatrixExt.MergeDiagonalBlocks(([1], 1), (m0, 4)));
            var mtGLSuper = Group.Generate("H", gl, mtGL.GetGenerators().Append(m1).ToArray());
            ct.RestrictionFromSuperGroup(mtGLSuper);
            
            // if (mtGL.Count() == 24)
            //     ct.SolveOrthogonality((2, 3.Range()));
            // else
            //     ct.SolveOrthogonality((2, 6.Range()));
        }

        Console.WriteLine($"Generators in {GLnC}");
        foreach (var (mat, k) in mtGL.GetGenerators().Select((mat, k) => (mat, k + 1)))
        {
            var tr = FG.PrettyPrintCnf(isoMt[mat].Trace.Simplify()).c;
            Console.WriteLine($"gen{k} of class {ct.Classes.GetClassName(mat)} Tr = {tr}");
            Console.WriteLine(isoMt[mat]);
        }

        Console.WriteLine();

        // var mtGLnC = Group.Generate(mtGL.Name, GLnC, gens.Keys.ToArray());
        // if (!mtGLnC.SetEquals(isoMt.Values))
        //     throw new();
        // if (!mtGL.IsIsomorphicTo(mtGLnC))
        //     throw new();

        ct.AllCharacters = ct.AllCharacters.Order().ToArray();

        // var isoMt = Group.IsomorphismMap(mtGL, mtGLnC, gens.ToDictionary(kv => kv.Value, kv => kv.Key));
        var map = ct.Classes.ToDictionary(cl => cl, cl => isoMt[cl].Trace);
        var chiMt = new Character<Mat>(ct.Classes, map.ToDictionary(kv => kv.Key, kv => (Cnf?)kv.Value.Simplify()));

        var chiDecomp = ct.AllCharacters.Select((chi, k) => (chi, k))
            .Select(e => (FG.InnerProduct(e.chi, chiMt).Simplify(), e.k))
            .Where(e => !e.Item1.IsZero())
            .ToArray();

        var isIrreductible = (FG.InnerProduct(chiMt, chiMt) - Cnf.CnfOne).IsZero();
        var isIsotypic = chiDecomp.All(e => e.Item1.Im.IsZero() && e.Item1.Re.IsPositiveInteger);

        if (isIrreductible)
        {
            Console.WriteLine("Irreductible Representation");
            Console.WriteLine($"Ꭓ.{chiDecomp[0].k + 1} = {chiMt}");
        }
        else
        {
            Console.WriteLine("Reductible Representation");
            Console.WriteLine($"ρ = {chiMt}");
            if (isIsotypic)
            {
                Console.WriteLine($"{chiDecomp.Length} Isotypic components");
                Console.WriteLine($"ρ = {chiDecomp.Select(e => $"{((e.Item1 - 1).IsZero() ? "" : e.Item1)}Ꭓ.{e.Item2 + 1}").Glue(" + ")}");
            }
            else
            {
                chiDecomp.Println(e => $"({FG.PrettyPrintCnf(e.Item1).c}, Ꭓ.{e.Item2 + 1})", $"Decomposition");
                Console.WriteLine();
                ct.DisplayCells(tableOnly: true);
                throw new();
            }
        }

        Console.WriteLine();
        ct.DisplayCells(tableOnly: true);

        var sum = chiDecomp.Select(e => e.Item1 * ct.AllCharacters[e.k]).Aggregate((chi0, chi1) => chi0 + chi1);
        if (!(sum - chiMt).IsZero())
            throw new();

        Logger.Level = lvl;
        Console.WriteLine();
    }

    public static void ExamplesAbelianByC2()
    {
        // Mab x: C2 Of order 24
        FG.AllAbelianGroupsOfOrder(12)
            .Select(ab => FG.Abelian(Group.AbelianGroupType(ab)))
            .SelectMany(ab => FG.AllSDPFilter(ab, FG.Abelian(2)))
            .Select(sdp => FG.WordGroup(sdp.Name, Graph.DefiningRelatorsOfGroup(sdp)))
            .Select(g => MatrixFormOfGroup(g))
            .Select(e => (e.matSubgrs, e.names))
            .DisplayNames();

        // Mab x: C3 Of order 24
        FG.AllAbelianGroupsOfOrder(8)
            .Select(ab => FG.Abelian(Group.AbelianGroupType(ab)))
            .SelectMany(ab => FG.AllSDPFilter(ab, FG.Abelian(3)))
            .Select(sdp => FG.WordGroup(sdp.Name, Graph.DefiningRelatorsOfGroup(sdp)))
            .Select(g => MatrixFormOfGroup(g))
            .Select(e => (e.matSubgrs, e.names))
            .DisplayNames();
    }

    public static void ExampleGroupOrderUpTo24()
    {
        MatrixFormGroupsOfOrder(minOrd: 1, maxOrd: 24);
    }

    public static void ExampleGroupOrder32()
    {
        MatrixFormGroupsOfOrder(ord: 32);
    }

    public static void ExampleGroupOrder48()
    {
        MatrixFormGroupsOfOrder(ord: 54);
    }

    public static void ExampleGroupOrderUpTo63()
    {
        Ring.MatrixDisplayForm = Ring.MatrixDisplay.OneLineArray;
        MatrixFormGroupsOfOrder(minOrd: 1, maxOrd: 63);
        // Missing:0 Found:319/319
        // # END Time:4m39s
    }

    public static void ExampleMetaCyclicGroupsRepresentations()
    {
        GlobalStopWatch.Restart();
        var maxOrd = 48;
        foreach (var mtGL in (maxOrd - 5).Range(6).SelectMany(o => FG.MetaCyclicSdpMat(o)))
        {
            var mtGLSubgrs = mtGL.AllSubgroups();
            var names = NamesTree.BuildName(mtGL);
            FG.DisplayName(mtGL, mtGLSubgrs, names, false, false, true, 20);
            GetCharacter(mtGL, mtGLSubgrs);
        }

        GlobalStopWatch.Show("END");
        Console.Beep();
    }

    public static void ExampleGroupsRepresentations()
    {
        Ring.DisplayPolynomial = MonomDisplay.StarCaret;
        GlobalStopWatch.Restart();
        var maxOrd = 48; // 24, 32, 63
        foreach (var (g, mtGL, matSubgrs, names) in FG.AllGroupsOfOrder(1, maxOrd).Select(sg => MatrixFormOfGroup(sg)))
        {
            FG.DisplayName(mtGL, matSubgrs, names, false, false, true, 20);
            GetCharacter(mtGL, matSubgrs);
        }

        GlobalStopWatch.Show("END");
        Console.Beep();
        
        // Representations of all groups of order up to 48
        // END Time:12m13s
        //
        // Irreductibles:141 groups
        // Reductibles with Isotypic components:109 groups
    }

    public static void ExampleIsotypicDecomposition()
    {
        GlobalStopWatch.Restart();
        // Ring.MatrixDisplayForm = Ring.MatrixDisplay.OneLineArray;
        var sl23 = FG.WordGroup("SL(2,3)", "a4, b3, ababab, a2ba2b-1");
        var c2sl23 = FG.WordGroup("C2 x SL(2,3)", "a4, b3, c2, ababab, caca-1, cbcb-1, a2ba2b-1");
        var sl23byc2 = FG.WordGroup("SL(2,3) x: C2", "a4, c3, a2b2, abab, acacac, cbc-1b-1");
        foreach (var (g, mtGL, matSubgrs, names) in new[] { sl23, c2sl23, sl23byc2 }.Select(sg => MatrixFormOfGroup(sg)))
        {
            FG.DisplayName(mtGL, matSubgrs, names, false, false, true, 20);
            GetCharacter(mtGL, matSubgrs);
        }

        GlobalStopWatch.Show("END");
        Console.Beep();
    }
}
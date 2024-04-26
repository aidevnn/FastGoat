using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Matrix;
using FastGoat.UserGroup.Perms;

namespace FastGoat.Examples;

public static class GroupMatrixForm
{
    static (int m, int n, int r)[] MetaCyclicSdp(int order)
    {
        return IntExt.Dividors(order).Where(d => d > 1)
            .SelectMany(m => FG.MetaCyclicSdpGetR(m, order / m).Select(r => (m, n: order / m, r)))
            .ToArray();
    }

    static ConcreteGroup<Mat> MetaCyclicGLnp_DiagByPerm(int m, int n, int r, int dim)
    {
        // Console.WriteLine($"Solve M({m}x:{n}){r} in GL({dim},{p})");
        var nks = new Dictionary<int, int[]>()
        {
            [2] = [2], [3] = [2, 3], [4] = [2, 3, 4], [5] = [3, 4, 5, 6], [6] = [3, 4, 5, 6, 8, 9]
        };
        foreach (var nk in nks[dim].Append(1).Where(nk => n % nk == 0).Descending())
        {
            var p = IntExt.Primes10000.First(p => (p - 1) % m == 0 && (p - 1) % (n / nk) == 0);
            var Fp = FG.UnInt(p);
            var gl = new GL(dim, p);
            var ordn = Fp.Where(e => n % Fp.ElementsOrders[e] == 0).Select(e => gl.At(gl.Neutral().Table, 0, e.K)).ToArray();
            var sn = new Sn(dim);
            var m1s = IntExt.Partitions32[dim].Where(l => l.Count == l.Distinct().Count() && l.Count <= 2)
                .OrderBy(l => l.Count)
                .Select(t => IntExt.PermAndCyclesFromType(t.Order().ToArray()))
                .Select(e =>
                {
                    var e0 = gl.Neutral().Table.Chunk(dim).ToArray();
                    var perm = sn.CreateElement(e.perm.Select(i => i + 1).ToArray());
                    var e1 = perm.Apply(e0);
                    var mat0 = gl.Create(e1.SelectMany(v => v).ToArray());
                    return ordn.Select(mat => gl.Op(mat0, mat))
                        .Where(mat => mat.IsOrder(n))
                        .Select(mat => (e.perm, e.cycles, mat));
                })
                .SelectMany(e => e);
            
            foreach (var (perm, cycles, m1) in m1s)
            {
                var m1i = gl.Invert(m1);
                var seq = cycles.Select(c => c.Length).Select(l =>
                {
                    var r0 = IntExt.PowMod(r, l, m);
                    var ordm = Fp.Where(e => m % Fp.ElementsOrders[e] == 0 && e.Pow(r0).Equals(e))
                        .OrderByDescending(e => Fp.ElementsOrders[e]);
                    return ordm.Select(a => l.Range(1).Select(k => a.Pow(IntExt.PowMod(r, k, m))).Reverse().ToArray());
                }).MultiLoop().Select(l => l.ToArray());
                foreach (var l in seq)
                {
                    var arr = new int[dim * dim];
                    foreach (var (sols, idxs) in l.Zip(cycles))
                    foreach (var (idx, sol) in idxs.Zip(sols))
                        arr[perm[idx] * (dim + 1)] = sol.K;

                    var m0 = gl.Create(arr);
                    if (m0.IsOrder(m) && gl.Op(m1i, gl.Op(m0, m1)).Equals(gl.Times(m0, r)))
                    {
                        var mtGL = Group.Generate($"M({m}x:{n}){r}", gl, m0, m1);
                        if (mtGL.Count() == m * n)
                        {
                            // Console.WriteLine($"Permutation Type[{cycles.Select(c => c.Length).Glue(",")}] in {gl}");
                            return mtGL;
                        }
                    }
                }
            }
        }

        return Group.Generate(new GL(1, 2));
    }

    static ConcreteGroup<Mat> MetaCyclicGL2p_Meth2(int m, int n, int r)
    {
        foreach (var p in IntExt.Primes10000.Where(p => m % p == 0 && (p - 1) % n == 0))
        {
            var Fp = FG.UnInt(p);
            var gl = new GL(2, p);

            var ordm = Fp.Where(e => m % Fp.ElementsOrders[e] == 0).ToArray();
            var ordn = Fp.Where(e => Fp.ElementsOrders[e] == n).ToArray();

            var m0s = ordm.Select(e => gl[e.K, 1, 0, e.K])
                .Where(mat => mat.IsOrder(m))
                .ToArray();

            var m1s = ordn.Grid2D(ordn.Append(Fp.Neutral())).Select(e => gl[e.t1.K, 0, 0, e.t2.K])
                .Where(mat => mat.IsOrder(n))
                .ToArray();

            foreach (var (m0, m1) in m0s.Grid2D(m1s).Where(e => gl.Op(gl.Invert(e.t2), gl.Op(e.t1, e.t2)).Equals(gl.Times(e.t1, r))))
            {
                var mtGl = Group.Generate($"M({m}x:{n}){r}", gl, m0, m1);
                if (mtGl.Count() == m * n)
                    return mtGl;
            }
        }

        return Group.Generate(new GL(1, 2));
    }

    static ConcreteGroup<KMatrix<ZnInt>> MetaCyclicGLnK_DiagByPerm(int m, int n, int r, int dim)
    {
        var distinctTypes = IntExt.Partitions32[dim].Where(l => l.Count == l.Distinct().Count()).Select(l => l.Order().ToArray())
            .OrderBy(l => l.Length).ToArray();
        var nks = distinctTypes.Select(l => l.Aggregate((a0, a1) => a0 * a1))
            .SelectMany(e => IntExt.Dividors(e).Append(e).Where(j => j != 1)).Append(n).ToHashSet();
        foreach (var p in nks.Select(nk => IntExt.Primes10000.First(p => (p - 1) % m == 0 && (p - 1) % nk == 0)).Distinct().Order())
        {
            var Fp = FG.UnInt(p);
            var GL = FG.GLnK($"F{p}", dim, Fp.Neutral());
            var sn = new Sn(dim);
            var o = GL.Neutral().Rows.Select(rw => rw.ToArray()).ToArray();

            var ordn = Fp.Where(e => n % Fp.ElementsOrders[e] == 0)
                .Select(e =>
                {
                    var diag = Ring.Diagonal(e.One, dim);
                    diag[0, 0] = e;
                    return new KMatrix<ZnInt>(diag);
                }).ToArray();

            var m1s = distinctTypes.Select(t => IntExt.PermAndCyclesFromType(t)).Select(e =>
            {
                var perm0 = sn.CreateElement(e.perm.Select(i => i + 1).ToArray());
                var o0 = perm0.Apply(o.Select(l => l.ToArray()).ToArray());
                var mat0 = o0.SelectMany(rw => rw).ToKMatrix(dim);
                return ordn.Select(mat1 => mat1 * mat0).Where(mat1 => mat1.IsOrder(n)).Select(mat1 => (e.perm, e.cycles, mat1));
            }).SelectMany(e => e);

            foreach (var (perm, cycles, m1) in m1s)
            {
                var sols = cycles.Select(c => c.Length).Select(l =>
                {
                    var r0 = IntExt.PowMod(r, l, m);
                    var ordm = Fp.Where(e => e.K != 1 && m % Fp.ElementsOrders[e] == 0 && e.Pow(r0).Equals(e))
                        .OrderByDescending(e => Fp.ElementsOrders[e]);

                    return ordm.Select(a => l.Range(1).Select(k => a.Pow(IntExt.PowMod(r, k, m))).Reverse().ToArray());
                }).MultiLoop().Select(l => l.ToArray()).FirstOrDefault(Array.Empty<ZnInt[]>());
                if (sols.Length == 0)
                    continue;

                var arr = new int[dim * dim];
                foreach (var (sol, idx) in sols.Zip(cycles).SelectMany(e => e.First.Zip(e.Second)))
                    arr[perm[idx] * (dim + 1)] = sol.K;

                var m0 = arr.Select(i => i * Fp.Neutral()).ToKMatrix(dim);
                if (m0.IsOrder(m) && (m1.Inv() * m0 * m1).Equals(m0.Pow(r)))
                {
                    var mtGL = Group.Generate($"M({m}x:{n}){r}", GL, m0, m1);
                    if (mtGL.Count() == m * n)
                        return mtGL;
                }
            }
        }

        return Group.Generate(FG.GLnK("F2", 1, ZnInt.ZnZero(2)));
    }
    
    static void AllGensOfMtCycSdpUpToOrder(int maxOrd, bool altGL2Meth = true)
    {
        GlobalStopWatch.Restart();
        Console.WriteLine("Start filtering MetaCyclic Groups");
        GlobalStopWatch.AddLap();
        var missing = new List<(int, int, int)>();
        var allMtCycSdp = (maxOrd - 5).Range(6).SelectMany(ord => MetaCyclicSdp(ord))
            .Select(e => (e, FG.MetaCyclicSdp(e.m, e.n, e.r).AllSubgroups())).ToDictionary(e => e.Item2, e => e.e);
        var isoMtCycSdp = allMtCycSdp.Keys.FilterIsomorphic().ToDictionary(e => e, e => allMtCycSdp[e]);
        GlobalStopWatch.Show("End filtering");

        GlobalStopWatch.AddLap();
        foreach (var (m0, e) in isoMtCycSdp)
        {
            var id = FG.FindIdGroup(m0.Parent, m0.Infos);
            var found = false;
            for (int dim = 2; dim <= 6; dim++)
            {
                var mtGL = MetaCyclicGLnp_DiagByPerm(e.m, e.n, e.r, dim);
                if (mtGL.Count() != 1)
                {
                    found = true;
                    DisplayGroup.HeadOrdersGenerators(mtGL);
                    break;
                }

                if (altGL2Meth && dim == 2)
                {
                    var mtGL2 = MetaCyclicGL2p_Meth2(e.m, e.n, e.r);
                    if (mtGL2.Count() != 1)
                    {
                        found = true;
                        DisplayGroup.HeadOrdersGenerators(mtGL2);
                        break;
                    }
                }
            }

            if (!found)
                missing.Add(e);
            else
            {
                if (id.Length != 0)
                {
                    Console.WriteLine(id[0].FullName);
                    Console.WriteLine();
                }
            }
        }

        var total = isoMtCycSdp.Count;
        missing.Println(e => $"M({e.Item1}x:{e.Item2}){e.Item3}", $"Missing:{missing.Count} Found:{total - missing.Count}/{total}");
        GlobalStopWatch.Show("Generators");
        GlobalStopWatch.Show("END");
        Console.Beep();
    }

    static void AllGensOfMtCycSdpUpToOrder(int maxOrd, int maxDim)
    {
        GlobalStopWatch.Restart();
        var missing = new List<(int, int, int)>();
        var allMtCycSdp = (maxOrd - 5).Range(6).SelectMany(ord => MetaCyclicSdp(ord)).ToArray();

        foreach (var e in allMtCycSdp)
        {
            var found = false;
            foreach (var dim in maxDim.Range(1).Where(d => d != 1 && (IntExt.Gcd(e.m, d) != 1 || IntExt.Gcd(e.m - 1, d) != 1)))
            {
                var mtGL = MetaCyclicGLnK_DiagByPerm(e.m, e.n, e.r, dim);
                if (mtGL.Count() != 1)
                {
                    found = true;
                    DisplayGroup.HeadOrdersGenerators(mtGL);
                    break;
                }
            }

            if (!found)
                missing.Add(e);
        }

        var total = allMtCycSdp.Length;
        missing.Println(e => $"M({e.Item1}x:{e.Item2}){e.Item3}", $"Missing:{missing.Count} Found:{total - missing.Count}/{total}");
        GlobalStopWatch.Show("Generators");
        Console.Beep();
    }

    static ConcreteGroup<Mat> AbelianGroup(params int[] seq)
    {
        var dim = seq.Length;
        if (dim > 5)
            throw new("Dim 5 Maximum Matrix");

        var p = IntExt.Primes10000.First(p => seq.All(o => (p - 1) % o == 0));
        var gl = new GL(dim, p);
        var a0 = IntExt.Solve_k_pow_m_equal_one_mod_n_strict(p, p - 1);
        var a = new ZnInt(p, a0);
        var Zp = Group.MulGroup($"Z({p})", a);
        var seq2 = seq.Select(o => Zp.ElementsOrders.First(e => e.Value == o).Key.K).ToArray();

        int[] Diag(int n, int k, int v) => n.Range().Grid2D()
            .Select(e => e.t1 == e.t2 ? (e.t1 == k ? v : 1) : 0).ToArray();

        var gens = seq2.Select((v, k) => gl.Create(Diag(dim, k, v))).ToArray();
        return Group.Generate(seq.Glue(" x ", "C{0}"), gl, gens);
    }

    public static void ExampleDihedalGL2p()
    {
        for (int n = 3; n < 33; n++)
        {
            var D2n = FG.DihedralGL2p(n);
            var D2pg = FG.Dihedral(n);
            D2pg.Name = $"{D2pg}pg";

            DisplayGroup.Generators(D2n);
            DisplayGroup.Generators(D2pg);

            if (!D2n.IsIsomorphicTo(D2pg))
                throw new();

            Console.WriteLine($"{D2n} IsIsomorphicTo {D2pg}");
            Console.WriteLine();
        }
    }

    public static void ExampleAbelian()
    {
        var maxP = 1;
        for (int n = 3; n < 64; n++)
        {
            var allAb = FG.AllAbelianGroupsOfOrder(n);
            Console.WriteLine($"############ Abelian groups or Order {n} ############");
            foreach (var ab in allAb)
            {
                var seq = ab.Neutral().Ei.Select(e => e.Mod).ToArray();
                var abMat = AbelianGroup(seq);

                DisplayGroup.Generators(abMat);
                DisplayGroup.Generators(ab);
                if (!ab.IsIsomorphicTo(abMat))
                    throw new();

                Console.WriteLine($"{abMat} IsIsomorphicTo {ab}");
                Console.WriteLine();

                maxP = int.Max(maxP, abMat.Neutral().GL.P);
            }

            Console.WriteLine();
        }

        Console.WriteLine(new { maxP });
    }

    public static void ExampleDicyclic()
    {
        for (int m = 2; m < 33; m++)
        {
            var Dic_m = FG.DicyclicGL2p(m);
            var Dic_m_wg = FG.DiCyclic(m);
            Dic_m_wg.Name = $"{Dic_m_wg}_wg";
            DisplayGroup.HeadGenerators(Dic_m);
            if (!Dic_m.IsIsomorphicTo(Dic_m_wg))
                throw new();

            Console.WriteLine($"{Dic_m} IsIsomorphicTo {Dic_m_wg}");
            Console.WriteLine();
        }
    }

    public static void ExampleSemiDihedralAndModularMax()
    {
        for (int n = 3; n < 9; n++)
        {
            var qd = FG.SemiDihedralGL2p(n);
            DisplayGroup.HeadGenerators(qd);
            if (!qd.IsIsomorphicTo(FG.SemiDihedral(n)))
                throw new();

            var mm = FG.ModularMaxGL2p(n);
            DisplayGroup.HeadGenerators(mm);
            if (!mm.IsIsomorphicTo(FG.ModularMax(n)))
                throw new();

            Console.WriteLine();
        }
    }

    public static void ExampleAllMetaCyclicSemiDirectProducts()
    {
        // Ring.MatrixDisplayForm = Ring.MatrixDisplay.SquareBracketNoFmt;
        Group.ActivedStorage(false);

        // AllGensOfMtCycSdpUpToOrder(32);
        // AllGensOfMtCycSdpUpToOrder(64);
        AllGensOfMtCycSdpUpToOrder(128);
        // Missing:0 Found:311/311
        // # Generators Time:1.026s
    }

    public static void ExampleAllMetaCyclicSdpUptoOrder256()
    {
        Group.ActivedStorage(false);
        // AllGensOfMtCycSdpUpToOrder(128, altGL2Meth: false);
        // Missing:1 Found:310/311
        // M(11x:10)2

        // AllGensOfMtCycSdpUpToOrder(maxOrd: 128, maxDim: 10);
        AllGensOfMtCycSdpUpToOrder(maxOrd: 256, maxDim: 12);
        // Missing:0 Found:1113/1113
        // # Generators Time:21.672s
    }
}
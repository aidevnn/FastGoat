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

    static bool IsOrder(Mat m, int o)
    {
        var gl = m.GL;
        var e = gl.Neutral();
        var mk = gl.Neutral();
        if (IntExt.PowMod(m.Det, o, gl.P) != 1)
            return false;

        for (int k = 0; k < o; k++)
        {
            if (k != 0 && mk.Equals(e))
                return false;

            mk = gl.Op(mk, m);
        }

        return mk.Equals(e);
    }

    static bool IsOrder(KMatrix<ZnInt> m, int o)
    {
        var e = m.One;
        var mk = m.One;
        if (!m.Det.Pow(o).Equals(m.KOne))
            return false;

        for (int k = 0; k < o; k++)
        {
            if (k != 0 && mk.Equals(e))
                return false;

            mk *= m;
        }

        return mk.Equals(e);
    }

    static (int[] perm, int[][] cycles) MatrixToPermutation(Mat mat)
    {
        var dim = mat.GL.N;
        if (dim.Range().Any(i => dim.Range().Count(j => mat.Table[i * dim + j] == 0) != dim - 1))
            throw new("Matrix is not a permutation");

        var perm = dim.Range();
        var set = perm.ToHashSet();
        var cycles = new List<int[]>();
        while (set.Any())
        {
            var i = set.Min();
            var lt = new List<int>();
            while (i != lt.FirstOrDefault(-1))
            {
                lt.Add(i);
                set.Remove(i);
                var i0 = i;
                i = dim.Range().First(j => mat.Table[i0 * dim + j] != 0);
                perm[i0] = i;
            }

            cycles.Add(lt.ToArray());
        }

        return (perm, cycles.ToArray());
    }

    static IEnumerable<Mat> GetGLnPermutations(int dim, int n, int p)
    {
        var Fp = FG.UnInt(p);
        var gl = new GL(dim, p);
        var ordn = Fp.Where(e => n % Fp.ElementsOrders[e] == 0).ToArray();
        return ordn.Select(e => e.K)
            .SelectMany(e =>
            {
                if (gl.N == 2)
                    return new[]
                    {
                        gl[0, e, 1, 0] // type[2]
                    };
                if (gl.N == 3)
                    return new[]
                    {
                        gl[0, e, 0, 0, 0, 1, 1, 0, 0], // type[3]
                        gl[e, 0, 0, 0, 0, 1, 0, 1, 0], // type[1,2]
                    };
                if (gl.N == 4)
                    return new[]
                    {
                        gl[0, e, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0], // type[4]
                        gl[e, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0], // type[1,3]
                        // gl[0, e, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0], // type[2,2]
                    };
                if (gl.N == 5)
                    return new[]
                    {
                        gl[0, e, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0], // type[5]
                        gl[e, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0], // type[1,4]
                        gl[0, e, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0], // type[2,3]
                    };
                if (gl.N == 6)
                    return new[]
                    {
                        gl[0, e, 0, 0, 0, 0,
                            0, 0, 1, 0, 0, 0,
                            0, 0, 0, 1, 0, 0,
                            0, 0, 0, 0, 1, 0,
                            0, 0, 0, 0, 0, 1,
                            1, 0, 0, 0, 0, 0], // type[6]
                        // gl[0, e, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0], // type[3,3]
                        // gl[0, e, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0], // type[2.4]
                        // gl[e, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0], // type[1,5]
                    };

                throw new();
            })
            .Where(mat => IsOrder(mat, n));
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
            var m1s = GetGLnPermutations(dim, n, p);
            foreach (var m1 in m1s)
            {
                var gl = m1.GL;
                var m1i = gl.Invert(m1);
                var (perm, cycles) = MatrixToPermutation(m1);
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
                    if (IsOrder(m0, m) && gl.Op(m1i, gl.Op(m0, m1)).Equals(gl.Times(m0, r)))
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
                .Where(mat => IsOrder(mat, m))
                .ToArray();

            var m1s = ordn.Grid2D(ordn.Append(Fp.Neutral())).Select(e => gl[e.t1.K, 0, 0, e.t2.K])
                .Where(mat => IsOrder(mat, n))
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
                return ordn.Select(mat1 => mat1 * mat0).Where(mat1 => IsOrder(mat1, n)).Select(mat1 => (e.perm, e.cycles, mat1));
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
                if (IsOrder(m0, m) && (m1.Inv() * m0 * m1).Equals(m0.Pow(r)))
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
                    if (id.Length != 0)
                    {
                        Console.WriteLine(id[0].FullName);
                        Console.WriteLine();
                    }

                    if (!mtGL.IsIsomorphicTo(m0.Parent))
                        throw new();

                    break;
                }

                if (altGL2Meth && dim == 2)
                {
                    var mtGL2 = MetaCyclicGL2p_Meth2(e.m, e.n, e.r);
                    if (mtGL2.Count() != 1)
                    {
                        found = true;
                        DisplayGroup.HeadOrdersGenerators(mtGL2);
                        if (id.Length != 0)
                        {
                            Console.WriteLine(id[0].FullName);
                            Console.WriteLine();
                        }

                        if (!mtGL2.IsIsomorphicTo(m0.Parent))
                            throw new();

                        break;
                    }
                }
            }

            if (!found)
                missing.Add(e);
        }

        var total = isoMtCycSdp.Count;
        missing.Println(e => $"M({e.Item1}x:{e.Item2}){e.Item3}", $"Missing:{missing.Count} Found:{total - missing.Count}/{total}");
        GlobalStopWatch.Show("End Gens");
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
        GlobalStopWatch.Show("End Gens");
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

    static int OrderMatOrth(ZnInt x0, ZnInt y0)
    {
        var (x1, y1) = (x0.One, y0.Zero);
        for (int i = 1; i < 1000; i++)
        {
            (x1, y1) = (x0 * x1 - y0 * y1, x0 * y1 + y0 * x1);
            if (x1.Equals(x1.One) && y1.IsZero())
                return i;
        }

        throw new("####################################");
    }

    static ConcreteGroup<Mat> DihedralGO2p(int n = 16)
    {
        int p = 2;
        ZnInt a, x = ZnInt.ZnZero(p), y = ZnInt.ZnZero(p);
        foreach (var p0 in IntExt.Primes10000.Where(p0 => (p0 - 1) % n == 0))
        {
            p = p0;
            var a0 = IntExt.Solve_k_pow_m_equal_one_mod_n_strict(p, p - 1);
            a = new ZnInt(p, a0);
            var Zp = Group.MulGroup($"F{p}", a);
            var square = Zp.Append(a.Zero).Select(x0 => (x: x0, x2: x0 * x0)).GroupBy(e => e.x2)
                .ToDictionary(e => e.Key, e => e.Select(f => f.x).ToArray());
            var dicSquare = Zp.ToDictionary(x0 => x0, x0 => square.ContainsKey(x0) ? square[x0] : []);
            dicSquare[a.Zero] = [];
            var XYs = Zp.Append(a.Zero)
                .Select(x0 => (x: x0, yList: dicSquare[1 - x0 * x0]))
                .Where(e => e.yList.Length != 0)
                .SelectMany(e => e.yList.Select(y0 => (e.x, y: y0)))
                .Distinct()
                .Select(e => (e.x, e.y))
                .Where(e => OrderMatOrth(e.x, e.y) == n)
                .OrderBy(e => e.x.K)
                .ToArray();

            if (XYs.Length != 0)
            {
                (x, y) = XYs[0];
                break;
            }
        }

        var gl = new GL(2, p);
        var m0 = gl[x.K, y.K, (-y).K, x.K];
        var m1 = gl[0, 1, 1, 0];

        return Group.Generate($"D{2 * n}", gl, m0, m1);
    }

    public static void ExampleDihedalGO2p()
    {
        for (int n = 3; n < 33; n++)
        {
            var D2n = DihedralGO2p(n);
            var D2pg = FG.Dihedral(n);
            D2pg.Name = $"{D2pg}pg";

            DisplayGroup.Generators(D2n);
            DisplayGroup.Generators(D2pg);

            DisplayGroup.AreIsomorphics(D2n, D2pg);
            Console.WriteLine();
        }
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
    }

    public static void ExampleGL7p()
    {
        Group.ActivedStorage(false);
        // AllGensOfMtCycSdpUpToOrder(128, altGL2Meth: false);
        // Missing:1 Found:310/311
        // M(11x:10)2
        
        // AllGensOfMtCycSdpUpToOrder(maxOrd: 128, maxDim: 10);
        AllGensOfMtCycSdpUpToOrder(maxOrd: 256, maxDim: 12);
    }
}
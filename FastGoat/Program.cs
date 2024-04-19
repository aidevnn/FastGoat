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
using GroupRegX = System.Text.RegularExpressions.Group;
using OrdMats = System.Collections.Generic.Dictionary<int, FastGoat.UserGroup.Matrix.Mat[]>;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");

(int m, int n, int r)[] MetaCyclicSdp(int order)
{
    var ms = IntExt.Dividors(order).Where(d => d > 1).ToArray();

    HashSet<(int m, int n, int r)> all = new();
    foreach (var m in ms)
    {
        var n = order / m;
        var rs = IntExt.SolveAll_k_pow_m_equal_one_mod_n(m, n);
        foreach (var r in rs)
            all.Add((m, n, r));
    }

    return all.ToArray();
}

(int m, int n, int r)[] FrobeniusSdp(int order)
{
    var ms = IntExt.Dividors(order).Where(d => d > 1 && d % 2 == 1).ToArray();

    List<(int m, int n, int r)> all = new();
    foreach (var m in ms)
    {
        var n = order / m;
        foreach (var r in FG.FrobeniusGetR(m, n))
            all.Add((m, n, r));
    }

    return all.ToArray();
}

ConcreteGroup<Mat> MetaCyclicGL2p_Meth1(int m, int n, int r)
{
    var mtSdp = FG.MetaCyclicSdp(m, n, r);
    foreach (var p in Primes10000.Where(p => (p - 1) % m == 0 && (p - 1) % n == 0).Take(10))
    {
        var gl = new GL(2, p);
        var gl3 = new GL(3, p);
        var x = Solve_k_pow_m_equal_one_mod_n_strict(p, p - 1);
        var Fp = Group.MulGroup($"F{p}", new ZnInt(p, x));
        var ordms = Fp.Where(e => Fp.ElementsOrders[e] == m).Select(e => e.K).ToArray();
        var ordns = Fp.Where(e => Fp.ElementsOrders[e] == n).Select(e => e.K).ToArray();

        var a0s = ordms.Grid2D().Where(e => PowMod(e.t1, r, p) == e.t2 && PowMod(e.t2, r, p) == e.t1).ToArray();
        var b0s = ordns.Grid2D().Where(e => Group.Cycle(gl, gl[0, e.t1, e.t2, 0]).Count == n).ToArray();

        if (a0s.Length == 0 || b0s.Length == 0)
            continue;

        foreach (var (a0, a1) in a0s)
        {
            var a2 = ordms.Except([a0, a1]).FirstOrDefault(e => Fp.ElementsOrders[new(p, e)] == m, a0);
            foreach (var (b0, b1) in b0s)
            {
                {
                    var m0 = gl[a0, 0, 0, a1];
                    var m1 = gl[0, b0, b1, 0];
                    var mtGl = Group.Generate($"M({m}x:{n}){r}", gl, m0, m1);
                    if (mtGl.IsIsomorphicTo(mtSdp))
                        return mtGl;
                }

                {
                    var m03 = gl3[a2, 0, 0, 0, a0, 0, 0, 0, a1];
                    var m13 = gl3[1, 0, 0, 0, 0, b0, 0, b1, 0];
                    var mtGl3 = Group.Generate($"M({m}x:{n}){r}", gl3, m03, m13);

                    if (mtGl3.IsIsomorphicTo(mtSdp))
                        return mtGl3;
                }

                {
                    var m03 = gl3[1, 0, 0, 0, a0, 0, 0, 0, a1];
                    var m13 = gl3[b0, 0, 0, 0, 0, b0, 0, b1, 0];
                    var mtGl3 = Group.Generate($"M({m}x:{n}){r}", gl3, m03, m13);

                    if (mtGl3.IsIsomorphicTo(mtSdp))
                        return mtGl3;
                }
            }
        }
    }

    var gl12 = new GL(1, 2);
    return Group.Generate("()", gl12, gl12.Neutral());
}

bool IsOrder(Mat m, int o)
{
    var gl = m.GL;
    var e = gl.Neutral();
    var mk = gl.Neutral();
    for (int k = 0; k < o; k++)
    {
        if (k != 0 && mk.Equals(e))
            return false;

        mk = gl.Op(mk, m);
    }

    return mk.Equals(e);
}

ConcreteGroup<Mat> MetaCyclicGL2p_Meth2(int m, int n, int r)
{
    foreach (var p in Primes10000.Where(p => m % p == 0 && (p - 1) % n == 0))
    {
        var Fp = FG.UnInt(p);
        var gl = new GL(2, p);
        var i2 = gl.Neutral();

        var ordm = Fp.Where(e => m % Fp.ElementsOrders[e] == 0).ToArray();
        var ordn = Fp.Where(e => Fp.ElementsOrders[e] == n).ToArray();

        var m0s = ordm.Grid3D(ordm, Fp.Prepend(Fp.Neutral().Zero).OrderBy(e => e.K))
            .Select(e => gl[e.t1.K, e.t3.K, 0, e.t2.K])
            .Where(mat => IsOrder(mat, m))
            .ToArray();

        var m1s = ordn.Grid2D(ordn.Append(Fp.Neutral())).Select(e => gl[e.t1.K, 0, 0, e.t2.K])
            .Where(mat => IsOrder(mat, n))
            .ToArray();

        var (m0, m1) = m0s.Grid2D(m1s)
            .FirstOrDefault(e => gl.Op(gl.Invert(e.t2), gl.Op(e.t1, e.t2)).Equals(gl.Times(e.t1, r)), (t1: i2, t2: i2));

        if (m0.Equals(i2) && m1.Equals(i2))
            continue;
        
        var mtGL = Group.Generate($"M({m}x:{n}){r}", gl, m0, m1);
        if (mtGL.IsIsomorphicTo(FG.MetaCyclicSdp(m, n, r)))
            return mtGL;
    }

    var gl12bs = new GL(1, 2);
    return Group.Generate("()", gl12bs, gl12bs.Neutral());
}

void Run(int maxOrd = 32, bool frob = false)
{
    GlobalStopWatch.Restart();
    Func<int, (int m, int n, int r)[]> f = frob ? FrobeniusSdp : MetaCyclicSdp;

    var missing = new List<(int, int, int)>();
    var allMtCycSdp = (maxOrd - 5).Range(6).SelectMany(ord => f(ord))
        .Select(e => (e, FG.MetaCyclicSdp(e.m, e.n, e.r).AllSubgroups())).ToDictionary(e => e.Item2, e => e.e);
    var isoMtCycSdp = allMtCycSdp.Keys.FilterIsomorphic().ToDictionary(e => e, e => allMtCycSdp[e]);

    foreach (var (_, e) in isoMtCycSdp)
    {
        var mtGL1 = MetaCyclicGL2p_Meth1(e.m, e.n, e.r);
        if (mtGL1.Count() != 1)
            DisplayGroup.HeadOrdersGenerators(mtGL1);
        else
        {
            var mtGL2 = MetaCyclicGL2p_Meth2(e.m, e.n, e.r);
            if (mtGL2.Count() != 1)
                DisplayGroup.HeadOrdersGenerators(mtGL2);
            else
                missing.Add(e);
        }
    }

    var total = isoMtCycSdp.Count;
    missing.Println(e => $"M({e.Item1}x:{e.Item2}){e.Item3}", $"Missing:{missing.Count} Found:{total - missing.Count}/{total}");
    Console.WriteLine($"var missing = new [] {{ {missing.Glue(", ")} }};");
    GlobalStopWatch.Show("END");
    Console.Beep();
}

{
    Ring.MatrixDisplayForm = Ring.MatrixDisplay.SquareBracketNoFmt;
    
    // Run();
    // Run(frob: true);
    
    Run(maxOrd:64);
    // Run(maxOrd:96);
}

// Missing:9 Found:102/111
//     M(9x:3)4
//     M(5x:8)2
//     M(9x:6)2
//     M(9x:6)4
//     M(5x:12)2
//     M(15x:4)2
//     M(7x:9)2
//     M(16x:4)3
//     M(16x:4)5
// var missing = new [] { (9, 3, 4), (5, 8, 2), (9, 6, 2), (9, 6, 4), (5, 12, 2), (15, 4, 2), (7, 9, 2), (16, 4, 3), (16, 4, 5) };
// # END Time:9.438s
// 
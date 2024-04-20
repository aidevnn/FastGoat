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

bool IsOrder(Mat m, int o)
{
    var gl = m.GL;
    var e = gl.Neutral();
    var mk = gl.Neutral();
    if (PowMod(m.Det, o, gl.P) != 1)
        return false;

    for (int k = 0; k < o; k++)
    {
        if (k != 0 && mk.Equals(e))
            return false;

        mk = gl.Op(mk, m);
    }

    return mk.Equals(e);
}

ConcreteGroup<Mat> MetaCyclicGL2p_Meth1(int m, int n, int r)
{
    var mtSdp = FG.MetaCyclicSdp(m, n, r);
    var p = Primes10000.First(p => (p - 1) % m == 0 && (p - 1) % n == 0);

    var gl = new GL(2, p);
    var Fp = FG.UnInt(p);
    var ordms = Fp.Where(e => Fp.ElementsOrders[e] == m).Select(e => e.K).ToArray();
    var ordns = Fp.Where(e => Fp.ElementsOrders[e] == n).Select(e => e.K).ToArray();

    var a0s = ordms.Grid2D().Where(e => PowMod(e.t1, r, p) == e.t2 && PowMod(e.t2, r, p) == e.t1).ToArray();
    var b0s = ordns.Grid2D().Where(e => IsOrder(gl[0, e.t1, e.t2, 0], n)).ToArray();

    foreach (var (a0, a1) in a0s)
    {
        foreach (var (b0, b1) in b0s)
        {
            {
                var m0 = gl[a0, 0, 0, a1];
                var m1 = gl[0, b0, b1, 0];
                var mtGl = Group.Generate($"M({m}x:{n}){r}", gl, m0, m1);
                if (mtGl.Count() == m * n)
                    return mtGl;

                // if (mtGl.IsIsomorphicTo(mtSdp))
                //     return mtGl;
            }
        }
    }

    return Group.Generate(new GL(1, 2));
}

ConcreteGroup<Mat> MetaCyclicGL2p_Meth2(int m, int n, int r)
{
    var mtSdp = FG.MetaCyclicSdp(m, n, r);
    foreach (var p in Primes10000.Where(p => m % p == 0 && (p - 1) % n == 0))
    {
        var Fp = FG.UnInt(p);
        var gl = new GL(2, p);

        var ordm = Fp.Where(e => m % Fp.ElementsOrders[e] == 0).ToArray();
        var ordn = Fp.Where(e => Fp.ElementsOrders[e] == n).ToArray();

        var m0s = ordm.Grid3D(ordm, Fp.Prepend(Fp.Neutral().Zero).OrderBy(e => e.K))
            .Select(e => gl[e.t1.K, e.t3.K, 0, e.t2.K])
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

            // if (mtGl.IsIsomorphicTo(mtSdp))
            //     return mtGl;
        }
    }

    return Group.Generate(new GL(1, 2));
}

IEnumerable<Mat> GL3MatrixPerm(GL gl, (ZnInt t1, ZnInt t2, ZnInt t3) e, int ord) => new[]
    {
        gl[e.t1.K, 0, 0, 0, 0, e.t2.K, 0, e.t3.K, 0],
        gl[0, e.t1.K, 0, 0, 0, e.t2.K, e.t3.K, 0, 0]
    }
    .Where(mat => IsOrder(mat, ord));

ConcreteGroup<Mat> MetaCyclicGL3p_Meth(int m, int n, int r)
{
    var mtSdp = FG.MetaCyclicSdp(m, n, r);
    var p = Primes10000.First(p => (p - 1) % m == 0 && (p - 1) % n == 0);
    var gl = new GL(3, p);
    var Fp = FG.UnInt(p);
    var ordm = Fp.Where(e => m % Fp.ElementsOrders[e] == 0).OrderBy(e => e.K).ToArray();
    var ordn = Fp.Where(e => n % Fp.ElementsOrders[e] == 0).OrderBy(e => e.K).ToArray();

    var m0s = ordm.Grid3D().Select(e => gl[e.t1.K, 0, 0, 0, e.t2.K, 0, 0, 0, e.t3.K])
        .Where(mat => IsOrder(mat, m))
        .ToArray();
    var m1s = ordn.Grid3D().SelectMany(e => GL3MatrixPerm(gl, e, n)).ToArray();

    foreach (var (m0, m1) in m0s.Grid2D(m1s).Where(e => gl.Op(gl.Invert(e.t2), gl.Op(e.t1, e.t2)).Equals(gl.Times(e.t1, r))))
    {
        var mtGl = Group.Generate($"M({m}x:{n}){r}", gl, m0, m1);
        if (mtGl.Count() == m * n)
            return mtGl;

        // if (mtGl.IsIsomorphicTo(mtSdp))
        //     return mtGl;
    }

    return Group.Generate(new GL(1, 2));
}

ConcreteGroup<Mat> MetaCyclicGL4p_Meth(int m, int n, int r)
{
    if (n % 4 != 0)
        return Group.Generate(new GL(1, 2));

    var p = Primes10000.First(p => (p - 1) % m == 0 && (p - 1) % (n / 4) == 0);

    var Fp = FG.UnInt(p);
    var Fp4 = Fp.Where(e => e.Pow(m).K == 1)
        .MultiLoop(4)
        .Select(l => l.ToArray()).Select(l => (t1: l[0], t2: l[1], t3: l[2], t4: l[3]))
        .OrderBy(l => l)
        .ToArray();
    
    var gl = new GL(4, p);
    var m1s = Fp.Where(e => Fp.ElementsOrders[e] == n / 4)
        .SelectMany(e => new[]
        {
            gl[0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, e.K, 0],
            gl[0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, e.K, 0, 0, 0]
        })
        .Where(mat => IsOrder(mat, n))
        .ToArray();

    var m0s = Fp4.Select(e => gl[e.t1.K, 0, 0, 0, 0, e.t2.K, 0, 0, 0, 0, e.t3.K, 0, 0, 0, 0, e.t4.K])
        .Where(mat => IsOrder(mat, m))
        .ToArray();

    foreach (var (m0, m1) in m0s.Grid2D(m1s).Where(e => gl.Op(gl.Invert(e.t2), gl.Op(e.t1, e.t2)).Equals(gl.Times(e.t1, r))))
    {
        var mtGl = Group.Generate($"M({m}x:{n}){r}", gl, m0, m1);
        if (mtGl.Count() == m * n)
            return mtGl;

        // if (mtGl.IsIsomorphicTo(mtSdp))
        //     return mtGl;
    }

    return Group.Generate(new GL(1, 2));
}

void Run(int maxOrd = 32, bool frob = false)
{
    GlobalStopWatch.Restart();
    Func<int, (int m, int n, int r)[]> f = frob ? FrobeniusSdp : MetaCyclicSdp;

    var missing = new List<(int, int, int)>();
    var allMtCycSdp = (maxOrd - 5).Range(6).SelectMany(ord => f(ord))
        .Select(e => (e, FG.MetaCyclicSdp(e.m, e.n, e.r).AllSubgroups())).ToDictionary(e => e.Item2, e => e.e);
    var isoMtCycSdp = allMtCycSdp.Keys.FilterIsomorphic().ToDictionary(e => e, e => allMtCycSdp[e]);

    foreach (var (mt, e) in isoMtCycSdp)
    {
        var idGap = FG.FindIdGroup(mt.Parent, mt.Infos)[0].FullName;
        var mtGLmeth1 = MetaCyclicGL2p_Meth1(e.m, e.n, e.r);
        if (mtGLmeth1.Count() != 1)
        {
            DisplayGroup.HeadOrdersGenerators(mtGLmeth1);
            Console.WriteLine(idGap);
            Console.WriteLine();
        }
        else
        {
            var mtGLmeth2 = MetaCyclicGL2p_Meth2(e.m, e.n, e.r);
            if (mtGLmeth2.Count() != 1)
            {
                DisplayGroup.HeadOrdersGenerators(mtGLmeth2);
                Console.WriteLine(idGap);
                Console.WriteLine();
            }
            else
            {
                var mtGLmeth3 = MetaCyclicGL3p_Meth(e.m, e.n, e.r);
                if (mtGLmeth3.Count() != 1)
                {
                    DisplayGroup.HeadOrdersGenerators(mtGLmeth3);
                    Console.WriteLine(idGap);
                    Console.WriteLine();
                }
                else
                {
                    var mtGLmeth4 = MetaCyclicGL4p_Meth(e.m, e.n, e.r);
                    if (mtGLmeth4.Count() != 1)
                    {
                        DisplayGroup.HeadOrdersGenerators(mtGLmeth4);
                        Console.WriteLine(idGap);
                        Console.WriteLine();
                    }
                    else
                        missing.Add(e);
                }
            }
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

// Missing:1 Found:110/111
//     M(9x:6)2
// var missing = new [] { (9, 6, 2) };
// # END Time:3.041s
// 

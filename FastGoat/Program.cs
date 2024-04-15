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

bool BruteForce<T>(ConcreteGroup<T> gl, int m, int n, int r) where T : struct, IElt<T>
{
    var m0s = gl.Where(e => gl.ElementsOrders[e] == m).ToArray();
    var m1s = gl.Where(e => gl.ElementsOrders[e] == n).ToArray();
    Console.WriteLine($"###### Start search generators of M({m}x:{n}){r} in {gl}");
    foreach (var m0 in m0s)
    {
        var m0r = gl.Times(m0, r);
        foreach (var m1 in m1s)
        {
            if (!gl.Op(gl.Invert(m1), gl.Op(m0, m1)).Equals(m0r))
                continue;

            var gmn = Group.Generate($"M({m}x:{n}){r}", gl, m0, m1);
            if (gmn.Count() != m * n)
                continue;

            DisplayGroup.HeadGenerators(gmn);
            Console.WriteLine();
            return true;
        }
    }

    return false;
}

{
    GlobalStopWatch.Restart();
    GlobalStopWatch.AddLap();
    var gls0 = new List<ConcreteGroup<Mat>>();
    foreach (var p in Primes10000.Take(8))
    {
        var (a, b) = FG.GLnpGenerators(2, p);
        var gl = Group.Generate(a.GL, a, b);
        gls0.Add(gl);
        DisplayGroup.HeadOrders(gl);
    }

    foreach (var (n, p) in new[] { (3, 2), (3, 3), (4, 2) })
    {
        var (a, b) = FG.GLnpGenerators(n, p);
        var gl = Group.Generate(a.GL, a, b);
        gls0.Add(gl);
        DisplayGroup.HeadOrders(gl);
    }

    foreach (var p in new[] { 5, 7, 11, 13, 17, 19, 23, 29 })
    {
        var go = FG.GO3p(p);
        gls0.Add(go);
        DisplayGroup.HeadOrders(go);
    }

    gls0 = gls0.OrderBy(g => g.Neutral().GL.N).ThenBy(g => g.Neutral().GL.P).ToList();
    
    GlobalStopWatch.Show();

    GlobalStopWatch.AddLap();
    var (nbFound, nbTotal) = (0, 0);

    for (int ord = 6; ord < 33; ord++)
    {
        foreach (var (m, n, r) in MetaCyclicSdp(ord))
        {
            var found = false;
            if ((n == 2 && r == m - 1) || (n == 4 && m % 2 == 1 && r == m - 1))
                continue;
            else
            {
                ++nbTotal;
                foreach (var gl in gls0)
                {
                    if (BruteForce(gl, m, n, r))
                    {
                        found = true;
                        ++nbFound;
                        break;
                    }
                }
            }

            if (!found)
                Console.WriteLine($"########################### generators of M({m}x:{n}){r} not found ###########################\n");
        }
    }

    GlobalStopWatch.Show($"Found:{nbFound} Total:{nbTotal}");
    Console.Beep();
}
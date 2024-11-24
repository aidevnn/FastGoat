using System.CodeDom;
using System.Collections;
using System.Diagnostics;
using System.Numerics;
using System.Reflection.Emit;
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
using FastGoat.UserGroup.DatabaseSmallGroups;
using FastGoat.UserGroup.EllCurve;
using FastGoat.UserGroup.Floats;
using FastGoat.UserGroup.GModuleN;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Matrix;
using FastGoat.UserGroup.GModuleN;
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

(int y, bool sol) ApproxSolver(int x, int n)
{
    var y2 = BigInteger.Pow(x, 3) + n;
    if (y2 < 0)
        return (0, false);

    var ya = BigInteger.Parse($"{double.Ceiling(double.Sqrt((double)y2))}");
    while (ya * ya < y2 + 1)
    {
        if (ya * ya == y2)
            return ((int)ya, true);

        ya = (ya + y2 + 1) / 2;
    }

    return (0, false);
}

{
    int[] B081119 =
    [
        0, 5, 2, 2, 2, 2, 0, 0, 7, 10, 2, 0, 4, 0, 0, 4, 2, 16, 2, 2, 0, 0, 2, 0, 8, 2, 2, 1, 4, 0, 2, 2, 0, 2, 0, 2, 8,
        6, 2, 0, 2, 2, 0, 2, 4, 0, 0, 0, 2, 2, 2, 0, 2, 0, 2, 2, 2, 6, 0, 0, 0, 0, 0, 4, 5, 8, 0, 0, 4, 0, 0, 2, 2, 12,
        0, 0, 2, 0, 0, 2, 8, 2, 2, 0, 0, 0, 0, 0, 0, 8, 0, 2, 2, 0, 2, 0, 0, 2, 2, 2, 12, 4, 0, 0, 0, 2, 2, 2, 8, 0, 0,
        0, 2, 12, 0, 0, 0, 2, 0, 2, 2, 4, 2, 0, 0, 1, 2, 2, 4, 6, 0, 2, 2, 0, 2, 0, 2, 0, 2, 0, 0, 6, 2, 2, 2, 8, 0, 0,
        4, 0, 2, 2, 2, 0, 2, 0, 2, 0, 0, 0, 0, 8, 0, 2, 6, 0, 0, 0, 4, 6, 2, 6, 2, 0, 2, 0, 0, 2, 0, 0, 0, 0, 0, 0, 2,
        2, 2, 0, 2, 2, 2, 0, 2, 0, 0, 2, 6, 4, 4, 2, 0, 0, 0, 0, 4, 0, 2, 0, 4, 0, 0, 0, 0, 0, 0, 0, 1, 10, 0, 0, 4, 0,
        0, 2, 2, 26, 2, 0, 0, 2, 0, 0, 4, 8, 0, 0, 2, 0, 0, 0, 0, 2, 0, 0, 0, 0, 2, 0, 4, 2, 0, 2, 10, 0, 0, 2, 2, 2, 0,
        0, 8, 0, 2, 0, 2, 2, 0, 0, 2, 4, 0, 2, 2, 0, 0, 2, 0, 0, 0, 0, 2, 6, 2, 2, 0, 0, 0, 0, 2, 6, 2, 0, 0, 0, 4, 2,
        2, 18, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 14, 0, 0, 0, 2, 2, 0, 2, 2, 2, 0, 0, 0, 2, 0, 0, 2,
        2, 2, 0, 2, 4, 0, 0, 0, 0, 0, 3, 2, 0, 4, 2, 0, 0, 4, 2, 2, 8, 0, 0, 2, 0, 0, 6, 10, 2, 2, 0, 0, 0, 2, 0, 2, 4,
        0, 0, 0, 2, 0, 0, 0, 8, 0, 2, 2, 2, 0, 0, 0, 2, 2, 0, 12, 2, 0, 0, 6, 2, 0, 0, 2, 0, 0, 2, 2, 2, 0, 0, 4, 0, 2,
        2, 2, 4, 0, 0, 2, 0, 6, 0, 2, 0, 0, 0, 2, 0, 0, 2, 2, 2, 0, 6, 0, 0, 0, 0, 0, 6, 0, 0, 0, 0, 0, 0, 2, 8, 2, 4,
        2, 0, 0, 0, 0, 8, 0, 2, 0, 0, 0, 0, 0, 2, 0, 0, 2, 0, 0, 0, 6, 2, 0, 0, 2, 2, 0, 0, 2, 0, 0, 0, 2, 0, 0, 0, 0,
        4, 0, 2, 2, 4, 2, 2, 0, 0, 0, 0, 4, 0, 0, 0, 2, 0, 2, 0, 2, 0, 2, 0, 2, 8, 2, 2, 0, 0, 0, 2, 9, 6, 0, 0, 6, 0,
        2, 2, 0, 8, 0, 0, 0, 2, 0, 0, 8, 2, 2, 2, 0, 0, 0, 0, 0, 8, 0, 2, 4, 0, 0, 0, 2, 2, 0, 0, 10, 2, 0, 2, 0, 0, 0,
        0, 4, 2, 0, 0, 0, 4, 0, 2, 0, 0, 2, 2, 14, 0, 0, 0, 0, 0, 2, 4, 10, 6, 0, 0, 0, 0, 2, 0, 2, 0, 0, 0, 2, 0, 0, 0,
        0, 6, 2, 0, 0, 0, 2, 2, 2, 0, 2, 2, 2, 2, 0, 0, 0, 2, 0, 0, 6, 0, 0, 0, 2, 2, 4, 0, 2, 0, 2, 0, 2, 6, 2, 0, 0,
        0, 0, 2, 2, 6, 0, 2, 0, 0, 2, 0, 4, 0, 0, 0, 2, 0, 0, 0, 0, 8, 0, 0, 2, 0, 2, 0, 4, 2, 0, 4, 0, 0, 0, 0, 0, 6,
        0, 2, 2, 0, 0, 0, 0, 2, 0, 2, 2, 2, 0, 0, 0, 8, 0, 0, 6, 0, 2, 2, 0, 2, 0, 0, 0, 0, 2, 0, 2, 0, 0, 0, 2, 4, 4,
        2, 2, 0, 2, 0, 2, 0, 2, 0, 2, 2, 0, 2, 4, 0, 0, 2, 2, 6, 0, 0, 0, 0, 0, 0, 2, 5, 6, 0, 0, 2, 0, 0, 0, 8, 2, 2,
        4, 0, 0, 2, 0, 8, 2, 0, 2, 0, 2, 0, 0, 0, 2, 0, 2, 4, 0, 2, 0, 0, 0, 0, 2, 2, 2, 0, 4, 0, 0, 0, 2, 0, 0, 4, 2,
        2, 2, 0, 0, 0, 0, 2, 8, 4, 0, 0, 0, 0, 0, 0, 8, 6, 0, 0, 0, 0, 0, 0, 0, 6, 0, 0, 4, 0, 2, 0, 2, 0, 2, 2, 2, 2,
        2, 0, 0, 2, 0, 0, 0, 0, 0, 0, 2, 0, 0, 2, 2, 4, 0, 0, 0, 2, 2, 0, 10, 0, 0, 0, 2, 6, 2, 0, 0, 0, 0, 0, 2, 6, 4,
        0, 0, 0, 2, 0, 2, 4, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 4, 0, 0, 2, 2, 18, 0, 0, 2, 0, 2, 0, 0, 0, 2, 0, 0, 4, 0, 0,
        0, 2, 2, 0, 10, 0, 0, 0, 0, 2, 2, 6, 2, 4, 0, 0, 0, 2, 2, 0, 2, 2, 0, 0, 2, 0, 0, 0, 0, 0, 0, 2, 2, 0, 0, 0, 0,
        2, 0, 4, 0, 0, 0, 0, 2, 0, 2, 0, 2, 0, 0, 2, 6, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 2, 6, 4, 0, 0, 0, 0, 0, 6, 2,
        2, 0, 6, 0, 2, 2, 2, 2, 0, 2, 0, 2, 0, 0, 2, 2, 0, 0, 0, 2, 0, 0, 0, 4, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 2, 2, 0,
        0, 5
    ];

    GlobalStopWatch.Restart();
    var set = new HashSet<(BigInteger, BigInteger, int)>();
    int[][] xmax10e4 = [[10000], [17, 24, 100, 141, 217, 388, 414, 513, 516, 521, 568, 649, 659, 740, 757, 836, 960, 985]];
    int[][] xmax10e5 = [[100000], [297, 377, 427, 885, 899]];
    int[][] xmax10e6 = [[1000000], [225, 353, 618]];
    var list = new List<int[][]>() { xmax10e4, xmax10e5, xmax10e6 };
    var nmax = 1000;
    foreach (int n in nmax.Range(1))
    {
        var xmin = (int)double.Ceiling(double.Pow(n, 1.0 / 3.0));
        var idx = list.FindIndex(e => e[1].Contains(n));
        var xmax = idx != -1 ? list[idx][0][0] : 1000;
        for (var x = -xmin - 1; x < xmax; x++)
        {
            var (y, info) = ApproxSolver(x, n);
            if (!info)
                continue;
            
            set.UnionWith([(x, y, n), (x, -y, n)]);
        }

        Console.WriteLine($"n = {n} ");
        Console.CursorTop--;
    }
    
    var missingSet = new List<int>();
    foreach (var s in set.GroupBy(e => e.Item3))
    {
        var n = s.Key;
        var sols = s.ToArray();
        var nb = sols.Length;
        var pts = sols.Select(e => (e.Item1, e.Item2)).Order().ToArray();
        var missing = B081119[n] - nb;
        Console.WriteLine($"n = {n,4}, Integral Points {nb,-3}/{B081119[n],3} => {{ {pts.Glue(", ")} }}");
        if (missing != 0)
            missingSet.Add(n);
    }

    if (missingSet.Count != 0)
        Console.WriteLine($"{missingSet.Count} Missing {{ {missingSet.Glue(", ")} }}");
    
    Console.WriteLine();
    GlobalStopWatch.Show(); // Time:6.191s
}
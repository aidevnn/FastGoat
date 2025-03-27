using System.Numerics;
using System.Reflection;
using FastGoat.Commons;
using FastGoat.Examples;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.EllCurve;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.LWE;
using FastGoat.UserGroup.Matrix;
using FastGoat.UserGroup.Padic;
using FastGoat.UserGroup.Polynoms;
using Rq = FastGoat.Structures.VecSpace.KPoly<FastGoat.UserGroup.Integers.Rational>;
using static FastGoat.Commons.IntExt;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");
Ring.DisplayPolynomial = MonomDisplay.StarCaret;

{
    GlobalStopWatch.Restart();
    var nMax = 20000;
    var n = 31; // 5, 34, 1254 
    var rk = EllipticCurves.EllRank(-n * n, 0);
    var ng = EllipticCurves.NagellLutz(-n * n, 0);
    ng.intPts.Println("Int Pts");
    
    var xs = ng.E.Where(e => !e.IsO).Select(e => e.X.Num)
        .Union(ng.intPts.Select(e => e.X.Num)).Where(e => e > 0).ToHashSet();
    Ring.DisplayPolynomial = MonomDisplay.StarCaret;

    // congruence problem for n
    // a/b^2 solution of Y^2=X^3 - n^2*X
    // (a^3 - n^2*b^4 * a)/b^6 = (c/b^3)^2
    // a^3 - n^2*b^4 * a = c^2
    // new congruence for n1 = (n*b^2)

    GlobalStopWatch.AddLap();
    foreach (var b in nMax.SeqLazy(1))
    {
        var a = BigInteger.Pow(n * b * b, 2);
        var set = nMax.SeqLazy(1).Select(x => (x: new BigInteger(x), y2: BigInteger.Pow(x, 3) - a * x))
            .Where(e => !xs.Contains(e.x) && e.y2 > 0)
            .Select(e => (e.x, e.y2, y: SqrtBigInt(e.y2)))
            .Where(e => BigInteger.Abs(e.y2 - e.y * e.y) == 0)
            .Take(rk)
            .ToArray();

        if (set.Any(e => e.y != 0))
        {
            Console.WriteLine($"b={b} a={a}");
            set.Println(e => $"({e.x}, {e.y}) => ({new Rational(e.x, b * b)}, {new Rational(e.y, b * b * b)})", "Pts");
            break;
        }
    }

    GlobalStopWatch.Show(); // Time:1.134s
}
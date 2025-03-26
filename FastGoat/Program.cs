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

// RENÉ SCHOOF
// Counting points on elliptic curves over finite fields
// Journal de Théorie des Nombres de Bordeaux, tome 7, no 1 (1995),
// p. 219-254
int SchoofEllPtsCount(BigInteger a, BigInteger b, int p)
{
    return p + 1 + p.Range().Select(x => LegendreJacobiBigint((BigInteger.ModPow(x, 3, p) + a * x + b) % p, p))
        .Sum(k => k <= 1 ? (int)k : -1);
}

// BULLETIN (New Series) OF THE
// AMERICAN MATHEMATICAL SOCIETY
// Volume 39, Number 4, Pages 455–474
// S 0273-0979(02)00952-7
// Article electronically published on July 8, 2002
// RANKS OF ELLIPTIC CURVES
// KARL RUBIN AND ALICE SILVERBERG
//
// RANKS “CHEAT SHEET”
// ALICE SILVERBERG
int EllRank(BigInteger a, BigInteger b, int n = 500)
{
    GlobalStopWatch.AddLap();
    var r = 1.0;
    var (sumX, sumY, sumX2, sumXY) = (0.0, 0.0, 0.0, 0.0);
    foreach (var p in Primes10000.Take(n))
    {
        r *= 1.0 * SchoofEllPtsCount(a, b, p) / p;

        var (x, y) = (double.Log(double.Log(p)), double.Log(r));
        sumX += x;
        sumY += y;
        sumX2 += x * x;
        sumXY += x * y;
    }

    var B = (sumY * sumX2 - sumX * sumXY) / (n * sumX2 - sumX * sumX);
    var A = (n * sumXY - sumX * sumY) / (n * sumX2 - sumX * sumX);
    var rk = (int)double.Round(A);
    Console.WriteLine($"Y = {A:f4} * X + {B:f4}");
    Console.WriteLine($"Rank(E[{a},{b}](Q)) = {rk}");
    GlobalStopWatch.Show();
    Console.WriteLine();

    return rk;
}

{
    GlobalStopWatch.Restart();

    // Rank 0
    EllRank(-432, 8208);
    EllRank(-675, 13662);
    EllRank(-27, 8694);

    // Rank 1
    EllRank(0, 3);
    EllRank(-36, 0);
    EllRank(-961, 0);

    // Rank 2
    EllRank(-3024, 46224);
    EllRank(-5292, -101520);
    EllRank(-2052, 34560);

    // Rank 0, 1, 2, 3, 4
    foreach (var d in new[] { 1, 5, 34, 1254, 29274 })
        EllRank(-d * d, 0);

    GlobalStopWatch.Show(); // Time:5.079s
}
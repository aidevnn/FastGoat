using System.Numerics;
using FastGoat.Commons;
using FastGoat.Examples;
using FastGoat.Structures;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.EllCurve;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Polynoms;
using static FastGoat.Commons.IntExt;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");
Ring.DisplayPolynomial = MonomDisplay.StarCaret;

void factorFq()
{
    var q = BigInteger.Pow(5, 24);
    var (X, a) = FG.EPolyXc(FG.FqX(q, 'a').F, 'a');
    var P = X.Pow(24) + 4 * X.Pow(22) + 3 * X.Pow(20) + 4 * X.Pow(18) + 2 * X.Pow(16) + 3 * X.Pow(14) + 2 * X.Pow(12) +
            X.Pow(10) + 4 * X.Pow(8) + 3 * X.Pow(6) + X.Pow(4) + 3 * X.Pow(2) + 2;

    GlobalStopWatch.AddLap();
    foreach (var f in IntFactorisation.CantorZassenhausAECF(P, a, q))
        Console.WriteLine(f);

    GlobalStopWatch.Show("CantorZassenhausAECF");
    Console.WriteLine();

    GlobalStopWatch.AddLap();
    foreach (var f in IntFactorisation.BerlekampProbabilisticAECF(P, a, q))
        Console.WriteLine(f);

    GlobalStopWatch.Show("BerlekampProbabilisticAECF");
    Console.WriteLine();
}

void benchFactorsFq()
{
    RngSeed(1452);
    var (p, n1, n2) = (5, 5, 10);
    var q1 = BigInteger.Pow(p, n1);
    var q2 = BigInteger.Pow(p, n2);
    var a1 = FG.FqX(q1, 'a');
    var a2 = FG.FqX(q2, 'a');

    var X1 = FG.KPoly('X', a1);
    var P1 = n1.SeqLazy().Select(_ => X1 - Ring.FastPow(a1, Rng.NextInt64((long)q1))).Aggregate((ai, aj) => ai * aj);
    Console.WriteLine($"Factors P = {P1}");
    Console.WriteLine($"    in F{p}^{n1} and F{p}^{n2}");
    var facts1 = IntFactorisation.CantorZassenhausAECF(P1, a1, q1).ToArray();
    var cv = facts1.ToDictionary(e => e, e => e.ToGF(q2));
    cv.Println(l => $"{l.Key} --> {l.Value}", "Correspondance");
    var facts2 = IntFactorisation.CantorZassenhausAECF(P1.ToGF(q2), a2, q2).ToArray();
    Console.WriteLine($"Check {facts2.ToHashSet().SetEquals(cv.Values)}");
    Console.WriteLine();

    var P2 = P1.ToGF(q2);
    GlobalStopWatch.Bench(5, "BK-AECF  ", () => IntFactorisation.BerlekampProbabilisticAECF(P2, a2, q2).Count());
    GlobalStopWatch.Bench(5, "BK-VShoup", () => IntFactorisation.BerlekampProbabilisticVShoup(P2, a2, q2).Count());
    GlobalStopWatch.Bench(5, "CZ-AECF  ", () => IntFactorisation.CantorZassenhausAECF(P2, a2, q2).Count());
    GlobalStopWatch.Bench(5, "CZ-VShoup", () => IntFactorisation.CantorZassenhausVShoup(P2, a2, q2).Count());
    GlobalStopWatch.Bench(5, "BK-AECF  ", () => IntFactorisation.BerlekampProbabilisticAECF(P2, a2, q2).Count());
    GlobalStopWatch.Bench(5, "BK-VShoup", () => IntFactorisation.BerlekampProbabilisticVShoup(P2, a2, q2).Count());
    GlobalStopWatch.Bench(5, "CZ-AECF  ", () => IntFactorisation.CantorZassenhausAECF(P2, a2, q2).Count());
    GlobalStopWatch.Bench(5, "CZ-VShoup", () => IntFactorisation.CantorZassenhausVShoup(P2, a2, q2).Count());

    GlobalStopWatch.Bench(10, "BK-AECF  ", () => IntFactorisation.BerlekampProbabilisticAECF(P2, a2, q2).Count());
    GlobalStopWatch.Bench(10, "BK-VShoup", () => IntFactorisation.BerlekampProbabilisticVShoup(P2, a2, q2).Count());
    GlobalStopWatch.Bench(10, "CZ-AECF  ", () => IntFactorisation.CantorZassenhausAECF(P2, a2, q2).Count());
    GlobalStopWatch.Bench(10, "CZ-VShoup", () => IntFactorisation.CantorZassenhausVShoup(P2, a2, q2).Count());
}

{
    GlobalStopWatch.Restart();
    
    factorFq();
    benchFactorsFq();
}
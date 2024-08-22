using System.CodeDom;
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
using FastGoat.UserGroup.DatabaseSmallGroups;
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

void PslqExprPi(int algo = 0)
{
    var n = 8;
    var O1 = 30;
    var pi = BigReal.Pi(2 * n + O1);
    // beta = 5/8*π^6 - 35/6*π^4 - 5*π^2
    //      = −16.69947371922907049618724340073146784130179174288144470245664281170485208378578722
    var beta = BigReal.FromBigIntegerAndExponent(BigInteger.Parse("-1669947371922907049618724340073146784"), 1, O1);

    var ai = n.Range().Select(k => pi.Pow(k).ToBigReal(O1)).ToKMatrix();
    ai.Coefs[0, n - 1] = beta;

    var gamma = 2 / BigReal.Sqrt(BigReal.FromBigInteger(3, O1));
    
    GlobalStopWatch.AddLap();
    var (coefs, name) = algo == 0
        ? (PSLQ.OnelevelMultipair(ai, gamma), "One level")
        : algo == 1
            ? (PSLQ.TwoLevelMultipairXP(ai, gamma), "Two level bigreal/bigreal")
            : (PSLQM2.TwoLevelMultipairXP(ai, gamma), "Two level decimal/bigreal");

    Console.WriteLine(coefs.ToKMatrix());
    GlobalStopWatch.Show($"{name} Multipair PSLQ");
    var P = FG.KPoly('π', coefs.SkipLast(1).ToArray());
    Console.WriteLine($"beta = {-P / coefs.Last()}");
    Console.WriteLine();
}

KPoly<Rational> LLLminPoly(int r, int s, int O)
{
    var n = r * s + 1; // Expected polynomial degree plus one 
    var alpha = BigReal.NthRoot(3, r, O) - BigReal.NthRoot(2, s, O); // a = 3^(1/r) - 2^(1/s)

    if (Logger.Level != LogLevel.Off)
        GlobalStopWatch.AddLap();
    
    var coefs = AlgebraicIntegerRelationLLL.AlphaBetaPolynomial(alpha, alpha.Pow(n - 1), n, O);
    var x = FG.QPoly('X');
    var P = x.Pow(n - 1) - coefs.Select((c, k) => c * x.Pow(k)).Aggregate((a, b) => a + b);
    if (Logger.Level != LogLevel.Off)
    {
        GlobalStopWatch.Show($"LLL min poly a = 3^(1/{r}) - 2^(1/{s})");
        Console.WriteLine($"P = {P} P(a) = {P.Substitute(alpha).ToBigReal(3 * O / 4)}");
        Console.WriteLine();
    }
    return P;
}

KPoly<Rational> PslqMinPoly(int r, int s, int O, int algo = 0)
{
    var n = r * s + 1; // Expected polynomial degree plus one
    var alpha = BigReal.NthRoot(3, r, O) - BigReal.NthRoot(2, s, O); // a = 3^(1/r) - 2^(1/s)
    var ai = n.Range().Select(k => alpha.Pow(k)).ToKMatrix();

    var gamma = 2 / BigReal.Sqrt(BigReal.FromBigInteger(3, O));
    
    GlobalStopWatch.AddLap();
    var (coefs, name) = algo == 0
        ? (PSLQ.OnelevelMultipair(ai, gamma), "One level")
        : algo == 1
            ? (PSLQ.TwoLevelMultipairXP(ai, gamma), "Two level bigreal/bigreal")
            : (PSLQM2.TwoLevelMultipairXP(ai, gamma), "Two level decimal/bigreal");

    Console.WriteLine(coefs.ToKMatrix());
    var P = FG.KPoly('X', coefs).Monic;
    GlobalStopWatch.Show($"{name} Multipair PSLQ min poly a = 3^(1/{r}) - 2^(1/{s})");
    if (Logger.Level != LogLevel.Off)
    {
        Console.WriteLine($"P = {P} and P(a) = {P.Substitute(alpha).ToBigReal(3 * O / 4)}");
        Console.WriteLine();
    }

    return P;
}

void RunMinPolyWithParams((int r, int s, int O)[] rsO, bool lll, params int[] algos)
{
    var x = FG.QPoly('X');
    KPoly<Rational> P1 = x.Zero, P2;
    foreach (var (r, s, O) in rsO)
    {
        Console.WriteLine($"Min poly of a = 3^(1/{r}) - 2^(1/{s}), A.I.R. algos on {O} digits");
        Console.WriteLine();

        if (lll)
            P1 = LLLminPoly(r, s, O);
        foreach (var algo in algos.Where(k => k >= 0 && k <= 2))
        {
            P2 = PslqMinPoly(r, s, O, algo);
            if (lll && !P1.Equals(P2))
                throw new();
        }
    }
}

void RunExprPi()
{
    GlobalStopWatch.Restart();
    Ring.DisplayPolynomial = MonomDisplay.StarCaret;
    Logger.Level = LogLevel.Level1;
    
    PslqExprPi();
    PslqExprPi(algo: 1);
    PslqExprPi(algo: 2);
}

void RunMinPoly()
{
    GlobalStopWatch.Restart();
    Ring.DisplayPolynomial = MonomDisplay.StarCaret;
    Logger.Level = LogLevel.Level1;
    
    var rsO = new[]
    {
        (2, 2, 20), (3, 3, 50), (2, 5, 70), (3, 4, 70), (2, 7, 90),
        (3, 5, 90), (4, 4, 90), (4, 5, 120), (5, 5, 180)
        
    };

    // LLL, PSLQM1 and PSLQM2
    // Run(rsO.ToArray(), lll: true, algos: [0, 1]);

    // PSLQM2 bigreal and PSLQM2 decimal
    RunMinPolyWithParams(rsO.ToArray(), lll: false, algos: [0, 1, 2]);
}

{
    RunExprPi();
    RunMinPoly();
    Console.Beep();
}

// Possible Solution step:57
// [0, 0, 120, 0, 140, 0, -15, 24]
// # One level Multipair PSLQ Time:184ms
// beta = 5/8*π^6 - 35/6*π^4 - 5*π^2
// 
// Possible Solution step:60
// [0, 0, 120, 0, 140, 0, -15, 24]
// # Two level bigreal/bigreal Multipair PSLQ Time:93ms
// beta = 5/8*π^6 - 35/6*π^4 - 5*π^2
// 
// Possible Solution step:60
// [0, 0, 120, 0, 140, 0, -15, 24]
// # Two level decimal/bigreal Multipair PSLQ Time:41ms
// beta = 5/8*π^6 - 35/6*π^4 - 5*π^2
// 
// ...
//
// Possible Solution step:558
// [1, 0, 0, 0, 0, -116255, 0, 0, 0, 0, -11240, 0, 0, 0, 0, -3760, 0, 0, 0, 0, 5, 0, 0, 0, 0, -1]
// # One level Multipair PSLQ min poly a = 3^(1/5) - 2^(1/5) Time:2m8s
// P = X^25 - 5*X^20 + 3760*X^15 + 11240*X^10 + 116255*X^5 - 1 and P(a) = 0
// 
// Possible Solution step:598
// [1, 0, 0, 0, 0, -116255, 0, 0, 0, 0, -11240, 0, 0, 0, 0, -3760, 0, 0, 0, 0, 5, 0, 0, 0, 0, -1]
// # Two level bigreal/bigreal Multipair PSLQ min poly a = 3^(1/5) - 2^(1/5) Time:16.857s
// P = X^25 - 5*X^20 + 3760*X^15 + 11240*X^10 + 116255*X^5 - 1 and P(a) = 0
// 
// Possible Solution step:581
// [1, 0, 0, 0, 0, -116255, 0, 0, 0, 0, -11240, 0, 0, 0, 0, -3760, 0, 0, 0, 0, 5, 0, 0, 0, 0, -1]
// # Two level decimal/bigreal Multipair PSLQ min poly a = 3^(1/5) - 2^(1/5) Time:5.225s
// P = X^25 - 5*X^20 + 3760*X^15 + 11240*X^10 + 116255*X^5 - 1 and P(a) = 0
// 

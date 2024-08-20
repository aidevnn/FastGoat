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

void testPSQL1()
{
    var n = 8;
    var O1 = 20;
    var O2 = O1 + n;
    var pi = BigReal.Pi(O2);
    var beta = BigReal.FromBigIntegerAndExponent(BigInteger.Parse("-1669947371922907049619"), 1, O2);

    var ai = n.Range().Select(k => pi.Pow(k)).ToKMatrix();
    ai.Coefs[0, n - 1] = beta;

    var gamma = 2 / BigReal.Sqrt(BigReal.FromBigInteger(3, O2));
    // var gamma = BigReal.FromBigInteger(3, O2) / 2;

    GlobalStopWatch.AddLap();
    var coefs = PSLQ.OnelevelMultipair(ai, gamma);
    Console.WriteLine(coefs.ToKMatrix());
    GlobalStopWatch.Show($"multipair PSLQ");
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

KPoly<Rational> PSLQonelvlminPoly(int r, int s, int O)
{
    var n = r * s + 1; // Expected polynomial degree plus one
    var alpha = BigReal.NthRoot(3, r, O) - BigReal.NthRoot(2, s, O); // a = 3^(1/r) - 2^(1/s)
    var ai = n.Range().Select(k => alpha.Pow(k)).ToKMatrix();

    var gamma = 2 / BigReal.Sqrt(BigReal.FromBigInteger(3, O));
    // var y = BigReal.FromBigInteger(3, O2) / 2;
    // var y = BigReal.FromBigInteger(2, O2);

    GlobalStopWatch.AddLap();
    var coefs = PSLQ.OnelevelMultipair(ai, gamma);
    var P = FG.KPoly('X', coefs).Monic;
    GlobalStopWatch.Show($"One level multipair PSLQ min poly a = 3^(1/{r}) - 2^(1/{s})");
    if (Logger.Level != LogLevel.Off)
    {
        Console.WriteLine($"P = {P} and P(a) = {P.Substitute(alpha).ToBigReal(3 * O / 4)}");
        Console.WriteLine();
    }

    return P;
}

void Run()
{
    GlobalStopWatch.Restart();
    Ring.DisplayPolynomial = MonomDisplay.StarCaret;

    testPSQL1();
    testPSQL1();
    
    // MinPoly of a = 3^(1/r) - 2^(1/s)
    List<(int r, int s, int O)> rsO = new()
    {
        (2, 2, 20),
        (3, 3, 40),
        (2, 5, 50),
        (3, 4, 70),
        (2, 7, 90),
        (3, 5, 90),
        (4, 4, 90),
        (4, 5, 120),
        (5, 5, 180),
    };

    Logger.Level = LogLevel.Level1;
    foreach (var (r, s, O) in rsO.SkipLast(4))
    {
        Console.WriteLine($"Min poly of a = 3^(1/{r}) - 2^(1/{s}) with LLL and multipair PSLQ algorithms on {O} digits");
        Console.WriteLine();

        var P1 = LLLminPoly(r, s, O);
        var P2 = PSLQonelvlminPoly(r, s, O);
        
        if (!P1.Equals(P2))
            throw new();
    }
    
    // Min poly of a = 3^(1/5) - 2^(1/5) with LLL and multipair PSLQ algorithms on 180 digits
    // 
    // # LLL  min poly a = 3^(1/5) - 2^(1/5) Time:1m19s
    // P = X^25 - 5*X^20 + 3760*X^15 + 11240*X^10 + 116255*X^5 - 1 P(a) = 0
    // 
    // # PSLQ multipair min poly a = 3^(1/5) - 2^(1/5) Time:1m46s
    // P = X^25 - 5*X^20 + 3760*X^15 + 11240*X^10 + 116255*X^5 - 1 and P(a) = 0
    // 
}

void QRtest()
{
    var o = BigReal.BrOne(20);
    var A = new[] { 12, -51, 4, 6, 167, -68, -4, 24, -41 }.Select(e => e * o).ToKMatrix(3);
    Console.WriteLine("A");
    Console.WriteLine(A);
    Console.WriteLine();
    
    var (Q, R) = PSLQ.QR(A);
    Console.WriteLine("Q");
    Console.WriteLine(Q);
    Console.WriteLine();
    Console.WriteLine("R");
    Console.WriteLine(R.Select(c => BigReal.Round(c, 4)).ToKMatrix(3));
    Console.WriteLine();

    Console.WriteLine("QRpslq");
    var B = PSLQ.LQpslq(A.T).T;
    Console.WriteLine(B.Select(c => BigReal.Round(c, 4)).ToKMatrix(3));
    Console.WriteLine();
    var C = PSLQ.LQpslq(A.Select(c => c.ToDble).ToKMatrix(3).T).T;
    Console.WriteLine(C.Select(c => Dble.Round(c, 4)).ToKMatrix(3));

    for (int i = 0; i < 20; i++)
    {
        var A0 = 12.Range().Select(_ => Rng.Next(-20, 21) * o).ToKMatrix(3);
        if (Rng.Next(5) == 0)
            A0.Coefs[0, 0] = A0.KZero;
        
        Console.WriteLine("A0");
        Console.WriteLine(A0);
        Console.WriteLine();
        var (_, R0) = PSLQ.QR(A0);
        Console.WriteLine("R0");
        Console.WriteLine(R0.Select(c => BigReal.Round(c, 4)).ToKMatrix(3));
        Console.WriteLine();
        
        var R1 = PSLQ.LQpslq(A0.T).T;
        Console.WriteLine("R1");
        Console.WriteLine(R1.Select(c => BigReal.Round(c, 4)).ToKMatrix(3));
        Console.WriteLine();

        var check = (R0 - R1).Aggregate(o.Zero, (acc, c) => acc + c * c).IsZero();
        Console.WriteLine($"Are Equal {check}");
        Console.WriteLine();
        if (!check)
            throw new();
    }
}

KPoly<Rational> PSLQtwolvlminPoly(int r, int s, int O)
{
    var n = r * s + 1; // Expected polynomial degree plus one
    var alpha = BigReal.NthRoot(3, r, O) - BigReal.NthRoot(2, s, O); // a = 3^(1/r) - 2^(1/s)
    var ai = n.Range().Select(k => alpha.Pow(k)).ToKMatrix();
    
    var gamma = 2 / BigReal.Sqrt(BigReal.FromBigInteger(3, O));
    
    GlobalStopWatch.AddLap();
    var coefs = PSLQ.TwoLevelMultipair(ai, gamma);
    Console.WriteLine(coefs.Select(c => c.RoundEven).ToKMatrix());
    var P = FG.KPoly('X', coefs).Monic;
    GlobalStopWatch.Show($"Two level multipair PSLQ min poly a = 3^(1/{r}) - 2^(1/{s})");
    if (Logger.Level != LogLevel.Off)
    {
        Console.WriteLine($"P = {P} and P(a) = {P.Substitute(alpha).ToBigReal(3 * O / 4)}");
        Console.WriteLine();
    }
    Console.WriteLine();
    return P;
}

{
    GlobalStopWatch.Restart();
    Ring.DisplayPolynomial = MonomDisplay.StarCaret;
    Logger.Level = LogLevel.Level1;
    
    testPSQL1();
    testPSQL1();
    testPSQL1();
    
    // MinPoly of a = 3^(1/r) - 2^(1/s)
    List<(int r, int s, int O)> rsO = new()
    {
        (2, 2, 20),
        (3, 3, 50),
        (2, 5, 60),
        (3, 4, 70),
        (2, 7, 90),
        (3, 5, 90),
        (4, 4, 90),
        (4, 5, 120),
        (5, 5, 180),
    };

    Logger.Level = LogLevel.Level1;
    foreach (var (r, s, O) in rsO)
    {
        Console.WriteLine(
            $"Min poly of a = 3^(1/{r}) - 2^(1/{s}) with LLL and two level multipair PSLQ algorithms on {O} digits");
        Console.WriteLine();
        
        var P1 = LLLminPoly(r, s, O);
        var P2 = PSLQtwolvlminPoly(r, s, O);
        
        if (!P1.Equals(P2))
            throw new();
    }

    Console.Beep();
}
// Possible Solution step:230
// [-1, 0, 0, 0, 3860, 0, 0, 0, 666, 0, 0, 0, 20, 0, 0, 0, -1]
// 
// # One level multipair PSLQ min poly a = 3^(1/4) - 2^(1/4) Time:10.477s
// P = X^16 - 20*X^12 - 666*X^8 - 3860*X^4 + 1 and P(a) = 0
// 
// Possible Solution step:273
// [1, 0, 0, 0, -3860, 0, 0, 0, -666, 0, 0, 0, -20, 0, 0, 0, 1]
// # Two level multipair PSLQ Time:2.815s
// 
// ...
// 
// Possible Solution
// [-1, 0, 0, 0, 0, 116255, 0, 0, 0, 0, 11240, 0, 0, 0, 0, 3760, 0, 0, 0, 0, -5, 0, 0, 0, 0, -12233]
// 
// # LLL min poly a = 3^(1/5) - 2^(1/5) Time:1m12s
// P = X^25 - 5*X^20 + 3760*X^15 + 11240*X^10 + 116255*X^5 - 1 P(a) = 0
// 
// Possible Solution step:586
// [-1, 0, 0, 0, 0, 116255, 0, 0, 0, 0, 11240, 0, 0, 0, 0, 3760, 0, 0, 0, 0, -5, 0, 0, 0, 0, 1]
// # Two level multipair PSLQ min poly a = 3^(1/5) - 2^(1/5) Time:17.084s
// P = X^25 - 5*X^20 + 3760*X^15 + 11240*X^10 + 116255*X^5 - 1 and P(a) = 0
// 
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

void PslqExprPi()
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
    Console.WriteLine(ai);
    
    GlobalStopWatch.AddLap();
    var coefs = PSLQM2.TwoLevelMultipair(ai, gamma);
    Console.WriteLine(coefs.ToKMatrix());
    GlobalStopWatch.Show("Two level Multipair PSLQ");
    var P = FG.KPoly('π', coefs.SkipLast(1).ToArray());
    Console.WriteLine($"beta = {-P / coefs.Last()}");
    Console.WriteLine();
}

KPoly<Rational> PslqMinPoly(int r, int s, int O)
{
    var n = r * s + 1; // Expected polynomial degree plus one
    var alpha = BigReal.NthRoot(3, r, O) - BigReal.NthRoot(2, s, O); // a = 3^(1/r) - 2^(1/s)
    var ai = n.Range().Select(k => alpha.Pow(k)).ToKMatrix();

    var gamma = 2 / BigReal.Sqrt(BigReal.FromBigInteger(3, O));
    
    GlobalStopWatch.AddLap();
    var coefs = PSLQM2.TwoLevelMultipair(ai, gamma);
    Console.WriteLine(coefs.ToKMatrix());
    var P = FG.KPoly('X', coefs).Monic;
    GlobalStopWatch.Show($"Two level Multipair PSLQ min poly a = 3^(1/{r}) - 2^(1/{s})");
    Console.WriteLine($"P = {P} and P(a) = {P.Substitute(alpha).ToBigReal(3 * O / 4)}");
    Console.WriteLine();

    return P;
}

{
    Ring.DisplayPolynomial = MonomDisplay.StarCaret;
    PslqExprPi();
    PslqMinPoly(2, 2, 30);
    PslqMinPoly(3, 3, 50);
    PslqMinPoly(4, 4, 90);
    PslqMinPoly(5, 5, 180);
    PslqMinPoly(5, 6, 250);
    Console.Beep();
}
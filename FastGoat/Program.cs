using System.Numerics;
using FastGoat.Commons;
using FastGoat.Examples;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.LWE;
using FastGoat.UserGroup.Polynoms;
using Rq = FastGoat.Structures.VecSpace.KPoly<FastGoat.UserGroup.Integers.Rational>;
using static FastGoat.Commons.IntExt;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");
// RngSeed(259128);
Ring.DisplayPolynomial = MonomDisplay.StarCaret;
RecomputeAllPrimesUpTo(5000000);

{
    for (int i = 2; i < 12; i++)
    {
        var n = 1 << i;
        var a = RLWE.Alpha(n);
        var p = RLWE.RlwePrime(n);
        var s = RLWE.Sigma(n, p);
        Console.WriteLine($"N={2 * n} n={n} p={p} s={s:f4} a={a:f4}");
    }
}
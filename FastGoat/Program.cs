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
RecomputeAllPrimesUpTo(200000);

{
    for (int i = 2; i < 12; i++)
    {
        var n = 1 << i;
        var N = 2 * n;
        var a = RLWE.Alpha(n);
        var p = RLWE.RlwePrime(n);
        var q = RLWE.FirstPrimeEqualOneMod(N * p);
        var s = RLWE.Sigma(n, p);
        Console.WriteLine($"N={2 * n} n={n} p={p} q={q} s={s:f4} a={a:f4}");
        Console.WriteLine($"        q={q.Mod(N)} mod N");
        Console.WriteLine($"         ={q.Mod(p)} mod p");
    }
}
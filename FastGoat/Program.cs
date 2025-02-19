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
    var seq1 = 200.SeqLazy().Where(p => IsPrime(p)).ToArray();
    var seq2 = Primes10000.Where(p => p < 200).ToArray();
    Console.WriteLine("seq1 = seq2 {0}", seq1.SequenceEqual(seq2));

    foreach (var e in 2000.SeqLazy(2))
    {
        var dec = PrimesDecomposition(e).ToArray();
        if (e != dec.Aggregate((ei, ej) => ei * ej) || dec.Any(ei => !Primes10000.Contains(ei)))
            throw new();
    }

    Console.WriteLine("Decomposition in primes factors of integer up to 2000");

    Console.WriteLine($"26172497720 = {PrimesDecompositionBigInt(26172497720).Glue(" * ")}");
    Console.WriteLine($"26172497920 = {PrimesDecompositionBigInt(26172497920).Glue(" * ")}");
    try
    {
        Console.WriteLine($"26172497921 = {PrimesDecompositionBigInt(26172497921).Glue(" * ")}");
    }
    catch (Exception e)
    {
        Console.WriteLine(e.Message);
    }

    try
    {
        Console.WriteLine($"52076511233 = {PrimesDecompositionBigInt(52076511233).Glue(" * ")}");
    }
    catch (Exception e)
    {
        Console.WriteLine(e.Message);
    }

    Console.WriteLine();
}

{
    for (int i = 2; i <= 12; i++)
    {
        var n = 1 << i;
        var a = RLWE.Alpha(n);
        var p = RLWE.RlwePrime(n);
        var q = RLWE.FirstPrimeEqualOneMod(2 * n * p);
        var s = RLWE.Sigma(n, p);
        Console.WriteLine($"n={n,-6} p={p,-6} q={q,-12} s={s,-8:f4} a={a:f4} log2(p)={BigInteger.Log2(p.Num)}");
        Console.WriteLine($"         p={p.Mod(2 * n)}      q={q.Mod(2 * n)} mod 2*n");
        Console.WriteLine($"                  q={q.Mod(p)} mod p");
        Console.WriteLine();
    }
}
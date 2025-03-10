using System.Numerics;
using FastGoat.Commons;
using FastGoat.Examples;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.LWE;
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
// RngSeed(259128);
Ring.DisplayPolynomial = MonomDisplay.StarCaret;
RecomputeAllPrimesUpTo(1000000);

int NthRootsANTV1Composite(int a, int r, int p)
{
    var ai = AmodP(a, p);
    if (r == 1 || ai <= 1)
        return ai;

    return PrimesDecomposition(r).Aggregate(ai, (acc, ri) => NthRootsANTV1Step(acc, ri, p));
}

int NthRootsANTV1Step(int a, int r, int p)
{
    if (Gcd(r, p - 1) == 1)
    {
        var invr = InvModPbez(r, p - 1);
        return PowMod(a, invr, p);
    }

    return NumberTheory.NthRootANTV1(a, r, p);
}

void testNthRoots()
{
    var p = 97;
    for (int r = 2; r < p / 2; r++)
    {
        for (int a = 0; a < p; a++)
        {
            try
            {
                var g = NthRootsANTV1Composite(a, r, p);
                var gr = PowMod(g, r, p);
                Console.WriteLine($"p={p} r={r,3} a={a,3} g={g,3} g^r={gr,3} check:{gr == a}");
            }
            catch (Exception e)
            {
                Console.WriteLine($"p={p} r={r,3} a={a,3} error msg:{e.Message}");
            }
        }

        Console.WriteLine();
    }
}

void test_CT_GS()
{
    // RLWE.NoiseOff();
    var N = 64;
    var n = N / 2;
    var level = 4;
    var t = RLWE.CiphertextModulusBGV(n);
    var primes = RLWE.SequencePrimesBGV(N, t, level).primes;
    var q = primes[0];
    var qL = primes.Aggregate((pi, pj) => pi * pj);
    var tables = RLWE.PrepareNTT(n, t, primes);
    Console.WriteLine($"N={N} n={n} t={t} primes=[{primes.Glue(", ")}] qL={qL}");

    for (int k = 0; k < 20; ++k)
    {
        var m = RLWE.GenUnif(n, qL);
        var c1 = RLWE.Rq2NTT(m, tables);
        var c2 = RLWE.CooleyTukey(m, tables);
        var mf = RLWE.GentlemanSande(c1, tables);
        if (!m.Equals(mf))
        {
            Console.WriteLine(m);
            Console.WriteLine(mf);
            throw new();
        }
        if (!c1.Equals(c2))
        {
            Console.WriteLine(c1.T);
            Console.WriteLine(c2.T);
            throw new();
        }
    }

    Console.WriteLine("Success");
}

{
    GlobalStopWatch.Restart();
    var lt = 100.SeqLazy().ToList();
    for (int l = 0; l < 4; l++)
    {
        var N = 1 << (l + 4);
        var n = N / 2;
        var level = 5;
        var t = RLWE.CiphertextModulusBGV(n);
        var primes = RLWE.SequencePrimesBGV(N, t, level).primes;
        var q = primes[0];
        var qL = primes.Aggregate((pi, pj) => pi * pj);
        var tables = RLWE.PrepareNTT(n, t, primes);
        Console.WriteLine($"N={N} n={n} t={t} primes=[{primes.Glue(", ")}] qL={qL}");

        var m = RLWE.GenUnif(n, qL);
        var c1 = RLWE.Rq2NTT(m, tables);
        GlobalStopWatch.Bench(10, "MatMul NTT ", () => lt.ForEach(_ => RLWE.Rq2NTT(m, tables)));
        GlobalStopWatch.Bench(10, "CT         ", () => lt.ForEach(_ => RLWE.CooleyTukey(m, tables)));
        GlobalStopWatch.Bench(10, "MatMul INTT", () => lt.ForEach(_ => RLWE.NTT2Rq(c1, tables)));
        GlobalStopWatch.Bench(10, "GS         ", () => lt.ForEach(_ => RLWE.GentlemanSande(c1, tables)));
        Console.WriteLine();
    }
}
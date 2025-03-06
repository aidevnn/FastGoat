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
RngSeed(259128);
Ring.DisplayPolynomial = MonomDisplay.StarCaret;
RecomputeAllPrimesUpTo(1000000);

void testRgswNTT()
{
    // RLWE.NoiseOff();
    var N = 16;
    var n = N / 2;
    var level = 3;
    var (pm, sk, t, primes, sp, pk, _) = RLWE.SetupBGV(N, level);
    var q = primes[0];
    var qL = pk.Q;
    var B = RLWE.RNSGadgetBase(primes);
    Console.WriteLine($"BGV level = {level}, Primes = [{primes.Glue(", ")}]");
    Console.WriteLine($"pm = {pm} T = {t} q = {q} sp = {sp} qL = {qL} sigma = {RLWE.Sigma(n, t):f4}");
    Console.WriteLine($"sk = {sk}");
    Console.WriteLine($"pk => {pk.Params}");
    Console.WriteLine();

    var tables = RLWE.PrepareNTT(n, t, primes);

    for (int k = 0; k < 5; ++k)
    {
        var m1 = RLWE.GenUnif(n, t);
        var m2 = RLWE.GenUnif(n, t);
        var cm1 = RLWE.EncryptBGV(m1, pk);
        var cm1ntt = RLWE.DecompRNS(cm1, primes).Select(e => RLWE.RLWECipher2NTT(e, tables)).ToVec();
        var (csm2ntt, cm2ntt) = RLWE.EncryptRgswNTTBGV(m2, pk, B, tables);

        var cm1m2ntt = RLWE.MulRgsw(cm1ntt, cm2ntt, csm2ntt);
        var d_m1m2ntt = RLWE.DecryptBGV(RLWE.NTTCipher2RLWE(cm1m2ntt, tables), sk);
        Console.WriteLine($"m1        = {m1}");
        Console.WriteLine($"m2        = {m2}");
        Console.WriteLine($"m1 * m2   = {(m1 * m2).ResModSigned(pm, t)}");
        Console.WriteLine($"          = {d_m1m2ntt}");
        Console.WriteLine($"          = {RLWE.ErrorsBGV(RLWE.NTTCipher2RLWE(cm1m2ntt, tables), sk).NormInf()}");
        Console.WriteLine();

        var cm1m2ntt_ = RLWE.DecompRNS(cm1m2ntt, primes);
        var cm1m2m2ntt = RLWE.MulRgsw(cm1m2ntt_, cm2ntt, csm2ntt);
        var d = RLWE.DecryptBGV(RLWE.NTTCipher2RLWE(cm1m2m2ntt, tables), sk);
        Console.WriteLine($"m1 * m2^2 = {(m1 * m2 * m2).ResModSigned(pm, t)}");
        Console.WriteLine($"          = {d}");
        Console.WriteLine($"          = {RLWE.ErrorsBGV(RLWE.NTTCipher2RLWE(cm1m2m2ntt, tables), sk).NormInf()}");
        Console.WriteLine();
    }
}

void testRotateNTT()
{
    var (n, p) = (16, 97);
    var pm = FG.QPoly().Pow(n) + 1;

    var w = new ZnBigInt(p, NumberTheory.NthRootsANTV1(p - 1, n, p).First());
    var nttInfos = RLWE.PrepareNTT(n, new(p), w);

    for (int i = 0; i < 5; i++)
    {
        var u = Rng.Next(p / 2) * RngSign;
        var m1 = RLWE.GenUnif(n, p);
        var xu = RLWE.XpowA(u, pm, new(p));
        var m2 = (m1 * (xu - 1)).ResModSigned(pm, p);
        var m3 = RLWE.RotateStep(m1, new(p), n, u);
        var m4ntt = RLWE.RotateStepNTT(RLWE.Rq2NTT(m1, nttInfos), nttInfos, u);
        var m4 = RLWE.NTT2Rq(m4ntt, nttInfos);
        Console.WriteLine($"m1 = {m1}");
        Console.WriteLine($"u  = {u,4}    X^u={xu}");
        Console.WriteLine($"m2 = {m2}");
        Console.WriteLine($"m3 = {m3}");
        Console.WriteLine($"m4 = {m4}");
        Console.WriteLine();
        Console.WriteLine();
        if (!m2.Equals(m3) || !m2.Equals(m4))
            throw new();
    }
}

void testBRNTT()
{
    // RLWE.NoiseOff();
    var N = 64;
    var n = N / 2;
    var level = 3;
    var (pm, sk, t, primes, sp, pk, _) = RLWE.SetupBGV(N, level);
    var q = primes[0];
    var qL = pk.Q;
    var B = RLWE.RNSGadgetBase(primes);
    Console.WriteLine($"BGV level = {level}, Primes = [{primes.Glue(", ")}]");
    Console.WriteLine($"pm = {pm} T = {t} q = {q} sp = {sp} qL = {qL} sigma = {RLWE.Sigma(n, t):f4}");
    Console.WriteLine($"sk = {sk}");
    Console.WriteLine($"pk => {pk.Params}");
    Console.WriteLine();

    var tables = RLWE.PrepareNTT(n, t, primes);

    GlobalStopWatch.AddLap();
    var brk = RLWE.BRKgswBGV(sk, pk, B);
    GlobalStopWatch.Show("Brk");

    GlobalStopWatch.AddLap();
    var brkNTT = RLWE.BRKgswNTTBGV(sk, pk, B, tables);
    GlobalStopWatch.Show("BrkNTT");

    Console.WriteLine();

    for (int k = 0; k < 10; ++k)
    {
        var ai = new Rational(Rng.Next(1, n + 1)).Signed(n);
        var bi = RLWE.GenUnif(n, n).CoefsExtended(n - 1);

        GlobalStopWatch.AddLap();
        var acc1 = RLWE.BlindRotategswBGV((ai, bi), brk, B, primes);
        GlobalStopWatch.Show("BR");
        var d_acc1 = RLWE.DecryptBGV(acc1, sk);

        GlobalStopWatch.AddLap();
        var acc2 = RLWE.BlindRotateNTTgswBGV((ai, bi), brkNTT, B, primes);
        GlobalStopWatch.Show("BRNTT");
        var d_acc2 = RLWE.DecryptBGV(acc2, sk);
        if (!d_acc1.Equals(d_acc2))
            throw new();

        Console.WriteLine();
    }
}

{
    testRgswNTT();
    testRotateNTT();
    testBRNTT();
}
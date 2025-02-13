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
// RecomputeAllPrimesUpTo(5000000);

bool CheckBR(Rational[] s, Rq pm, (Rational ai, Rational[] bi) ab, Rq actual, Rational q)
{
    var n = pm.Degree;
    var x = pm.X;
    var (ai, bi) = ab;
    var c = n / 4;
    var f = (2 * c + 1).SeqLazy(-c).Select(j => j * RLWE.XpowA(j, pm, new(2 * n)))
        .Aggregate(x.Zero, (acc, v) => acc + v).ResModSigned(pm, 2 * n);

    Console.WriteLine($"f           :{f}");
    Console.WriteLine($"f           :[{f.CoefsExtended(n - 1).Glue(", ", "{0,4}")}]");
    Console.WriteLine($"ai          :{ai}");
    Console.WriteLine($"bi          :[{bi.Glue(", ", "{0,4}")}]");
    Console.WriteLine($"blind rotate:{actual}");
    Console.WriteLine($"            :[{actual.CoefsExtended(n - 1).Glue(", ", "{0,4}")}]");

    // Testing result
    var u = (ai - bi.Zip(s).Aggregate(q.Zero, (sum, e) => sum + (e.First * e.Second)));
    var expected = (RLWE.XpowA((int)u.Num, pm, q) * f).ResModSigned(pm, q);
    var check = expected.Equals(actual);
    Console.WriteLine($"u= a - <b,s>:{u}");
    Console.WriteLine($"f*X^{u,-4}    :{expected}");
    Console.WriteLine($"f*X^{u,-4}    :[{expected.CoefsExtended(n - 1).Glue(", ", "{0,4}")}]");
    Console.WriteLine($"{(check ? "success" : "fail")}");
    Console.WriteLine();

    return check;
}

void testRgswMul()
{
    var N = 32;
    var n = N / 2;
    var level = 1;
    var (pm, sk, t, primes, sp, pk, swks) = RLWE.SetupBGV(N, level);
    var q = primes[0];
    var qL = pk.Q;
    var u = (int)BigInteger.Log2(t.Num);
    var B = new Rational(BigInteger.Pow(2, u));
    Console.WriteLine($"BGV level = {level}, Gadget Base = {B}");
    Console.WriteLine($"pm = {pm} T = {t} q = {q} sp = {sp} qL = {qL}");
    Console.WriteLine($"sk = {sk}");
    Console.WriteLine($"pk => {pk.Params}");
    Console.WriteLine();

    for (int k = 0; k < 5; ++k)
    {
        var m1 = RLWE.GenUnif(n, t);
        var cm1 = RLWE.EncryptBGV(m1, pk);
        var m2 = RLWE.GenUnif(n, t);
        // var m2 = RLWE.GenXpow(n);
        var _cm2 = RLWE.EncryptBGV(m2, pk);
        var (csm2, cm2) = RLWE.EncryptRgswBGV(m2, pk, B);
        var m1m2 = (m1 * m2).ResModSigned(pm, t);
        var cm1m2gsw = RLWE.MulRgsw(cm1, cm2, csm2, B);
        var cm1m2rlk = RLWE.MulRelinBGV(cm1, _cm2, swks[qL].rlk).ModSwitch(qL);
        var dm1m2gsw = RLWE.DecryptBGV(cm1m2gsw, sk);
        var dm1m2rlk = RLWE.DecryptBGV(cm1m2rlk, sk);
        Console.WriteLine($"m1      = {m1}");
        Console.WriteLine($"m2      = {m2}");
        Console.WriteLine($"m1 * m2 = {m1m2}");
        Console.WriteLine($"        = {dm1m2gsw}");
        Console.WriteLine($"        = {dm1m2rlk}");
        Console.WriteLine($"egsw    = {RLWE.ErrorsBGV(cm1m2gsw, sk).NormInf()}");
        Console.WriteLine($"erlk    = {RLWE.ErrorsBGV(cm1m2rlk, sk).NormInf()}");
        Console.WriteLine();
        if (!dm1m2gsw.Equals(m1m2))
            throw new();
    }
}

void testBR()
{
    var N = 32;
    var n = N / 2;
    var level = 1;
    var (pm, sk, t, primes, sp, pk, swks) = RLWE.SetupBGV(N, level);
    var q = primes[0];
    var qL = pk.Q;
    var u = (int)BigInteger.Log2(t.Num);
    var B = new Rational(BigInteger.Pow(2, u));
    Console.WriteLine($"BGV level = {level}, Gadget Base = {B}");
    Console.WriteLine($"pm = {pm} T = {t} q = {q} sp = {sp} qL = {qL}");
    Console.WriteLine($"sk = {sk}");
    Console.WriteLine($"pk => {pk.Params}");
    Console.WriteLine();

    var brk = RLWE.BRKgswBGV(sk, pk, B);
    var rlwe0 = RLWE.EncryptRgswBGV(pm.One, pk, B, noiseMode: false);

    var s = n.SeqLazy().Select(i => sk[i].Signed(t)).ToVec();

    var nbFails = 0;
    for (int k = 0; k < 5; ++k)
    {
        var ai = new Rational(Rng.Next(1, n + 1)).Signed(n);
        var bi = RLWE.GenUnif(n, n).CoefsExtended(n - 1);

        var acc = RLWE.BlindRotategswBGV((ai, bi), B, rlwe0, brk);
        var actual = RLWE.DecryptBGV(acc, sk);

        if (!CheckBR(s.ToArray(), pm, (ai, bi), actual, n * t.One))
            Console.WriteLine($"fail:{++nbFails}");
    }

    Console.WriteLine($"nbFails:{nbFails}");
    Console.WriteLine();
}

void testBootstrapping()
{
    RLWE.NoiseOff();
    GlobalStopWatch.Restart();
    var N = 32;
    var n = N / 2;
    var level = 3;

    Console.WriteLine("####       Start        ####");
    var (pm, sk, t, primes, sp, pk, swks) = RLWE.SetupBGV(N, level);
    var q = primes[0];
    var qL = pk.Q;

    var u = (int)BigInteger.Log2(t.Num);
    var B = new Rational(BigInteger.Pow(2, u));
    var rlk = swks[qL].rlk;
    var brk = RLWE.BRKgswBGV(sk, pk, B);
    var ak = N.SeqLazy()
        .Select(j => RLWE.SWKBGV(pm, sk, sk.Substitute(pm.X.Pow(j)).ResModSigned(pm, t), t, qL, sp.Pow(level)))
        .ToArray();

    Console.WriteLine($"BGV level = {level}, Gadget Base = {B}");
    Console.WriteLine($"pm = {pm} T = {t} q = {q} sp = {sp} qL = {qL}");
    Console.WriteLine($"sk = {sk}");
    Console.WriteLine($"pk => {pk.Params}");

    Console.WriteLine();

    Console.WriteLine("#### Bootstrapping test ####");
    for (int k = 0; k < 15; ++k)
    {
        GlobalStopWatch.AddLap();
        var m = RLWE.GenUnif(n, t);
        var m2 = (m * m).ResModSigned(pm, t);
        var cm = RLWE.EncryptBGV(m, pk); // level qL
        var ct = RLWE.MulRelinBGV(cm, cm, rlk).ModSwitch(q); // level q0

        var (ctboot, ctsm) = RLWE.Bootstrapping(ct, B, pk, ak, brk);
        Console.WriteLine($"ct     {ct.Params}");
        Console.WriteLine($"ctboot {ctboot.Params}");

        var m2boot = RLWE.DecryptBGV(ctboot, sk);
        Console.WriteLine($"m       = {m}");
        Console.WriteLine($"m^2     = {m2}");
        Console.WriteLine($"ctboot  = {m2boot}");
        Console.WriteLine($"eboot   = {RLWE.ErrorsBGV(ctboot, sk).NormInf()}");
        Console.WriteLine($"emodsw  = {RLWE.ErrorsBGV(ct.ModSwitch(qL), sk).NormInf()}");
        Console.WriteLine();
        if (!m2.Equals(m2boot))
            throw new();

        var c2 = RLWE.MulRelinBGV(ctboot, ctboot, rlk).ModSwitch(swks[qL].nextMod);
        var d2 = RLWE.DecryptBGV(c2, sk);
        var m4 = (m2 * m2).ResModSigned(pm, t);
        Console.WriteLine($"m^4 = {m4}");
        Console.WriteLine($"    = {d2}");

        GlobalStopWatch.Show();
        Console.WriteLine();
        if (!d2.Equals(m4))
            throw new();
    }

    Console.WriteLine("####        End         ####");
    GlobalStopWatch.Show();
}

{
    // testRgswMul();
    // testBR();
    testBootstrapping();

    // BGV level = 3, Gadget Base = 64
    // pm = x^16 + 1 T = 97 q = 43457 sp = 90599 qL = 12815095448073162241
    // sk = x^15 - x^14 - x^12 + x^10 - x^9 - x^8 + x^7 - x^5 + x^4 - x^3 - x^2 - 1
    // pk => RLWECipher Q:12815095448073162241    T:97    PM:x^16 + 1
    // 
    // ct     RLWECipher Q:43457    T:97    PM:x^16 + 1
    // ctboot RLWECipher Q:12815095448073162241    T:97    PM:x^16 + 1
    // m       = 10*x^15 + 38*x^14 - 36*x^13 - 21*x^12 - x^11 - 17*x^10 - 21*x^9 - 44*x^8 - 11*x^7 - 36*x^6 - 36*x^5 - 27*x^4 + 34*x^3 + 9*x^2 - 8*x - 16
    // m^2     = -22*x^15 - 26*x^14 - 8*x^13 - 22*x^12 - 22*x^11 - 10*x^10 - 34*x^9 + 28*x^8 + 42*x^7 - 48*x^6 - 41*x^5 + 13*x^4 + 31*x^3 + 32*x^2 - 13*x - 29
    // ctboot  = -22*x^15 - 26*x^14 - 8*x^13 - 22*x^12 - 22*x^11 - 10*x^10 - 34*x^9 + 28*x^8 + 42*x^7 - 48*x^6 - 41*x^5 + 13*x^4 + 31*x^3 + 32*x^2 - 13*x - 29
    // eboot   = 223715756
    // emodsw  = 78146220257711971
    // 
    // m^4 = 23*x^15 - 33*x^14 - 36*x^13 - 12*x^12 - 13*x^11 + 43*x^10 - 9*x^9 + 45*x^8 - 14*x^7 + 10*x^6 - 11*x^5 + 31*x^4 - 36*x^3 + x^2 - 3*x - 19
    //     = 23*x^15 - 33*x^14 - 36*x^13 - 12*x^12 - 13*x^11 + 43*x^10 - 9*x^9 + 45*x^8 - 14*x^7 + 10*x^6 - 11*x^5 + 31*x^4 - 36*x^3 + x^2 - 3*x - 19
    // #  Time:5.799s
    // 
}
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
RngSeed(259128);
Ring.DisplayPolynomial = MonomDisplay.StarCaret;

(Rational pi, Rational invi)[] RNSGadgetBase(Rational[] primes)
{
    var qL = primes.Aggregate((pi, pj) => pi * pj);
    return primes
        .Select(e => (e, (new Rational(InvModPbezbigint((qL / e).Num, e.Num)).Signed(e) * (qL / e)).Signed(qL)))
        .ToArray();
}

Vec<RLWECipher> DecompRNS(RLWECipher cipher, Rational[] primes)
{
    var (pm, t, _) = cipher.PM_T_Q;
    var qL = primes.Aggregate((pi, pj) => pi * pj);
    return primes.Select(pi => new RLWECipher(cipher.A.CoefsModSigned(pi), cipher.B.CoefsModSigned(pi), pm, t, qL))
        .ToVec();
}

Vec<RLWECipher> PKBGV(Rq pm, Rq sk, Rational t, Rational[] primes)
{
    return primes.Select(pi => RLWE.PKBGV(pm, sk, t, pi)).ToVec();
}

(Vec<RLWECipher> csm, Vec<RLWECipher> cm) EncryptRgswBGV(Rq m, RLWECipher pk, (Rational pi, Rational invi)[] B,
    bool noiseMode = true)
{
    var (pm, t, q) = pk.PM_T_Q;
    var z = pm.Zero;
    RLWECipher ctZero() => RLWE.EncryptBGV(z, pk, noiseMode);
    var cm = B.Select(e => ctZero() + (m * e.invi, z, pm, t, q)).ToVec();
    var csm = B.Select(e => ctZero() + (z, -m * e.invi, pm, t, q)).ToVec();
    return (csm, cm);
}

Vec<RLWECipher> MulRgsw(Vec<RLWECipher> c1, Vec<RLWECipher> cm, Vec<RLWECipher> csm)
{
    var c1ac2m = c1.Zip(cm).Select(e => e.First.A * e.Second).ToVec();
    var c1bc2sm = c1.Zip(csm).Select(e => e.First.B * e.Second).ToVec();
    return c1ac2m - c1bc2sm;
}

RLWECipher RecompRNS(Vec<RLWECipher> ciphers, (Rational pi, Rational invi)[] B)
{
    var (pm, t, _) = ciphers[0].PM_T_Q;
    var qL = B.Select(e => e.pi).Aggregate(Rational.KOne(), (acc, pi) => acc * pi);
    var a = ciphers.Zip(B).Select(e => e.First.A * e.Second.invi).Aggregate((ci, cj) => ci + cj).CoefsModSigned(qL);
    var b = ciphers.Zip(B).Select(e => e.First.B * e.Second.invi).Aggregate((ci, cj) => ci + cj).CoefsModSigned(qL);
    return new RLWECipher(a, b, pm, t, qL);
}

void testRNSGadgetBase()
{
    var (N, level) = (16, 5);
    var rlwe = new RLWE(N, level, bootstrappingMode: false);
    rlwe.Show();
    var (n, pm, sk, t, q, pk, rlk) = rlwe;
    var qL = pk.Q;

    var bs = RNSGadgetBase(rlwe.Primes);
    bs.Println("bs");
    for (int i = 0; i < 10; i++)
    {
        var m = RLWE.GenUnif(n, qL);
        Console.WriteLine($"m = {m}");
        var dec = rlwe.Primes.Select(qi => m.CoefsModSigned(qi)).ToVec();
        dec.Println();
        var m2 = bs.Zip(dec).Select(e => e.First.invi * e.Second).Aggregate((ci, cj) => ci + cj).CoefsModSigned(qL);
        Console.WriteLine($" = {m2}");
        Console.WriteLine(m - m2);
        Console.WriteLine();
    }
}

void testEncryptDecrypt()
{
    // Weak parameters

    var (N, level) = (16, 5);
    var n = N / 2;
    var (pm, sk, t, primes, sp, pk, rlks) = RLWE.SetupBGV(N, level);
    var q = primes[0];
    var qL = pk.Q;

    Console.WriteLine($"BGV level = {level}");
    Console.WriteLine($"pm = {pm} T = {t} q = {q} sp = {sp} qL = {qL} Sigma = {RLWE.Sigma(n, t):f4}");
    Console.WriteLine($"sk = {sk}");
    Console.WriteLine($"pk => {pk.Params}");
    Console.WriteLine($"primes = [{primes.Glue(", ")}]");
    Console.WriteLine();

    var rnsBs = RNSGadgetBase(primes);
    for (int k = 0; k < 5; ++k)
    {
        var m1 = RLWE.GenUnif(n, t);
        var cm1 = RLWE.EncryptBGV(m1, pk);
        var _cm1 = DecompRNS(cm1, primes);
        var cm2 = RecompRNS(_cm1, rnsBs);
        var m2 = RLWE.DecryptBGV(cm2, sk);
        Console.WriteLine($"m1 = {m1}");
        Console.WriteLine($"m2 = {m2}");
        Console.WriteLine();
        if (!m1.Equals(m2))
            throw new();
    }
}

// void testMulRgsw()
{
    // Weak parameters
    
    var (N, level) = (16, 5);
    var n = N / 2;
    var (pm, sk, t, primes, sp, pk, rlks) = RLWE.SetupBGV(N, level);
    var q = primes[0];
    var qL = pk.Q;

    Console.WriteLine($"BGV level = {level}");
    Console.WriteLine($"pm = {pm} T = {t} q = {q} sp = {sp} qL = {qL} Sigma = {RLWE.Sigma(n, t):f4}");
    Console.WriteLine($"sk = {sk}");
    Console.WriteLine($"pk => {pk.Params}");
    Console.WriteLine($"primesRNS = [{primes.Glue(", ")}]");
    Console.WriteLine();

    var rnsBs = RNSGadgetBase(primes);
    rnsBs.Println("rnsBs");
    for (int k = 0; k < 5; ++k)
    {
        var m1 = RLWE.GenUnif(n, t);
        var cm1 = RLWE.EncryptBGV(m1, pk);
        var _cm1 = DecompRNS(cm1, primes);
        var m2 = RLWE.GenUnif(n, t);
        var _cm2 = RLWE.EncryptBGV(m2, pk);
        var (csm2, cm2) = EncryptRgswBGV(m2, pk, rnsBs);
        var m1m2 = (m1 * m2).ResModSigned(pm, t);
        var cm1m2gsw = MulRgsw(_cm1, cm2, csm2).Sum();
        var cm1m2rlk = RLWE.MulRelinBGV(cm1, _cm2, rlks[qL].rlk).ModSwitch(qL);
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
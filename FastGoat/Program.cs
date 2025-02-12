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

Vec<Rq> DecompRq(Rq a, Rational bs, Rational mod)
{
    var size = (int)(BigInteger.Log10(mod.Num) / BigInteger.Log10(bs.Num)) + 1;
    var a0 = a.CoefsModSigned(mod);
    var queue = new Queue<Rq>();
    for (int i = 0; i < size; ++i)
    {
        var a1 = (a0 / bs).TruncPoly();
        var r = a0 - bs * a1;
        a0 = a1;
        queue.Enqueue(r);
    }

    return queue.ToVec();
}

(Vec<RLWECipher> csm, Vec<RLWECipher> cm) EncryptRgswBGV(Rq m, RLWECipher pk, Rational B, bool noiseMode = true)
{
    var (pm, t, q) = pk.PM_T_Q;
    var size = (int)(BigInteger.Log10(q.Num) / BigInteger.Log10(B.Num)) + 1;
    var z = pm.Zero;
    RLWECipher ctZero() => RLWE.EncryptBGV(z, pk, noiseMode);
    var cm = size.SeqLazy().Select(i => ctZero() + (m * B.Pow(i), z, pm, t, q)).ToVec();
    var csm = size.SeqLazy().Select(i => ctZero() + (z, -m * B.Pow(i), pm, t, q)).ToVec();
    return (csm, cm);
}

RLWECipher MulRgsw(RLWECipher c1, Vec<RLWECipher> cm, Vec<RLWECipher> csm, Rational B)
{
    var q = c1.Q;
    var c1ac2m = DecompRq(c1.A, B, q).Zip(cm).Select(e => e.First * e.Second).ToVec().Sum();
    var c1bc2sm = DecompRq(c1.B, B, q).Zip(csm).Select(e => e.First * e.Second).ToVec().Sum();
    return c1ac2m - c1bc2sm;
}

Rq XpowA(int a, Rq pm, Rational q)
{
    var x = pm.X;
    if (a == 0)
        return x.One;

    var n = pm.Degree;
    var sgn = a > 0 || a % n == 0 ? 0 : 1;
    return (x.Pow(IntExt.AmodP(a, n)) * (-1).Pow((a / n + sgn) % 2)).CoefsModSigned(q);
}

(RLWECipher ct1, RLWECipher ctprep) CtPrep(RLWECipher ct)
{
    var (pm, t, q) = ct.PM_T_Q;
    var N = new Rational(2 * pm.Degree);
    var a1 = (N * ct.A).CoefsModSigned(q);
    var b1 = (N * ct.B).CoefsModSigned(q);
    var a2 = (N * ct.A - a1) / q;
    var b2 = (N * ct.B - b1) / q;
    return ((a1, b1, pm, t, q), (a2, b2, pm, t, N));
}

((Vec<RLWECipher> csm, Vec<RLWECipher> cm) plus, (Vec<RLWECipher> csm, Vec<RLWECipher> cm) minus)[] 
    BRKgswBGV(Rq sk, RLWECipher pk, Rational B)
{
    var pm = pk.PM;
    var enc = (int k) => EncryptRgswBGV(pm.One * k, pk, B);
    var n = pm.Degree;
    return n.SeqLazy().Select(i => sk[i])
        .Select(c => (plus: c.IsOne() ? enc(1) : enc(0), minus: (-c).IsOne() ? enc(1) : enc(0)))
        .ToArray();
}

RLWECipher BlindRotategswBGV((Rational ai, Rational[] bi) ab, Rq f, Rational B, 
    (Vec<RLWECipher> csm, Vec<RLWECipher> cm) rlwe0,
    ((Vec<RLWECipher> csm, Vec<RLWECipher> cm) plus, (Vec<RLWECipher> csm, Vec<RLWECipher> cm) minus)[] brk)
{
    var (pm, t, qL) = brk[0].minus.cm[0].PM_T_Q;
    var n = pm.Degree;
    var N = new Rational(n);
    var x = pm.X;

    var beta = ab.bi;
    var alpha = (int)ab.ai.Num;
    var xalpha = XpowA(alpha, pm, N);
    
    var (encSOne0, encOne0) = rlwe0;
    RLWECipher acc = ((f * xalpha).ResModSigned(pm, N), x.Zero, pm, t, qL);
    
    for (int i = 0; i < n; i++)
    {
        var (encSi_plus, encSi_minus) = brk[i];
        var ai = (int)beta[i].Opp().Num;

        var exai = (XpowA(ai, pm, N) - 1);
        var cxai = encSi_plus.cm.Select(e => exai * e).ToVec();
        var csxai = encSi_plus.csm.Select(e => exai * e).ToVec();
        
        var ex_ai = (XpowA(-ai, pm, N) - 1);
        var cx_ai = encSi_minus.cm.Select(e => ex_ai * e).ToVec();
        var csx_ai = encSi_minus.csm.Select(e => ex_ai * e).ToVec();

        var acci = encOne0 + cxai + cx_ai;
        var sacci = encSOne0 + csxai + csx_ai;
        acc = MulRgsw(acc, acci, sacci, B);
    }

    return acc;
}

bool CheckBR(Rational[] s, Rq pm, (Rational ai, Rational[] bi) ab, Rq f, Rq actual, Rational q)
{
    var n = pm.Degree;
    var (ai, bi) = ab;
    Console.WriteLine($"f           :{f}");
    Console.WriteLine($"f           :[{f.CoefsExtended(n - 1).Glue(", ", "{0,4}")}]");
    Console.WriteLine($"ai          :{ai}");
    Console.WriteLine($"bi          :[{bi.Glue(", ", "{0,4}")}]");
    Console.WriteLine($"blind rotate:{actual}");
    Console.WriteLine($"            :[{actual.CoefsExtended(n - 1).Glue(", ", "{0,4}")}]");

    // Testing result
    var u = (ai - bi.Zip(s).Aggregate(q.Zero, (sum, e) => sum + (e.First * e.Second)));
    var expected = (XpowA((int)u.Num, pm, q) * f).ResModSigned(pm, q);
    var check = expected.Equals(actual);
    Console.WriteLine($"u= a - <b,s>:{u}");
    Console.WriteLine($"f*X^{u,-4}    :{expected}");
    Console.WriteLine($"f*X^{u,-4}    :[{expected.CoefsExtended(n - 1).Glue(", ", "{0,4}")}]");
    Console.WriteLine($"{(check ? "success" : "fail")}");
    Console.WriteLine();

    return check;
}

RLWECipher RepackingBGV(int n, RLWECipher[] accs, RLWECipher[] autSk)
{
    var acc0 = accs[0];
    var (pm, t, q) = acc0.PM_T_Q;
    var d = pm.Degree;
    var x = pm.X;
    var CT = Ring.Matrix(acc0.One, d, d + 1);
    for (int i = 0; i < n; i++)
    {
        CT[i, n] = accs[i];
        for (int j = 1; j < d / n; j++)
            CT[i, n] = CT[i, n] + XpowA(n * j, pm, new(d)) * CT[i + j * n, n];
    }

    for (int k = n; k > 1; k /= 2)
    {
        for (int i = 0; i < k / 2; i++)
        {
            CT[i, k / 2] = CT[i, k] + x.Pow(k / 2) * CT[i + k / 2, k];
            var K = 1 + 2 * d / k;
            var crot = RLWE.AutoMorphBGV(CT[i, k] - x.Pow(k / 2) * CT[i + k / 2, k], K, autSk[K]).ModSwitch(q);
            CT[i, k / 2] = CT[i, k / 2] + crot;
        }
    }

    return CT[0, 1];
}

(RLWECipher ctboot, RLWECipher ctsm) Bootstrapping(RLWECipher ct, Rational B, RLWECipher pk, RLWECipher[] skAut,
    ((Vec<RLWECipher> csm, Vec<RLWECipher> cm) plus, (Vec<RLWECipher> csm, Vec<RLWECipher> cm) minus)[] brk)
{
    var (pm, t, qL) = pk.PM_T_Q;
    var n = pm.Degree;

    // 1. Extract
    var (ct1, ctprep) = CtPrep(ct);
    var extract = RLWE.Extract(ctprep);

    // 2. BlindRotate
    var c = n / 4;
    var f = (2 * c + 1).SeqLazy(-c).Select(j => j * XpowA(j, pm, t)).Aggregate((v0, v1) => v0 + v1).ResModSigned(pm, t);

    var rlwe0 = EncryptRgswBGV(pm.One, pk, B, noiseMode: false);
    var ni = new Rational(InvModPbezbigint(n, qL.Num)).Signed(qL);
    var seqBR = new List<RLWECipher>();
    foreach (var ab in extract)
        seqBR.Add(ni * BlindRotategswBGV(ab, f, B, rlwe0, brk));

    // Step 3. Repacking
    var ctsm = RepackingBGV(n, seqBR.ToArray(), skAut);
    var ctboot = ctsm - ct1;
    var h = (t / (2 * n)).Trunc;
    return (ctboot * h, ctsm);
}

void testRgswMul()
{
    var N = 32;
    var n = N / 2;
    var level = 1;
    var (pm, sk, t, primes, pk, swks) = RLWE.SetupBGV(N, level);
    var q = primes[0];
    var sp = Rational.NthRoot(swks[pk.Q].skPow[0].Q / pk.Q, level);
    var qL = pk.Q;
    var u = (int)BigInteger.Log2(t.Num);
    var B = new Rational(BigInteger.Pow(2, u));
    Console.WriteLine($"BGV level = {level}, Gadget Base = {B}");
    Console.WriteLine($"pm = {pm} T = {t} q = {q} sp = {sp} qL = {qL}");
    Console.WriteLine($"sk = {sk}");
    Console.WriteLine($"pk => {pk.Params}");
    Console.WriteLine();
    
    for(int k = 0; k < 5; ++k)
    {
        var m1 = RLWE.GenUnif(n, t);
        var cm1 = RLWE.EncryptBGV(m1, pk);
        var m2 = RLWE.GenUnif(n, t);
        // var m2 = RLWE.GenXpow(n);
        var _cm2 = RLWE.EncryptBGV(m2, pk);
        var (csm2, cm2) = EncryptRgswBGV(m2, pk, B);
        var m1m2 = (m1 * m2).ResModSigned(pm, t);
        var cm1m2gsw = MulRgsw(cm1, cm2, csm2, B);
        var cm1m2rlk = RLWE.MulRelinBGV(cm1, _cm2, swks[qL].skPow[2]).ModSwitch(qL);
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
    var (pm, sk, t, primes, pk, swks) = RLWE.SetupBGV(N, level);
    var q = primes[0];
    var sp = Rational.NthRoot(swks[pk.Q].skPow[0].Q / pk.Q, level);
    var qL = pk.Q;
    var u = (int)BigInteger.Log2(t.Num);
    var B = new Rational(BigInteger.Pow(2, u));
    Console.WriteLine($"BGV level = {level}, Gadget Base = {B}");
    Console.WriteLine($"pm = {pm} T = {t} q = {q} sp = {sp} qL = {qL}");
    Console.WriteLine($"sk = {sk}");
    Console.WriteLine($"pk => {pk.Params}");
    Console.WriteLine();
    
    var brk = BRKgswBGV(sk, pk, B);
    var rlwe0 = EncryptRgswBGV(pm.One, pk, B, noiseMode: false);

    var x = pm.X;
    var c = n / 2 - 1;
    var f = (2 * c + 1).SeqLazy(-c).Select(j => j * XpowA(j, pm, t))
        .Aggregate(x.Zero, (acc, v) => acc + v).ResModSigned(pm, t);

    var s = n.SeqLazy().Select(i => sk[i].Signed(t)).ToVec();

    var nbFails = 0;
    for (int k = 0; k < 5; ++k)
    {
        var ai = new Rational(Rng.Next(1, n + 1)).Signed(n);
        var bi = RLWE.GenUnif(n, n).CoefsExtended(n - 1);

        var acc = BlindRotategswBGV((ai, bi), f, B, rlwe0, brk);
        var actual = RLWE.DecryptBGV(acc, sk);
        
        if (!CheckBR(s.ToArray(), pm, (ai, bi), f, actual, n * t.One)) 
            Console.WriteLine($"fail:{++nbFails}");

    }

    Console.WriteLine($"nbFails:{nbFails}");
    Console.WriteLine();
}

void testBootstrapping()
{
    GlobalStopWatch.Restart();
    var N = 32;
    var n = N / 2;
    var level = 3;
    var (pm, sk, t, primes, pk, swks) = RLWE.SetupBGV(N, level);
    var q = primes[0];
    var sp = Rational.NthRoot(swks[pk.Q].skPow[0].Q / pk.Q, level);
    var qL = pk.Q;
    var u = (int)BigInteger.Log2(t.Num);
    var B = new Rational(BigInteger.Pow(2, u));
    Console.WriteLine($"BGV level = {level} Gadget Base = {B}");
    Console.WriteLine($"pm = {pm} T = {t} q = {q} sp = {sp} qL = {qL}");
    Console.WriteLine($"sk = {sk}");
    Console.WriteLine($"pk => {pk.Params}");

    Console.WriteLine();

    var rlk = swks[qL].skPow[2];
    var brk = BRKgswBGV(sk, pk, B);
    var ak = N.SeqLazy()
        .Select(j => RLWE.SWKBGV(pm, sk, sk.Substitute(pm.X.Pow(j)).ResModSigned(pm, t), t, qL, sp.Pow(level)))
        .ToArray();

    for (int k = 0; k < 5; ++k)
    {
        GlobalStopWatch.AddLap();
        var m = RLWE.GenUnif(n, t);
        var m2 = (m * m).ResModSigned(pm, t);
        var cm = RLWE.EncryptBGV(m, pk); // level qL
        var ct = RLWE.MulRelinBGV(cm, cm, rlk).ModSwitch(q); // level q0

        var (ctboot, ctsm) = Bootstrapping(ct, B, pk, ak, brk);
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
}

{
    testRgswMul();
    testBR();
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
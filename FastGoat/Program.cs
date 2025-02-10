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

Rq InvRq(Rq a, Rq f, Rational t)
{
    var a0 = a.ResModSigned(f, t);
    var ia = Ring.FastBezout(a0, f).Item1.ResModSigned(f, t);
    var q = (a0 * ia).ResModSigned(f, t);
    if (q.Degree != 0 || q.IsZero())
        throw new($"a={a} not invertible");

    var iq = InvModPbezbigint(q[0].Num, t.Num);
    var iaf = (new Rational(iq) * ia).CoefsModSigned(t);
    return iaf;
}

(Rq[] invs, Rq[] facts) CrtBase(int N, Rational t)
{
    if (!int.IsPow2(N))
        throw new();
    if (!Primes10000.Contains((int)t.Num) || !t.Mod(N).IsOne())
        throw new();
    var n = N / 2;
    var t0 = (int)t.Num;
    var sols = SolveAll_k_pow_m_equal_one_mod_n_strict_long(t0, N)
        .Select(e => new Rational(e).Signed(t0))
        .Order()
        .ToArray();

    var pm = FG.QPoly().Pow(n) + 1;
    var facts = sols.Select(e => pm.X - e).ToArray();
    var invs = facts.Select(mi => (quo: pm / mi, mi)).Select(e => e.quo * InvRq(e.quo, e.mi, t)).ToArray();
    facts.Println($"pm = {pm}");
    Console.WriteLine();

    return (invs, facts);
}

Rq SolveCRTRq(Rq[] ais, Rq[] invs, Rq pm, Rational t)
{
    return ais.Zip(invs).Aggregate(pm.Zero, (acc, ej) => (acc + ej.First * ej.Second).ResModSigned(pm, t));
}

Rq[] InvCRTRq(Rq x, Rq[] mods, Rational t) => mods.Select(mod => x.ResModSigned(mod, t)).ToArray();

void testCRT()
{
    RecomputeAllPrimesUpTo(1000000);
    var x = FG.QPoly();
    for (int i = 2; i < 9; ++i)
    {
        var N = 1 << i;
        var n = N / 2;
        var t0 = Primes10000.First(t0 => t0 % N == 1);
        var t = new Rational(t0);
        var pm = x.Pow(n) + 1;

        Console.WriteLine($"N = {N} pm = {pm} t = {t0}");
        var crtBase = CrtBase(N, t);

        for (int j = 0; j < 5; j++)
        {
            var m1 = RLWE.GenUnif(n, t);
            var crt1 = InvCRTRq(m1, crtBase.facts, t).ToArray();
            crt1.Zip(crtBase.facts).Println(e => $"m = {e.First,-6} mod ({e.Second}) mod {t}", $"m1 = {m1}");
            var rm1 = SolveCRTRq(crt1, crtBase.invs, pm, t);
            Console.WriteLine($"  => {rm1}");
            Console.WriteLine();

            var m2 = RLWE.GenUnif(n, t);
            var crt2 = InvCRTRq(m2, crtBase.facts, t).ToArray();
            crt2.Zip(crtBase.facts).Println(e => $"m = {e.First,-6} mod ({e.Second}) mod {t}", $"m2 = {m2}");
            var rm2 = SolveCRTRq(crt2, crtBase.invs, pm, t);
            Console.WriteLine($"  => {rm2}");
            Console.WriteLine();

            var m1m2 = (m1 * m2).ResModSigned(pm, t);
            var crt12a = InvCRTRq(m1m2, crtBase.facts, t).ToArray();
            crt12a.Zip(crtBase.facts).Println(e => $"m = {e.First,-6} mod ({e.Second}) mod {t}", $"m1m2a = {m1m2}");
            var rm1m2a = SolveCRTRq(crt12a, crtBase.invs, pm, t);
            Console.WriteLine($"  => {rm1m2a}");
            Console.WriteLine();

            var crt12b = crt1.Zip(crt2).Zip(crtBase.facts)
                .Select(e => (e.First.First * e.First.Second).ResModSigned(e.Second, t)).ToArray();
            crt12b.Zip(crtBase.facts).Println(e => $"m = {e.First,-6} mod ({e.Second}) mod {t}", $"m1m2b = {m1m2}");
            var rm1m2b = SolveCRTRq(crt12b, crtBase.invs, pm, t);
            Console.WriteLine($"  => {rm1m2b}");
            Console.WriteLine();

            if (!m1.Equals(rm1) || !m2.Equals(rm2) || !m1m2.Equals(rm1m2a) || !m1m2.Equals(rm1m2b))
                throw new("Fail");
        }
    }
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

void testKeySwitch()
{
    var (N, level) = (32, 1);
    var n = N / 2;
    var (pm, sk, t, primes, pk, swks) = RLWE.SetupBGV(N, level);
    var q = primes[0];
    var qL = pk.Q;
    var spqL = swks[pk.Q].skPow[0].Q;
    var sp = Rational.NthRoot(spqL / qL, level);
    Console.WriteLine($"pm = {pm} T = {t} Primes = [{primes.Glue(", ")}] q = {q} qL = {qL} sp = {sp}");
    Console.WriteLine($"sk = {sk}");
    Console.WriteLine($"pk => {pk.Params}");

    var swk = swks[qL].skPow[1];
    var nextMod = swks[qL].nextMod;

    for (int i = 0; i < 50; ++i)
    {
        var m = RLWE.GenUnif(n, t);
        var cm = RLWE.EncryptBGV(m, pk);
        var cm1 = RLWE.SwitchKeyBGV(cm, swk).ModSwitch(nextMod);
        var d_m1 = RLWE.DecryptBGV(cm1, sk);
        Console.WriteLine($"m     = {m}");
        Console.WriteLine($"      = {d_m1}");
        Console.WriteLine($"em    = {RLWE.ErrorsBGV(cm, sk)} level {cm.Q}");
        Console.WriteLine($"em1   = {RLWE.ErrorsBGV(cm1, sk)} level {cm1.Q}");
        Console.WriteLine();
        if (!m.Equals(d_m1))
            throw new();
    }
}

void testAutoMorph()
{
    var (N, level) = (32, 1);
    var n = N / 2;
    var (pm, sk, t, primes, pk, swks) = RLWE.SetupBGV(N, level);
    var q = primes[0];
    var qL = pk.Q;
    var spqL = swks[pk.Q].skPow[0].Q;
    var sp = Rational.NthRoot(spqL / qL, level);
    Console.WriteLine($"pm = {pm} T = {t} Primes = [{primes.Glue(", ")}] q = {q} qL = {qL} sp = {sp}");
    Console.WriteLine($"sk = {sk}");
    Console.WriteLine($"pk => {pk.Params}");

    var x = pm.X;
    var seqMods = (level + 1).SeqLazy(1).Select(i => primes.Take(i).Aggregate((pi, pj) => pi * pj)).ToArray();
    var ak = RLWE.LeveledAutMorphSwitchKeyGenBGV(level, pm, sk, t, seqMods, sp);
    foreach (var k in IntExt.Coprimes(2 * n))
    {
        var xk = x.Pow(k);
        for (var i = 0; i < 5; ++i)
        {
            var m = RLWE.GenUnif(n, t);
            var mk = m.Substitute(xk).ResModSigned(pm, t);
            Console.WriteLine($"m       :{m}");
            Console.WriteLine($"{$"m(x^{k})",-8}:{mk}");

            var cm = RLWE.EncryptBGV(m, pk);
            var ck = RLWE.AutoMorphBGV(cm, k, ak[qL].autSk[k]).ModSwitch(ak[qL].nextMod);
            var dk = RLWE.DecryptBGV(ck, sk);
            Console.WriteLine($"        :{dk}");
            Console.WriteLine($"err     :{RLWE.ErrorsBGV(ck, sk)} level {ck.Q}");
            if (!dk.Equals(mk))
                throw new($"k:{k} step[{i}] {dk.Div(mk)}");

            Console.WriteLine();
        }
    }
}

{
    // testCRT();
    // testKeySwitch();
    // testAutoMorph();
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

void testXpowAlpha()
{
    var (n, q) = (8, 97);
    var x = FG.QPoly();
    var pm = x.Pow(n) + 1;
    for (int a = 0; a < q; a++)
    {
        var xi0 = x.Pow(a).ResModSigned(pm, q);
        var xi1 = XpowA(a, pm, new(q));
        Console.WriteLine($"i = {a} x^i = {xi1}");
        if (!xi0.Equals(xi1))
            throw new($"expected:{xi0} actual:{xi1}");
    }
}

(RLWECipher plus, RLWECipher minus)[] BRKBGV(Rq sk, RLWECipher pk)
{
    var pm = pk.PM;
    var enc = (int k) => RLWE.EncryptBGV(pm.One * k, pk);
    var n = pm.Degree;

    Console.WriteLine();
    return n.SeqLazy().Select(i => sk[i])
        .Select(c => (plus: c.IsOne() ? enc(1) : enc(0), minus: (-c).IsOne() ? enc(1) : enc(0)))
        .ToArray();
}

RLWECipher BlindRotateBGV((Rational ai, Rational[] bi) ab, Rq f, (RLWECipher plus, RLWECipher minus)[] brk,
    RLWECipher rlk)
{
    var (pm, t, q) = brk[0].minus.PM_T_Q;
    var n = pm.Degree;
    var N = new Rational(n);
    var x = pm.X;

    var beta = ab.bi;
    var alpha = (int)ab.ai.Num;
    var xalpha = XpowA(alpha, pm, N);
    RLWECipher encOne = (x.One, x.Zero, pm, t, q);

    var acc = ((f * xalpha).ResModSigned(pm, N) * encOne);
    for (int i = 0; i < n; i++)
    {
        var (encSi_plus, encSi_minus) = brk[i];
        var ai = (int)beta[i].Opp().Num;

        var exai = (XpowA(ai, pm, N) - 1).CoefsModSigned(N);
        var cxai = exai * encSi_plus;
        var ex_ai = (XpowA(-ai, pm, N) - 1).CoefsModSigned(N);
        var cx_ai = ex_ai * encSi_minus;
        
        var acci = encOne + cxai + cx_ai;
        acc = RLWE.MulRelinBGV(acc, acci, rlk).ModSwitch(q);
    }
    
    var ni = new Rational(InvModPbezbigint(n, t.Num));
    return (ni * acc);
}

void RPShow(Rq sk, RLWECipher[,] mat, string header)
{
    Console.WriteLine(header);
    var (nbRows, nbCols) = (mat.GetLength(0), mat.GetLength(1));
    var gr = nbRows.Range().Grid2D(nbCols.Range()).ToArray();
    var matQ = gr.Select(e => mat[e.t1, e.t2].Q).ToKMatrix(nbRows);
    Console.WriteLine(matQ);
    Console.WriteLine();
    var matMsg = gr.Select(e => RLWE.DecryptBGV(mat[e.t1, e.t2], sk)).ToKMatrix(nbRows);
    Console.WriteLine(matMsg);
    Console.WriteLine();
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

(RLWECipher ctboot, RLWECipher ctsm) Bootstrapping(RLWECipher ct, RLWECipher pk, RLWECipher rlk, RLWECipher[] skAut,
    (RLWECipher plus, RLWECipher minus)[] brk)
{
    var (pm, t, qL) = pk.PM_T_Q;
    var n = pm.Degree;
    
    // 1. Extract
    var (ct1, ctprep) = CtPrep(ct);
    var extract = RLWE.Extract(ctprep);
    var c = n / 2 - 1;
    var ni = new Rational(InvModPbezbigint(n, qL.Num));
    var f = (2 * c + 1).SeqLazy(-c).Where(j => j != 0).Select(j => ni * ct.Q * j * XpowA(j, pm, t))
        .Aggregate((v0, v1) => v0 + v1).ResModSigned(pm, t);

    var seqBR = new List<RLWECipher>();
    foreach (var ab in extract)
        seqBR.Add(BlindRotateBGV(ab, f, brk, rlk));

    // Step 3. Repacking
    var ctsm = RepackingBGV(n, seqBR.ToArray(), skAut);
    var ctboot = ctsm - ct1;
    return (ctboot, ctsm);
}

void testBR()
{
    RecomputeAllPrimesUpTo(500000);
    var N = 32;
    var n = N / 2;
    var level = n;
    var (pm, sk, t, primes, pk, swks) = RLWE.SetupBGV(N, level);
    var q = primes[0];
    var sp = Rational.NthRoot(swks[pk.Q].skPow[0].Q / pk.Q, level);
    var qL = pk.Q;
    Console.WriteLine($"BGV level = {level}");
    Console.WriteLine($"pm = {pm} T = {t} q = {q} sp = {sp} qL = {qL}");
    Console.WriteLine($"sk = {sk}");
    Console.WriteLine($"pk => {pk.Params}");
    foreach (var swk in swks)
        Console.WriteLine($"swk[{swk.Key}] => {swk.Value.skPow[0].Params}");

    Console.WriteLine();

    var rlk = swks[qL].skPow[2];
    var brk = BRKBGV(sk, pk);
    var seqMods = (level + 1).SeqLazy(1).Select(i => primes.Take(i).Aggregate((pi, pj) => pi * pj)).ToArray();
    var aks = RLWE.LeveledAutMorphSwitchKeyGenBGV(level, pm, sk, t, seqMods, sp);

    var x = pm.X;
    var c = n / 2 - 1;
    var f = (2 * c + 1).SeqLazy(-c).Select(j => j * XpowA(j, pm, t))
        .Aggregate(x.Zero, (acc, v) => acc + v).ResModSigned(pm, t);

    var s = n.SeqLazy().Select(i => sk[i].Signed(t)).ToVec();

    var nbFails = 0;
    for (int k = 0; k < 50; ++k)
    {
        var ai = new Rational(Rng.Next(1, n + 1)).Signed(n);
        var bi = RLWE.GenUnif(n, n).CoefsExtended(n - 1);

        var acc = n * BlindRotateBGV((ai, bi), f, brk, rlk);
        var actual = RLWE.DecryptBGV(acc, sk);

        var u = CheckBR(s.ToArray(), pm, (ai, bi), f, actual, n * t.One);
        if (!u) Console.WriteLine($"fail:{++nbFails}");

        Console.WriteLine();
    }

    Console.WriteLine($"nbFails:{nbFails}");
}

void testBootstrapping()
{
    RecomputeAllPrimesUpTo(500000);
    var N = 16;
    var n = N / 2;
    var level = 8;
    var (pm, sk, t, primes, pk, swks) = RLWE.SetupBGV(N, 17, level);
    var q = primes[0];
    var sp = Rational.NthRoot(swks[pk.Q].skPow[0].Q / pk.Q, level);
    var qL = pk.Q;
    Console.WriteLine($"BGV level = {level}");
    Console.WriteLine($"pm = {pm} T = {t} q = {q} sp = {sp} qL = {qL}");
    Console.WriteLine($"sk = {sk}");
    Console.WriteLine($"pk => {pk.Params}");
    foreach (var swk in swks)
        Console.WriteLine($"swk[{swk.Key}] => {swk.Value.skPow[0].Params}");

    Console.WriteLine();

    var rlk = swks[qL].skPow[2];
    var brk = BRKBGV(sk, pk);
    var seqMods = (level + 1).SeqLazy(1).Select(i => primes.Take(i).Aggregate((pi, pj) => pi * pj)).ToArray();
    var aks = RLWE.LeveledAutMorphSwitchKeyGenBGV(level, pm, sk, t, seqMods, sp);

    for (int k = 0; k < 50; ++k)
    {
        var m = RLWE.GenUnif(n, t);
        var cm = RLWE.EncryptBGV(m, pk); // level qL
        var ct = cm.ModSwitch(q); // level q0

        var (ctboot, ctsm) = Bootstrapping(ct, pk, rlk, aks[qL].autSk, brk);
        Console.WriteLine($"pm = {pm} T = {t} q = {q} sp = {sp} qL = {qL}");
        Console.WriteLine($"pk     {pk.Params}");
        Console.WriteLine($"ct     {ct.Params}");
        Console.WriteLine($"ctboot {ctboot.Params}");
        
        var mboot = RLWE.DecryptBGV(ctboot, sk);
        Console.WriteLine($"m       = {m}");
        Console.WriteLine($"ctboot  = {mboot}");
        Console.WriteLine($"eboot   = {RLWE.ErrorsBGV(ctboot, sk).NormInf()}");
        Console.WriteLine($"emodsw  = {RLWE.ErrorsBGV(ct.ModSwitch(qL), sk).NormInf()}");
        Console.WriteLine();
        if (!m.Equals(mboot))
            throw new();
    }
}

{
    // testBR();
    // f           :x^15 + 2*x^14 + 3*x^13 + 4*x^12 + 5*x^11 + 6*x^10 + 7*x^9 + 7*x^7 + 6*x^6 + 5*x^5 + 4*x^4 + 3*x^3 + 2*x^2 + x
    // f           :[   0,    1,    2,    3,    4,    5,    6,    7,    0,    7,    6,    5,    4,    3,    2,    1]
    // ai          :-5
    // bi          :[   4,   -2,    0,    3,    2,   -7,    1,    6,   -5,    6,    2,    3,   -1,   -1,   -1,    8]
    // blind rotate:6*x^15 + 7*x^14 + 7*x^12 + 6*x^11 + 5*x^10 + 4*x^9 + 3*x^8 + 2*x^7 + x^6 - x^4 - 2*x^3 - 3*x^2 - 4*x - 5
    //             :[  -5,   -4,   -3,   -2,   -1,    0,    1,    2,    3,    4,    5,    6,    7,    0,    7,    6]
    // u= a - <b,s>:5
    // f*X^5       :6*x^15 + 7*x^14 + 7*x^12 + 6*x^11 + 5*x^10 + 4*x^9 + 3*x^8 + 2*x^7 + x^6 - x^4 - 2*x^3 - 3*x^2 - 4*x - 5
    // f*X^5       :[  -5,   -4,   -3,   -2,   -1,    0,    1,    2,    3,    4,    5,    6,    7,    0,    7,    6]
    // success
    // 
    
    testBootstrapping(); // lucky parameters
    // pm = x^8 + 1 T = 17 q = 1361 sp = 12343 qL = 35599582020475615427283552359810689
    // pk     RLWECipher Q:35599582020475615427283552359810689    T:17    PM:x^8 + 1
    // ct     RLWECipher Q:1361    T:17    PM:x^8 + 1
    // ctboot RLWECipher Q:35599582020475615427283552359810689    T:17    PM:x^8 + 1
    // m       = -4*x^7 + 5*x^6 + 8*x^5 - 2*x^4 + 5*x^3 + 5*x^2 - 6*x + 3
    // ctboot  = -4*x^7 + 5*x^6 + 8*x^5 - 2*x^4 + 5*x^3 + 5*x^2 - 6*x + 3
    // eboot   = 5068021424691139538
    // emodsw  = 575452464695417736517441698689073
    // 
}
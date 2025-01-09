using System.Numerics;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.UserGroup.Floats;
using FastGoat.UserGroup.Integers;
using Rq = FastGoat.Structures.VecSpace.KPoly<FastGoat.UserGroup.Integers.Rational>;

namespace FastGoat.UserGroup.Lattice;

public static class FHE
{
    public static bool NoiseMode { get; set; } = true;
    public static void NoiseOn() => NoiseMode = true;
    public static void NoiseOff() => NoiseMode = false;

    public static Rq GenDiscrGauss(int n, double s = 3.0)
    {
        return DistributionExt.DiscreteGaussianSample(n, s).ToKPoly(Rational.KZero());
    }

    public static Rq GenTernary(int n)
    {
        return DistributionExt.DiceSample(n, [-1, 0, 1]).ToKPoly(Rational.KZero());
    }

    public static Rq GenUnif(int n, int q)
    {
        return DistributionExt.DiceSample(n, 0, q - 1).Select(e => new Rational(e)).ToKPoly();
    }

    public static Rq GenUnif(int n, BigInteger q)
    {
        return DistributionExt.DiceSampleBigInt(n, 0, q - 1).Select(e => new Rational(e)).ToKPoly();
    }

    public static Rq GenUnif(int n, Rational q) => GenUnif(n, q.Num);

    public static (Rq pm, Rq sk) SKBGV(int n)
    {
        var pm = FG.QPoly().Pow(n) + 1;
        var sk = 10000.SeqLazy().Select(_ => GenTernary(n))
            .First(s => !s[n - 1].IsZero() && s.Coefs.Count(e => e.IsZero()) <= n / 4);
        return (pm, sk);
    }

    public static RLWECipher PKBGV(Rq pm, Rq sk, Rational t, Rational q)
    {
        var n = pm.Degree;
        var e = GenDiscrGauss(n);
        var c1 = GenUnif(n, q);
        var c0 = (t * e + c1 * sk).ResMod(pm, q);
        return (c0, c1, pm, q);
    }

    public static RLWECipher RLKBGV(Rq pm, Rq sk, Rational t, Rational q)
    {
        var n = pm.Degree;
        var e = GenDiscrGauss(n);
        var c1 = GenUnif(n, q.Num);
        var c0 = (t * e + c1 * sk + sk.Pow(2)).ResMod(pm, q);
        return (c0, c1, pm, q);
    }

    public static (Rq pm, Rq sk, Rational t, Rational q, RLWECipher pk, RLWECipher rlk) KeyGenBGV(int n, int t, int q)
    {
        var (pm, sk) = SKBGV(n);
        var (T, Q) = (new Rational(t), new Rational(q));
        var pk = PKBGV(pm, sk, T, Q);
        var rlk = RLKBGV(pm, sk, T, Q);
        return (pm, sk, T, Q, pk, rlk);
    }

    public static RLWECipher EncryptBGV(Rq m, Rq pm, Rational t, Rational q, RLWECipher pk, bool noise = true)
    {
        var n = pm.Degree;
        var ea = GenDiscrGauss(n);
        var eb = GenDiscrGauss(n);
        var u = GenTernary(n);
        if (!NoiseMode || !noise)
        {
            ea = eb = eb.Zero;
            u = u.One;
        }

        var a = (u * pk.A + m + t * ea).ResMod(pm, q);
        var b = (u * pk.B + t * eb).ResMod(pm, q);
        return (a, b, pm, q);
    }

    public static (RLWECipher csm, RLWECipher cm) EncryptRgswBGV(Rq m, Rq pm, Rational t, Rational q, RLWECipher pk, 
        bool noise = true)
    {
        var cm = EncryptBGV(m, pm, t, q, pk);
        var ct = EncryptBGV(pm.Zero, pm, t, q, pk, noise);
        var csm = (ct.A, (ct.B - m).CoefsMod(q), pm, q);
        return (csm, cm);
    }

    public static Rq DecryptBGV(RLWECipher cipher, Rq pm, Rq sk, Rational t)
    {
        return (cipher.A - sk * cipher.B).ResMod(pm, t);
    }

    public static RLWECipher AddBGV(RLWECipher cipher0, RLWECipher cipher1, Rational q)
    {
        var a = (cipher0.A + cipher1.A).CoefsMod(q);
        var b = (cipher0.B + cipher1.B).CoefsMod(q);
        return (a, b, cipher0.PM, q);
    }

    public static RLWECipher SubBGV(RLWECipher cipher0, RLWECipher cipher1, Rational q)
    {
        var a = (cipher0.A - cipher1.A).CoefsMod(q);
        var b = (cipher0.B - cipher1.B).CoefsMod(q);
        return (a, b, cipher0.PM, q);
    }

    public static RLWECipher KMulBGV(RLWECipher cipher, Rational k, Rational q)
    {
        var a = (cipher.A * k).CoefsMod(q);
        var b = (cipher.B * k).CoefsMod(q);
        return (a, b, cipher.PM, q);
    }

    public static RLWECipher KMulBGV(RLWECipher cipher, Rq k, Rq pm, Rational q)
    {
        var a = (cipher.A * k).ResMod(pm, q);
        var b = (cipher.B * k).ResMod(pm, q);
        return (a, b, pm, q);
    }

    public static RLWECipher MulNaiveBGV(RLWECipher cipher0, RLWECipher cipher1, Rq pm, Rational q)
    {
        var a = (cipher0.A * cipher1.A).ResMod(pm, q);
        var b = (cipher0.B * cipher1.B).ResMod(pm, q);
        return (a, b, pm, q);
    }

    public static RLWECipher MulRelinBGV(RLWECipher cipher0, RLWECipher cipher1, Rq pm, Rational q, RLWECipher rlk)
    {
        var d0 = (cipher0.A * cipher1.A).ResMod(pm, q);
        var d1 = (cipher0.A * cipher1.B + cipher0.B * cipher1.A).ResMod(pm, q);
        var d2 = (cipher0.B * cipher1.B).ResMod(pm, q);

        var a = (d0 + d2 * rlk.A).ResMod(pm, q);
        var b = (d1 + d2 * rlk.B).ResMod(pm, q);
        return (a, b, pm, q);
    }
    
    public static (RLWECipher es2, RLWECipher es) ESKBGV(Rq pm, Rq sk, Rational t, Rational q, RLWECipher pk)
    {
        return (EncryptBGV(sk.Pow(2).ResMod(pm, t), pm, t, q, pk), EncryptBGV(sk, pm, t, q, pk));
    }
    
    public static RLWECipher MulSwkBGV(RLWECipher cm1, RLWECipher cm2, Rq pm, Rational q, RLWECipher es, RLWECipher es2)
    {
        var csm2a = (es.A * cm2.A - es2.A * cm2.B).ResMod(pm, q);
        var csm2b = (es.B * cm2.A - es2.B * cm2.B).ResMod(pm, q);
        var cm1m2a = (cm2.A * cm1.A - csm2a * cm1.B).ResMod(pm, q);
        var cm1m2b = (cm2.B * cm1.A - csm2b * cm1.B).ResMod(pm, q);
        return (cm1m2a, cm1m2b, pm, q);
    }
    
    public static RLWECipher SWKBGV(Rq pm, Rq sk, Rational t, Rational q, RLWECipher pk)
    {
        return EncryptBGV(sk, pm, t, q, pk);
    }

    public static RLWECipher SwitchKeysBGV(RLWECipher cipher, Rq pm, Rational t, Rational q, RLWECipher pk, RLWECipher swk)
    {
        return SubBGV((cipher.A, cipher.A.Zero, pm, q), KMulBGV(swk, cipher.B, pm, q), q);
    }

    public static Dictionary<int, RLWECipher> AKBGV(Rq pm, Rq sk, Rational t, Rational q, RLWECipher pk)
    {
        var n = pm.Degree;
        var x = pm.X;
        return FG.UnInt(2 * n).Order()
            .Select(i => (i.K, EncryptBGV(sk.Substitute(x.Pow(i.K)).ResMod(pm, t), pm, t, q, pk)))
            .ToDictionary(e => e.K, e => e.Item2);
    }

    public static RLWECipher AutoMorphBGV(RLWECipher cipher, int k, Rq pm, Rational t, Rational q, RLWECipher pk,
         Dictionary<int, RLWECipher> ak)
    {
        var xk = pm.X.Pow(k);
        var a = cipher.A.Substitute(xk).ResMod(pm, q);
        var b = cipher.B.Substitute(xk).ResMod(pm, q);
        return SwitchKeysBGV((a, b, pm, q), pm, t, q, pk, ak[k]);
    }

    public static (RLWECipher plus, RLWECipher minus)[] BRKBGV(Rq pm, Rq sk, Rational t, Rational q, RLWECipher pk)
    {
        var enc = (int k) => EncryptBGV(pm.One * k, pm, t, q, pk);
        var n = pm.Degree;
        return n.SeqLazy().Select(i => sk[i])
            .Select(c => (plus: ((c - 1).IsZero() ? enc(1) : enc(0)), minus: ((c + 1).IsZero() ? enc(1) : enc(0))))
            .ToArray();
    }

    public static ((RLWECipher sp, RLWECipher p) plus, (RLWECipher sm, RLWECipher m) minus)[] 
        BRKgswBGV(Rq pm, Rq sk, Rational t, Rational q, RLWECipher pk)
    {
        var enc = (int k) => EncryptRgswBGV(pm.One * k, pm, t, q, pk);
        var n = pm.Degree;
        return n.SeqLazy().Select(i => sk[i])
            .Select(c => (plus: ((c - 1).IsZero() ? enc(1) : enc(0)), minus: ((c + 1).IsZero() ? enc(1) : enc(0))))
            .ToArray();
    }

    public static Rq XpowA(int a, Rq pm, Rational q)
    {
        var x = pm.X;
        if (a == 0)
            return x.One;
        
        var n = pm.Degree;
        var sgn = a > 0 || a % n == 0 ? 0 : 1;
        return (x.Pow(IntExt.AmodP(a, n)) * (-1).Pow(a / n + sgn)).CoefsMod(q);
    }

    public static RLWECipher BlindRotateBGV((Rational ai, Rational[] bi) ab, Rq f, Rq pm, Rational q, RLWECipher pk,
        RLWECipher rlk, (RLWECipher plus, RLWECipher minus)[] brk)
    {
        var n = pm.Degree;
        var x = pm.X;

        var beta = ab.bi;
        var alpha = (int)ab.ai.Num;
        var xalpha = XpowA(alpha, pm, q);
        RLWECipher encOne = (x.One, x.Zero, pm, q);
        RLWECipher acc = ((f * xalpha).ResMod(pm, q), x.Zero, pm, q);

        for (int i = 0; i < n; i++)
        {
            var (encSi_plus, encSi_minus) = brk[i];

            var ai = (int)beta[i].Opp().Mod(q).Num;
            
            var exai = (XpowA(ai, pm, q) - 1).CoefsMod(q);
            var cxai = KMulBGV(encSi_plus, exai, pm, q);
            var ex_ai = (XpowA(-ai, pm, q) - 1).CoefsMod(q);
            var cx_ai = KMulBGV(encSi_minus, ex_ai, pm, q);

            var c0 = (encOne.A + cxai.A + cx_ai.A).ResMod(pm, q);
            var c1 = (encOne.B + cxai.B + cx_ai.B).ResMod(pm, q);
            
            acc = MulRelinBGV(acc, (c0, c1, pm, q), pm, q, rlk);
        }

        return acc;
    }

    public static RLWECipher BlindRotateRgswBGV((Rational ai, Rational[] bi) ab, Rq f, Rq pm, Rational q, RLWECipher pk,
         ((RLWECipher sp, RLWECipher p) plus, (RLWECipher sm, RLWECipher m) minus)[] brk)
    {
        var n = pm.Degree;
        var x = pm.X;

        var beta = ab.bi;
        var alpha = (int)ab.ai.Num;
        var xalpha = XpowA(alpha, pm, q);
        (RLWECipher s, RLWECipher o) encOne = ((pm.Zero, -pm.One, pm, q), (pm.One, pm.Zero, pm, q));
        RLWECipher acc = ((f * xalpha).ResMod(pm, q), x.Zero, pm, q);

        for (int i = 0; i < n; i++)
        {
            var ((encSi_sp, encSi_p), (encSi_sm, encSi_m)) = brk[i];

            var ai = (int)beta[i].Opp().Mod(q).Num;
            
            var exai = (XpowA(ai, pm, q) - 1).CoefsMod(q);
            var cxai_sp = KMulBGV(encSi_sp, exai, pm, q);
            var cxai_p = KMulBGV(encSi_p, exai, pm, q);
            
            var ex_ai = (XpowA(-ai, pm, q) - 1).CoefsMod(q);
            var cx_ai_sm = KMulBGV(encSi_sm, ex_ai, pm, q);
            var cx_ai_m = KMulBGV(encSi_m, ex_ai, pm, q);
            
            var sa = (encOne.s.A + cxai_sp.A + cx_ai_sm.A).ResMod(pm, q);
            var sb = (encOne.s.B + cxai_sp.B + cx_ai_sm.B).ResMod(pm, q);
            var a = (encOne.o.A + cxai_p.A + cx_ai_m.A).ResMod(pm, q);
            var b = (encOne.o.B + cxai_p.B + cx_ai_m.B).ResMod(pm, q);
            
            acc = SubBGV(KMulBGV((a, b, pm, q), acc.A, pm, q), KMulBGV((sa, sb, pm, q), acc.B, pm, q), q);
            
        }

        return acc;
    }
    
    static Rational[] ExtractArr(int i, Rq poly, int n)
    {
        return n.SeqLazy(start: i, step: -1).Select(j => j >= 0 ? poly[j] : -poly[n + j]).ToArray();
    }

    public static (Rational ai, Rational[] bi)[] Extract(RLWECipher e, int n)
    {
        return n.SeqLazy().Select(i => (ai: e.A[i], bi: ExtractArr(i, e.B, n).ToArray())).ToArray();
    }
    
    public static RLWECipher RepackingBGV(RLWECipher[] accs, int n, Rq pm, Rational t, Rational q, RLWECipher pk,
        Dictionary<int, RLWECipher> ak)
    {
        var N = pm.Degree;
        var x = pm.X;
        var cOne = new RLWECipher(x.One, x.Zero, pm, q);
        var CT = Ring.Matrix(cOne, N, N + 1);
        for (int i = 0; i < n; i++)
        {
            CT[i, n] = accs[i];
            for (int j = 1; j < N / n; j++)
            {
                // Console.WriteLine($"[i, n - 1]{(i, n - 1)} [i, j * n]{(i, j * n)}");
                var exnj = XpowA(n * j, pm, t);
                var c = KMulBGV(CT[i, i + j * n], exnj, pm, q);
                CT[i, n] = AddBGV(c, CT[i, n], q);
            }
        }

        for (int k = n; k > 1; k /= 2)
        {
            for (int i = 0; i < k / 2; i++)
            {
                // Console.WriteLine($"[i, k]{(i, k)} [i + k / 2, k]{(i + k / 2, k)} [i, k / 2]{(i, k / 2)}");
                var expowk2 = XpowA(k / 2, pm, t);
                var c0 = KMulBGV(CT[i + k / 2, k], expowk2, pm, q);
                CT[i, k / 2] = AddBGV(CT[i, k], c0, q);

                var c1 = KMulBGV(CT[i + k / 2, k], expowk2, pm, q);
                var crot = AutoMorphBGV(SubBGV(CT[i, k], c1, q), 1 + 2 * N / k, pm, t, q, pk, ak);
                CT[i, k / 2] = AddBGV(CT[i, k / 2], crot, q);
            }
        }

        return CT[0, 1];
    }
    
    public static RLWECipher Bootstrapping(RLWECipher ct, Rq pm, Rational q, Rational Q, RLWECipher pk, RLWECipher rlk,
        (RLWECipher plus, RLWECipher minus)[] brk, Dictionary<int, RLWECipher> ak, int fact = 2)
    {
        var n = pm.Degree;
        var q1 = q / fact;
        if (!q1.IsInteger())
            throw new();
    
        var ct1 = ct.CoefsMod(q1);
        var ctprep = new RLWECipher((ct.A - ct1.A) / q1, (ct.B - ct1.B) / q1, pm, q).CoefsMod(new(fact));
    
        // Step 1. Extraction
        var extract = Extract(ctprep, n);

        var delta = double.Sqrt(n * 4.0);
        var gamma = q / 4 - q1 / 2 * delta;
        var c = (int)(0.5 * (delta + 1) + gamma * q1.Inv());
        var f = (2 * c + 1).SeqLazy(-c).Where(j => j != 0).Select(j => -q1 * j * XpowA(j, pm, q))
            .Aggregate((v0, v1) => v0 + v1).ResMod(pm, q);
    
        // Step 2. Blind Rotate
        var seqBR = new List<RLWECipher>();
        foreach (var ab in extract)
            seqBR.Add(BlindRotateBGV(ab, f, pm, Q, pk, rlk, brk));

        // Step 3. Repacking
        var seqBR0 = seqBR.Select(cipher => new RLWECipher(cipher.A[0] * pm.One, cipher.B[0] * pm.One, pm, q)).ToArray();
        var ctsm = RepackingBGV(seqBR0, n, pm, q1, Q, pk, ak);
        
        return AddBGV(ctsm, ct1, Q);
    }

    public static RLWECipher BootstrappingRgsw(RLWECipher ct, Rq pm, Rational q, Rational Q, RLWECipher pk,
        ((RLWECipher sp, RLWECipher p) plus, (RLWECipher sm, RLWECipher m) minus)[] brk, Dictionary<int, RLWECipher> ak, 
        int fact = 2)
    {
        var n = pm.Degree;
        var q1 = q / fact;
        if (!q1.IsInteger())
            throw new();
    
        var ct1 = ct.CoefsMod(q1);
        var ctprep = new RLWECipher((ct.A - ct1.A) / q1, (ct.B - ct1.B) / q1, pm, q).CoefsMod(new(fact));
    
        // Step 1. Extraction
        var extract = Extract(ctprep, n);

        var delta = double.Sqrt(n * 4.0);
        var gamma = q / 4 - q1 / 2 * delta;
        var c = (int)(0.5 * (delta + 1) + gamma * q1.Inv());
        var f = (2 * c + 1).SeqLazy(-c).Where(j => j != 0).Select(j => -q1 * j * XpowA(j, pm, q))
            .Aggregate((v0, v1) => v0 + v1).ResMod(pm, q);
    
        // Step 2. Blind Rotate
        var seqBR = new List<RLWECipher>();
        foreach (var ab in extract)
            seqBR.Add(BlindRotateRgswBGV(ab, f, pm, Q, pk, brk));

        // Step 3. Repacking
        var seqBR0 = seqBR.Select(cipher => new RLWECipher(cipher.A[0] * pm.One, cipher.B[0] * pm.One, pm, q)).ToArray();
        var ctsm = RepackingBGV(seqBR0, n, pm, q1, Q, pk, ak);
        ct1.Show("ct1");
        ctsm.Show("ctsm");
        
        return AddBGV(ctsm, ct1, Q);
    }

    public static Cplx[] Odd2nthRootsOfUnity(int n)
        => n.SeqLazy().Select(i => new Cplx(Complex.FromPolarCoordinates(1.0, (2 * i + 1) * Double.Pi / n))).ToArray();

    public static Rational InnerProd(Rational[] m0, Rational[] m1, Rational mod)
        => m0.Zip(m1).Aggregate(Rational.KZero(), (acc, e) => (acc + e.First * e.Second).Mod(mod));

    public static RLWECipher Rescale(RLWECipher cipher, Rational q)
    {
        return ((cipher.A / q).RoundPoly(), (cipher.B / q).RoundPoly(), cipher.PM, q);
    }
    public static Rq Signed(Rq poly, Rational q) => poly.Coefs.Select(c => 2 * c > q ? c - q : c).ToKPoly();

    public static double NormCan(Rq poly, Cplx[] roots, Rational q)
    {
        var poly2 = Signed(poly, q);
        return roots.Max(r => poly2.Substitute(r).Magnitude);
    }

    public static double Delta(RLWECipher e, Rq pm, Rational q)
    {
        var nab = NormInf((e.A * e.B).ResMod(pm, q), q);
        var nanb = NormInf(e.A, q) * NormInf(e.B, q);
        return nab / nanb;
    }

    public static double NormCan(RLWECipher e, Cplx[] roots, Rational q) =>
        double.Max(NormCan(e.A, roots, q), NormCan(e.B, roots, q));

    public static double NormInf(Rq poly, Rational q) => poly.Coefs.Select(c => 2 * c > q ? q - c : c).Max();
    public static double NormInf(RLWECipher e, Rational q) => double.Max(NormInf(e.A, q), NormInf(e.B, q));
    public static void Show(Rq pm, Rq sk, Rational t, Rational q, RLWECipher pk, RLWECipher rlk)
    {
        Console.WriteLine($"N = {pm.Degree} t = {t} q = {q}");

        Console.WriteLine("Private Key");
        Console.WriteLine(sk);
        pk.Show("Public Key");
        rlk.Show("Relinearisation Key");
        Console.WriteLine();
    }
}
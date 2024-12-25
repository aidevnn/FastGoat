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
        return DistributionExt.DiscreteGaussianSample(n, 0, s, n / s).ToKPoly(Rational.KZero());
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
        var sk = 10000.SeqLazy().Select(_ => GenTernary(n)).First(s => !s[n - 1].IsZero() && 
                                                                       s.Coefs.Count(e => e.IsZero()) <= n / 4);
        return (pm, sk);
    }

    public static BGVCipher PKBGV(Rq pm, Rq sk, Rational t, Rational q)
    {
        var n = pm.Degree;
        var e = GenDiscrGauss(n);
        var c1 = GenUnif(n, q);
        var c0 = (t * e + c1 * sk).ResMod(pm, q);
        return (c0, c1);
    }

    public static BGVCipher RLKBGV(Rq pm, Rq sk, Rational t, Rational q)
    {
        var n = pm.Degree;
        var e = GenDiscrGauss(n);
        var c1 = GenUnif(n, q.Num);
        var c0 = (t * e + c1 * sk + sk.Pow(2)).ResMod(pm, q);
        return (c0, c1);
    }

    public static (Rq pm, Rq sk, Rational t, Rational q, BGVCipher pk, BGVCipher rlk) KeyGenBGV(int n, int t, int q)
    {
        var (pm, sk) = SKBGV(n);
        var (T, Q) = (new Rational(t), new Rational(q));
        var pk = PKBGV(pm, sk, T, Q);
        var rlk = RLKBGV(pm, sk, T, Q);
        return (pm, sk, T, Q, pk, rlk);
    }

    public static BGVCipher EncryptBGV(Rq m, Rq pm, Rational t, Rational q, BGVCipher pk, bool noise = true)
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
        return (a, b);
    }

    public static (BGVCipher csm, BGVCipher cm) EncryptRgswBGV(Rq m, Rq pm, Rational t, Rational q, BGVCipher pk, 
        bool noise = true)
    {
        var cm = EncryptBGV(m, pm, t, q, pk);
        var ct = EncryptBGV(pm.Zero, pm, t, q, pk, noise);
        var csm = (ct.A, (ct.B - m).CoefsMod(q));
        return (csm, cm);
    }

    public static Rq DecryptBGV(BGVCipher cipher, Rq pm, Rq sk, Rational t)
    {
        return (cipher.A - sk * cipher.B).ResMod(pm, t);
    }

    public static BGVCipher AddBGV(BGVCipher cipher0, BGVCipher cipher1, Rational q)
    {
        var a = (cipher0.A + cipher1.A).CoefsMod(q);
        var b = (cipher0.B + cipher1.B).CoefsMod(q);
        return (a, b);
    }

    public static BGVCipher SubBGV(BGVCipher cipher0, BGVCipher cipher1, Rational q)
    {
        var a = (cipher0.A - cipher1.A).CoefsMod(q);
        var b = (cipher0.B - cipher1.B).CoefsMod(q);
        return (a, b);
    }

    public static BGVCipher KMulBGV(BGVCipher cipher, Rational k, Rational q)
    {
        var a = (cipher.A * k).CoefsMod(q);
        var b = (cipher.B * k).CoefsMod(q);
        return (a, b);
    }

    public static BGVCipher KMulBGV(BGVCipher cipher, Rq k, Rq pm, Rational q)
    {
        var a = (cipher.A * k).ResMod(pm, q);
        var b = (cipher.B * k).ResMod(pm, q);
        return (a, b);
    }

    public static BGVCipher MulNaiveBGV(BGVCipher cipher0, BGVCipher cipher1, Rq pm, Rational q)
    {
        var a = (cipher0.A * cipher1.A).ResMod(pm, q);
        var b = (cipher0.B * cipher1.B).ResMod(pm, q);
        return (a, b);
    }

    public static BGVCipher MulRelinBGV(BGVCipher cipher0, BGVCipher cipher1, Rq pm, Rational q, BGVCipher rlk)
    {
        var d0 = (cipher0.A * cipher1.A).ResMod(pm, q);
        var d1 = (cipher0.A * cipher1.B + cipher0.B * cipher1.A).ResMod(pm, q);
        var d2 = (cipher0.B * cipher1.B).ResMod(pm, q);

        var a = (d0 + d2 * rlk.A).ResMod(pm, q);
        var b = (d1 + d2 * rlk.B).ResMod(pm, q);
        return (a, b);
    }

    public static BGVCipher SWKBGV(Rq pm, Rq sk, Rational t, Rational q, BGVCipher pk)
    {
        return EncryptBGV(sk, pm, t, q, pk);
    }

    public static BGVCipher SwitchKeysBGV(BGVCipher cipher, Rq pm, Rational t, Rational q, BGVCipher pk, BGVCipher swk)
    {
        return SubBGV((cipher.A, cipher.A.Zero), KMulBGV(swk, cipher.B, pm, q), q);
    }

    public static Dictionary<int, BGVCipher> AKBGV(Rq pm, Rq sk, Rational t, Rational q, BGVCipher pk)
    {
        var n = pm.Degree;
        var x = pm.X;
        return FG.UnInt(2 * n).Order()
            .Select(i => (i.K, EncryptBGV(sk.Substitute(x.Pow(i.K)).ResMod(pm, t), pm, t, q, pk)))
            .ToDictionary(e => e.K, e => e.Item2);
    }

    public static BGVCipher AutoMorphBGV(BGVCipher cipher, int k, Rq pm, Rational t, Rational q, BGVCipher pk,
         Dictionary<int, BGVCipher> ak)
    {
        var xk = pm.X.Pow(k);
        var a = cipher.A.Substitute(xk).ResMod(pm, q);
        var b = cipher.B.Substitute(xk).ResMod(pm, q);
        return SwitchKeysBGV((a, b), pm, t, q, pk, ak[k]);
    }

    public static (BGVCipher plus, BGVCipher minus)[] BRKBGV(Rq pm, Rq sk, Rational t, Rational q, BGVCipher pk)
    {
        var enc = (int k) => EncryptBGV(pm.One * k, pm, t, q, pk);
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

    public static BGVCipher BlindRotateBGV((Rational ai, Rational[] bi) ab, Rq f, Rq pm, Rational q, BGVCipher pk,
        BGVCipher rlk, (BGVCipher plus, BGVCipher minus)[] brk)
    {
        var n = pm.Degree;
        var x = pm.X;

        var beta = ab.bi;
        var alpha = (int)ab.ai.Num;
        var xalpha = XpowA(alpha, pm, q);
        BGVCipher encOne = (x.One, x.Zero);
        BGVCipher acc = ((f * xalpha).ResMod(pm, q), x.Zero);

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
            
            acc = MulRelinBGV(acc, (c0, c1), pm, q, rlk);
        }

        return acc;
    }

    static Rational[] ExtractArr(int i, Rq poly, int n)
    {
        return n.SeqLazy(start: i, step: -1).Select(j => j >= 0 ? poly[j] : -poly[n + j]).ToArray();
    }

    public static (Rational ai, Rational[] bi)[] Extract(BGVCipher e, int n)
    {
        return n.SeqLazy().Select(i => (ai: e.A[i], bi: ExtractArr(i, e.B, n).ToArray())).ToArray();
    }
    
    public static BGVCipher RepackingBGV(BGVCipher[] accs, int n, Rq pm, Rational t, Rational q, BGVCipher pk,
        Dictionary<int, BGVCipher> ak)
    {
        var N = pm.Degree;
        var x = pm.X;
        var cOne = new BGVCipher(x.One, x.One);
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
    
    public static BGVCipher Bootstrapping(BGVCipher ct, Rq pm, Rational q, Rational Q, BGVCipher pk, BGVCipher rlk,
        (BGVCipher plus, BGVCipher minus)[] brk, Dictionary<int, BGVCipher> ak, int fact = 2)
    {
        var n = pm.Degree;
        var q1 = q / fact;
        if (!q1.IsInteger())
            throw new();
    
        var ct1 = ct.CoefsMod(q1);
        var ctprep = new BGVCipher((ct.A - ct1.A) / q1, (ct.B - ct1.B) / q1).CoefsMod(new(fact));
    
        // Step 1. Extraction
        var extract = Extract(ctprep, n);

        var delta = double.Sqrt(n * 4.0);
        var gamma = q / 4 - q1 / 2 * delta;
        var c = (int)(0.5 * (delta + 1) + gamma * q1.Inv());
        var f = (2 * c + 1).SeqLazy(-c).Where(j => j != 0).Select(j => -q1 * j * XpowA(j, pm, q))
            .Aggregate((v0, v1) => v0 + v1).ResMod(pm, q);
    
        // Step 2. Blind Rotate
        var seqBR = new List<BGVCipher>();
        foreach (var ab in extract)
            seqBR.Add(BlindRotateBGV(ab, f, pm, Q, pk, rlk, brk));

        // Step 3. Repacking
        var seqBR0 = seqBR.Select(cipher => new BGVCipher(cipher.A[0] * pm.One, cipher.B[0] * pm.One)).ToArray();
        var ctsm = RepackingBGV(seqBR0, n, pm, q1, Q, pk, ak);
        
        return AddBGV(ctsm, ct1, Q);
    }

    public static Cplx[] Odd2nthRootsOfUnity(int n)
        => n.SeqLazy().Select(i => new Cplx(Complex.FromPolarCoordinates(1.0, (2 * i + 1) * Double.Pi / n))).ToArray();

    public static Rational InnerProd(Rational[] m0, Rational[] m1, Rational mod)
        => m0.Zip(m1).Aggregate(Rational.KZero(), (acc, e) => (acc + e.First * e.Second).Mod(mod));

    public static BGVCipher Rescale(BGVCipher cipher, Rational q)
    {
        return ((cipher.A / q).RoundPoly(), (cipher.B / q).RoundPoly());
    }
    public static Rq Signed(Rq poly, Rational q) => poly.Coefs.Select(c => 2 * c > q ? c - q : c).ToKPoly();

    public static double NormCan(Rq poly, Cplx[] roots, Rational q)
    {
        var poly2 = Signed(poly, q);
        return roots.Max(r => poly2.Substitute(r).Magnitude);
    }

    public static double Delta(BGVCipher e, Rq pm, Rational q)
    {
        var nab = NormInf((e.A * e.B).ResMod(pm, q), q);
        var nanb = NormInf(e.A, q) * NormInf(e.B, q);
        return nab / nanb;
    }

    public static double NormCan(BGVCipher e, Cplx[] roots, Rational q) =>
        double.Max(NormCan(e.A, roots, q), NormCan(e.B, roots, q));

    public static double NormInf(Rq poly, Rational q) => poly.Coefs.Select(c => 2 * c > q ? q - c : c).Max();
    public static double NormInf(BGVCipher e, Rational q) => double.Max(NormInf(e.A, q), NormInf(e.B, q));
    public static void Show(Rq pm, Rq sk, Rational t, Rational q, BGVCipher pk, BGVCipher rlk)
    {
        Console.WriteLine($"N = {pm.Degree} t = {t} q = {q}");

        Console.WriteLine("Private Key");
        Console.WriteLine(sk);
        Console.WriteLine("Public Key");
        Console.WriteLine(pk.A);
        Console.WriteLine(pk.B);
        Console.WriteLine("Relinearisation Key");
        Console.WriteLine(rlk.A);
        Console.WriteLine(rlk.B);
        Console.WriteLine();
    }
}
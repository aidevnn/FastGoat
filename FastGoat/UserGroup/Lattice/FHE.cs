using System.Numerics;
using FastGoat.Commons;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Floats;
using FastGoat.UserGroup.Integers;
using Rq = FastGoat.Structures.VecSpace.KPoly<FastGoat.UserGroup.Integers.Rational>;

namespace FastGoat.UserGroup.Lattice;

public static class FHE
{
    public static bool NoiseMode { get; set; } = true;

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
        var c0 = (t * e + c1 * sk).ResMod(pm).CoefsMod(q);
        return (c0, c1);
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

        var a = (u * pk.A + m + t * ea).ResMod(pm).CoefsMod(q);
        var b = (u * pk.B + t * eb).ResMod(pm).CoefsMod(q);
        return (a, b);
    }

    public static Rq DecryptBGV(BGVCipher cipher, Rq pm, Rq sk, Rational t)
    {
        return (cipher.A - sk * cipher.B).ResMod(pm).CoefsMod(t);
    }

    public static (Rq pm, Rq sk, Rational t, Rational q, BGVCipher pk, BGVCipher rlk) KeyGenBGV(int n, int t, int q)
    {
        var (pm, sk) = SKBGV(n);
        var (T, Q) = (new Rational(t), new Rational(q));
        var pk = PKBGV(pm, sk, T, Q);
        var rlk = RLKBGV(pm, sk, T, Q);
        return (pm, sk, T, Q, pk, rlk);
    }

    public static BGVCipher AddBGV(BGVCipher cipher0, BGVCipher cipher1, Rational q)
    {
        var a = (cipher0.A + cipher1.A).CoefsMod(q);
        var b = (cipher0.B + cipher1.B).CoefsMod(q);
        return (a, b);
    }

    public static BGVCipher KAddBGV(BGVCipher cipher, Rq poly, Rational q)
    {
        var a = (cipher.A + poly).CoefsMod(q);
        var b = (cipher.B).CoefsMod(q);
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
        var a = (cipher.A * k).ResMod(pm).CoefsMod(q);
        var b = (cipher.B * k).ResMod(pm).CoefsMod(q);
        return (a, b);
    }

    public static BGVCipher RLKBGV(Rq pm, Rq sk, Rational t, Rational q)
    {
        var n = pm.Degree;
        var e = GenDiscrGauss(n);
        var c1 = GenUnif(n, q.Num);
        var c0 = (t * e + c1 * sk + sk.Pow(2)).ResMod(pm).CoefsMod(q);
        return (c0, c1);
    }

    public static BGVCipher MulBGV(BGVCipher cipher0, BGVCipher cipher1, Rq pm, Rational t, Rational q, BGVCipher rlk)
    {
        var d0 = (cipher0.A * cipher1.A).ResMod(pm).CoefsMod(q);
        var d1 = (cipher0.A * cipher1.B + cipher0.B * cipher1.A).ResMod(pm).CoefsMod(q);
        var d2 = (cipher0.B * cipher1.B).ResMod(pm).CoefsMod(q);

        var a = (d0 + d2 * rlk.A).ResMod(pm).CoefsMod(q);
        var b = (d1 + d2 * rlk.B).ResMod(pm).CoefsMod(q);
        return (a, b);
    }

    public static BGVCipher SWKBGV(Rq pm, Rq sk, Rational t, Rational q, BGVCipher pk)
    {
        return EncryptBGV(sk, pm, t, q, pk);
    }

    public static BGVCipher SwitchKeysBGV(BGVCipher cipher, Rq pm, Rational t, Rational q, BGVCipher pk, BGVCipher rlk,
        BGVCipher swk)
    {
        return SubBGV((cipher.A, cipher.A.Zero), MulBGV(EncryptBGV(cipher.B, pm, t, q, pk), swk, pm, t, q, rlk), q);
    }

    public static Dictionary<int, BGVCipher> AKBGV(Rq pm, Rq sk, Rational t, Rational q, BGVCipher pk)
    {
        var n = pm.Degree;
        var x = pm.X;
        return FG.UnInt(2 * n).Order()
            .Select(i => (i.K, EncryptBGV(sk.Substitute(x.Pow(i.K)).ResMod(pm).CoefsMod(t), pm, t, q, pk)))
            .ToDictionary(e => e.K, e => e.Item2);
    }

    public static BGVCipher AutoMorphBGV(BGVCipher cipher, int k, Rq pm, Rational t, Rational q, BGVCipher pk,
        BGVCipher rlk, Dictionary<int, BGVCipher> ak)
    {
        var xk = pm.X.Pow(k);
        var a = cipher.A.Substitute(xk).ResMod(pm).CoefsMod(q);
        var b = cipher.B.Substitute(xk).ResMod(pm).CoefsMod(q);
        return SwitchKeysBGV((a, b), pm, t, q, pk, rlk, ak[k]);
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

    public static BGVCipher BlindRotateBGV((Rational ai, Rational[] bi) ab, Rq f, Rq pm, Rational t, Rational q,
        BGVCipher pk, BGVCipher rlk, (BGVCipher plus, BGVCipher minus)[] brk)
    {
        var n = pm.Degree;
        var x = pm.X;

        var beta = ab.bi.Select(c => c.Opp().Mod(q)).ToArray();
        var alpha = (int)ab.ai.Num;
        var xalpha = XpowA(alpha, pm, t);
        BGVCipher encOne = (x.One,x.Zero);
        BGVCipher acc = ((f * xalpha).ResMod(pm).CoefsMod(q), x.Zero);

        for (int i = 0; i < n; i++)
        {
            var (encSi_plus, encSi_minus) = brk[i];

            var ai = (int)beta[i].Mod(t).Num;
            
            var exai = (XpowA(ai, pm, t) - 1).CoefsMod(q);
            var cxai = KMulBGV(encSi_plus, exai, pm, q);
            var ex_ai = (XpowA(-ai, pm, t) - 1).CoefsMod(q);
            var cx_ai = KMulBGV(encSi_minus, ex_ai, pm, q);

            var c0 = (encOne.A + cxai.A + cx_ai.A).ResMod(pm).CoefsMod(q);
            var c1 = (encOne.B + cxai.B + cx_ai.B).ResMod(pm).CoefsMod(q);
            acc = MulBGV(acc, (c0, c1), pm, t, q, rlk);
        }

        return acc;
    }

    public static Rq Scale(Rq poly, Rq pm, Rational fact, Rational Q)
    {
        var q = Q / fact;
        if (!q.IsInteger())
            throw new($"q = Q/fact = {q} is not integer");
        
        var delta = poly.CoefsMod(fact);
        return (delta + (poly - delta) / fact).CoefsMod(q);
    }

    public static BGVCipher Scale(BGVCipher cipher, Rq pm, Rational fact, Rational Q)
        => (Scale(cipher.A, pm, fact, Q), Scale(cipher.B, pm, fact, Q));

    static Rational[] ExtractArr(int i, Rq poly, int n)
    {
        return n.SeqLazy(start: i, step: -1).Select(j => j >= 0 ? poly[j] : -poly[n + j]).ToArray();
    }

    public static (Rational ai, Rational[] bi)[] Extract(BGVCipher e, int n)
    {
        return n.SeqLazy().Select(i => (ai: e.A[i], bi: ExtractArr(i, e.B, n).ToArray())).ToArray();
    }

    public static Cplx[] Odd2nthRootsOfUnity(int n)
        => n.SeqLazy().Select(i => new Cplx(Complex.FromPolarCoordinates(1.0, (2 * i + 1) * Double.Pi / n))).ToArray();

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
        var nab = NormInf((e.A * e.B).ResMod(pm).CoefsMod(q), q);
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
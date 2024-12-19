using FastGoat.Commons;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Integers;

namespace FastGoat.UserGroup.Lattice;

public record BGVParams(
    int N,
    Rational P,
    Rational Q,
    KPoly<Rational> PM,
    (KPoly<Rational> pk0, KPoly<Rational> pk1) PK,
    (KPoly<Rational> ek0, KPoly<Rational> ek1) EK,
    Dictionary<int, (KPoly<Rational> ek0, KPoly<Rational> ek1)> AK)
{
    public (int N, Rational P, Rational Q, KPoly<Rational> PM) Params_NPQ_PM => (N, P, Q, PM);
}

public struct BGVPublic
{
    public int N  { get; }
    public Rational P { get; }
    public Rational Q { get; }
    public KPoly<Rational> PM { get; }
    public (KPoly<Rational> pk0, KPoly<Rational> pk1) PK  { get; }

    public (KPoly<Rational> ek0, KPoly<Rational> ek1) EK  { get; }
    public Dictionary<int, (KPoly<Rational> ct0, KPoly<Rational> ct1)> AK { get; }

    public (int N, Rational P, Rational Q, KPoly<Rational> PM) Params_NPQ_PM => (N, P, Q, PM);

    public BGVPublic(
        int n,
        Rational p,
        Rational q,
        KPoly<Rational> pm,
        (KPoly<Rational> pk0, KPoly<Rational> pk1) pk,
        (KPoly<Rational> ek0, KPoly<Rational> ek1) ek,
        Dictionary<int, (KPoly<Rational> ek0, KPoly<Rational> ek1)> ak)
    {
        (N, P, Q, PM) = (n, p, q, pm);
        (PK, EK, AK) = (pk, ek, ak);
    }

    public (KPoly<Rational> ct0, KPoly<Rational> ct1) Encrypt(KPoly<Rational> m) => Encrypt(m, N, P, Q, PM, PK);
    
    public (KPoly<Rational> ct0, KPoly<Rational> ct1) RGSW(KPoly<Rational> m)
    {
        var (ct0, ct1) = Encrypt(m.Zero);
        return (ct0, (-m + ct1).CoefsMod(Q));
    }

    public (KPoly<Rational> ct0, KPoly<Rational> ct1) SwitchKey((KPoly<Rational> ct0, KPoly<Rational> ct1) cipher,
        (KPoly<Rational> ct0, KPoly<Rational> ct1) swk)
    {
        return Sub((cipher.ct0, cipher.ct0.Zero), Mul(Encrypt(cipher.ct1), swk));
    }

    public (KPoly<Rational> ct0, KPoly<Rational> ct1) EvalAuto((KPoly<Rational> ct0, KPoly<Rational> ct1) cipher, int k)
    {
        var xk = PM.X.Pow(k);
        var ct0 = cipher.ct0.Substitute(xk).ResMod(PM).CoefsMod(Q);
        var ct1 = cipher.ct1.Substitute(xk).ResMod(PM).CoefsMod(Q);
        var s1sk = AK[k];
        return SwitchKey((ct0,ct1), s1sk);
    }

    public (KPoly<Rational>ct0, KPoly<Rational>ct1) Add((KPoly<Rational>ct0, KPoly<Rational>ct1) cipher0,
        (KPoly<Rational>ct0, KPoly<Rational>ct1) cipher1)
    {
        var ct0 = (cipher0.ct0 + cipher1.ct0).ResMod(PM).CoefsMod(Q);
        var ct1 = (cipher0.ct1 + cipher1.ct1).ResMod(PM).CoefsMod(Q);
        return (ct0, ct1);
    }
    public (KPoly<Rational>ct0, KPoly<Rational>ct1) Sub((KPoly<Rational>ct0, KPoly<Rational>ct1) cipher0,
        (KPoly<Rational>ct0, KPoly<Rational>ct1) cipher1)
    {
        var ct0 = (cipher0.ct0 - cipher1.ct0).ResMod(PM).CoefsMod(Q);
        var ct1 = (cipher0.ct1 - cipher1.ct1).ResMod(PM).CoefsMod(Q);
        return (ct0, ct1);
    }

    public (KPoly<Rational>ct0, KPoly<Rational>ct1) KMul((KPoly<Rational>ct0, KPoly<Rational>ct1) cipher0, Rational k)
    {
        var ct0 = (cipher0.ct0 * k).ResMod(PM).CoefsMod(Q);
        var ct1 = (cipher0.ct1 * k).ResMod(PM).CoefsMod(Q);
        return (ct0, ct1);
    }

    public (KPoly<Rational>ct0, KPoly<Rational>ct1) Mul((KPoly<Rational>ct0, KPoly<Rational>ct1) cipher0,
        (KPoly<Rational>ct0, KPoly<Rational>ct1) cipher1)
        => Mul(cipher0, cipher1, EK);

    public (KPoly<Rational> ct0, KPoly<Rational> ct1) Pow((KPoly<Rational>ct0, KPoly<Rational>ct1) cipher, int k)
    {
        if (k < 1)
            throw new();

        if (k == 1)
            return cipher;

        return Mul(cipher, Pow(cipher, k - 1));
    }

    public (KPoly<Rational>ct0, KPoly<Rational>ct1) Mul((KPoly<Rational>ct0, KPoly<Rational>ct1) cipher0,
        (KPoly<Rational>ct0, KPoly<Rational>ct1) cipher1, (KPoly<Rational>ek0, KPoly<Rational>ek1) ek)
    {
        var (multi_ct0, multi_ct1, multi_ct2) = Tensor(cipher0, cipher1);
        var (ek0, ek1) = ek;
        var multi_ct20 = (multi_ct0 + multi_ct2 * ek0).ResMod(PM).CoefsMod(Q);
        var multi_ct21 = (multi_ct1 + multi_ct2 * ek1).ResMod(PM).CoefsMod(Q);
        return (multi_ct20, multi_ct21);
    }

    public (KPoly<Rational> d0, KPoly<Rational> d1, KPoly<Rational> d2) 
        Tensor((KPoly<Rational>ct0, KPoly<Rational>ct1) cipher0, (KPoly<Rational>ct0, KPoly<Rational>ct1) cipher1)
    {
        var (m1_ct0, m1_ct1) = cipher0;
        var (m2_ct0, m2_ct1) = cipher1;

        var multi_ct0 = (m1_ct0 * m2_ct0).ResMod(PM).CoefsMod(Q);
        var multi_ct1 = (m1_ct0 * m2_ct1 + m1_ct1 * m2_ct0).ResMod(PM).CoefsMod(Q);
        var multi_ct2 = (m1_ct1 * m2_ct1).ResMod(PM).CoefsMod(Q);
        return (multi_ct0, multi_ct1, multi_ct2);
    }

    public KPoly<Rational> Scale(KPoly<Rational> poly, Rational q)
    {
        if (q > Q)
            throw new();
        
        var delta = (-poly).CoefsMod(q);
        return ((poly + delta) / q).ResMod(PM).CoefsMod(Q);
    }

    public static (KPoly<Rational> ct0, KPoly<Rational> ct1) Encrypt(KPoly<Rational> m, int N, Rational P, Rational Q, 
        KPoly<Rational> PM, (KPoly<Rational> pk0, KPoly<Rational> pk1) PK)
    {
        var e0 = GenDiscrGauss(N, 3.0);
        var e1 = GenDiscrGauss(N, 3.0);
        var u = GenTernary(N);
        
        var ct0 = (PK.pk0 * u + P * e0 + m).ResMod(PM).CoefsMod(Q);
        var ct1 = (PK.pk1 * u + P * e1).ResMod(PM).CoefsMod(Q);
        return (ct0, ct1);
    }

    public static KPoly<Rational> GenDiscrGauss(int n, double s)
    {
        return DistributionExt.DiscreteGaussianSample(n, 0, s, n / s).ToKPoly(Rational.KZero());
    }

    public static KPoly<Rational> GenTernary(int n)
    {
        return DistributionExt.DiceSample(n, [-1, 0, 1]).ToKPoly(Rational.KZero());
    }

    public static KPoly<Rational> GenUnif(int n, int q)
    {
        return DistributionExt.DiceSample(n, 0, q - 1).Select(e => new Rational(e)).ToKPoly();
    }
}
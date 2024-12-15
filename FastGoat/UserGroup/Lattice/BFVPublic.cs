using FastGoat.Commons;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Integers;

namespace FastGoat.UserGroup.Lattice;

public struct BFVPublic
{
    public int N { get; }
    public int P { get; }
    public int T { get; }
    public int Q { get; }
    public KPoly<Rational> PM { get; }
    public (KPoly<Rational> pk0, KPoly<Rational> pk1) PK { get; }

    public (KPoly<Rational> ek0, KPoly<Rational> ek1) EK { get; }

    public Rational FloorQoverP { get; }

    public BFVPublic(int n, int p, int q, int t, KPoly<Rational> pm, (KPoly<Rational> pk0, KPoly<Rational> pk1) pk,
        (KPoly<Rational> ek0, KPoly<Rational> ek1) ek)
    {
        (N, P, Q, T) = (n, p, q, t);
        (PM, PK, EK) = (pm, pk, ek);
        FloorQoverP = new Rational(Q, P).Floor;
    }
    
    public (KPoly<Rational>ct0, KPoly<Rational>ct1) Encrypt(KPoly<Rational> m1)
    {
        var e0 = GenDiscrGauss(N, 3.0);
        var u = GenTernary(N);
        var tmp0 = (PK.pk0 * u + e0 + FloorQoverP * m1).Div(PM).rem;
        var ct0 = ModQ(tmp0).RoundPoly();

        var e1 = GenDiscrGauss(N, 3.0);
        var tmp1 = (PK.pk1 * u + e1).Div(PM).rem;
        var ct1 = ModQ(tmp1);
        return (ct0, ct1);
    }

    public (KPoly<Rational>ct0, KPoly<Rational>ct1) Add((KPoly<Rational>ct0, KPoly<Rational>ct1) cipher0,
        (KPoly<Rational>ct0, KPoly<Rational>ct1) cipher1)
    {
        var ct0 = ModQ((cipher0.ct0 + cipher1.ct0).Div(PM).rem);
        var ct1 = ModQ((cipher0.ct1 + cipher1.ct1).Div(PM).rem);
        return (ct0, ct1);
    }

    public (KPoly<Rational>ct0, KPoly<Rational>ct1) Mul((KPoly<Rational>ct0, KPoly<Rational>ct1) cipher0,
        (KPoly<Rational>ct0, KPoly<Rational>ct1) cipher1)
    {
        var (m1_ct0, m1_ct1) = cipher0;
        var (m2_ct0, m2_ct1) = cipher1;

        var fact = new Rational(P, Q);
        var multi_ct0 = ModQ((fact * m1_ct0 * m2_ct0).Div(PM).rem).RoundPoly();
        var multi_ct1 = ModQ((fact * (m1_ct0 * m2_ct1 + m1_ct1 * m2_ct0)).Div(PM).rem).RoundPoly();
        var multi_ct2 = ModQ((fact * m1_ct1 * m2_ct1).Div(PM).rem).RoundPoly();

        var (ek0, ek1) = EK;
        var multi_ct20 = ModQ((multi_ct2 * ek0 / T).RoundPoly());
        var multi_ct21 = ModQ((multi_ct2 * ek1 / T).RoundPoly());
        return (multi_ct0 + multi_ct20, multi_ct1 + multi_ct21);
    }

    KPoly<Rational> ModQ(KPoly<Rational> poly)=> CoefsMod(poly, Q);
    
    public static KPoly<Rational> CoefsMod(KPoly<Rational> poly, int q)
        => poly.Coefs.Select(c => c - (c / q).Floor * q).ToKPoly();

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
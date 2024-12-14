using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Integers;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");

void firstBFV()
{
    // RngSeed(87456);
    Ring.DisplayPolynomial = MonomDisplay.StarCaret;
    var (n, q) = (34, 424242);
    var bfv = new BFV(n, q);
    bfv.Show();
    for (int i = 0; i < 100; i++)
    {
        var mi = BFVPublic.GenUnif(n, bfv.P);
        var ct = bfv.BFVpublic.Encrypt(mi);
        var mf = bfv.Decrypt(ct);
        Console.WriteLine(mi);
        Console.WriteLine(mf);
        Console.WriteLine();
        if (!(mf - mi).IsZero())
            throw new($"step[{i}]");
    }
}

void HEAddBFV()
{
    // RngSeed(87456);
    Ring.DisplayPolynomial = MonomDisplay.StarCaret;
    var (n, q) = (16, 424242);
    var bfv = new BFV(n, q);
    bfv.Show();
    var bfvPub = bfv.BFVpublic;

    for (int i = 0; i < 100; ++i)
    {
        var m1 = BFVPublic.GenUnif(n, bfv.P);
        var m2 = BFVPublic.GenUnif(n, bfv.P);
        var m1m2 = BFVPublic.CoefsMod((m1 + m2).Div(bfvPub.PM).rem, bfv.P);
        var ct1 = bfvPub.Encrypt(m1);
        var ct2 = bfvPub.Encrypt(m2);

        var decrypt = bfv.Decrypt(bfvPub.Add(ct1, ct2));
        Console.WriteLine($"m1 + m2:{m1m2}");
        Console.WriteLine($"decrypt:{decrypt}");
        Console.WriteLine();
        if (!(decrypt - m1m2).IsZero())
            throw new($"step[{i}]");
    }
}

void HEMulBFV()
{
    IntExt.RngSeed(87456);
    Ring.DisplayPolynomial = MonomDisplay.StarCaret;
    var (n, p, q) = (16, 7, 165000);
    var bfv = new BFV(n, p, q);
    bfv.Show();
    var bfvPub = bfv.BFVpublic;

    for (int i = 0; i < 1000; i++)
    {
        var m1 = BFVPublic.GenUnif(n, bfv.P);
        var m2 = BFVPublic.GenUnif(n, bfv.P);
        var m1m2 = BFVPublic.CoefsMod((m1 * m2).Div(bfvPub.PM).rem, bfv.P);
        var ct1 = bfvPub.Encrypt(m1);
        var ct2 = bfvPub.Encrypt(m2);

        var decrypt = bfv.Decrypt(bfvPub.Mul(ct1, ct2));
        Console.WriteLine($"m1 * m2:{m1m2}");
        Console.WriteLine($"decrypt:{decrypt}");
        Console.WriteLine();
        if (!(decrypt - m1m2).IsZero())
            throw new($"step[{i}]");
    }
}

{
    firstBFV();
    HEAddBFV();
    HEMulBFV();
}

// HomomorphicEncryption
// from
// https://CRAN.R-project.org/package=HomomorphicEncryption
// https://github.com/bquast/HomomorphicEncryption/
public class BFV
{
    public KPoly<Rational> SK { get; }
    public int N { get; }
    public int P { get; }
    public int T { get; }
    public int Q { get; }
    public BFVPublic BFVpublic { get; }

    public BFV(int n, int q) : this(n, n / 2 - 1, q)
    {
    }

    public BFV(int n, int p, int q)
    {
        (N, Q, P) = (n, q, p);
        SK = DistributionExt.DiceSample(N, [-1, 0, 1]).ToKPoly(Rational.KZero());

        var t = T = (int)double.Sqrt(n * p * q); // magik T
        var x = FG.QPoly();
        var pm = x.Pow(n) + 1;
        var e0 = DistributionExt.DiscreteGaussianSample(N, 0, 3.0).ToKPoly(Rational.KZero());
        var pk1 = DistributionExt.DiceSample(N, 0, Q - 1).Select(e => new Rational(e)).ToKPoly();
        var pk0 = BFVPublic.CoefsMod(-(pk1 * SK + e0).Div(pm).rem, Q);
        var e1 = DistributionExt.DiscreteGaussianSample(N, 0, 3.0).ToKPoly(Rational.KZero());
        var ek1 = DistributionExt.DiceSample(N, 0, T * Q - 1).Select(e => new Rational(e)).ToKPoly();
        var ek0 = BFVPublic.CoefsMod(-(ek1 * SK + e1) + T * SK.Pow(2), T * Q);
        BFVpublic = new(n, p, q, t, pm, (pk0, pk1), (ek0, ek1));
    }

    public KPoly<Rational> Decrypt((KPoly<Rational>ct0, KPoly<Rational>ct1) cypher)
    {
        var tmp = cypher.ct1 * SK + cypher.ct0;
        return BFVPublic.CoefsMod((P * tmp.Div(BFVpublic.PM).rem / Q).RoundPoly(), P);
    }

    public void Show()
    {
        Console.WriteLine($"N = {N} Q = {Q} T = {T} P = {P}");
        Console.WriteLine("Private Key");
        Console.WriteLine(SK);
        Console.WriteLine("Public Key");
        Console.WriteLine(BFVpublic.PK.pk0);
        Console.WriteLine(BFVpublic.PK.pk1);
        Console.WriteLine();
        Console.WriteLine("Evaluation Key");
        Console.WriteLine(BFVpublic.EK.ek0);
        Console.WriteLine(BFVpublic.EK.ek1);
        Console.WriteLine();
    }
}

public class BFVPublic
{
    public int N { get; }
    public int P { get; }
    public int T { get; }
    public int Q { get; }
    public KPoly<Rational> PM { get; }
    public (KPoly<Rational> pk0, KPoly<Rational> pk1) PK { get; }

    public (KPoly<Rational> ek0, KPoly<Rational> ek1) EK { get; }

    public Rational FloorQoverP { get; }
    public Rational PoverQ => new(P, Q);

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
        var u = GenTernary(N - 1);
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
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
    for (int i = 0; i < 5000; i++)
    {
        var mi = BFV.GenUnif(n, bfv.P);
        var ct = bfv.Encrypt(mi);
        var mf = bfv.Decrypt(ct);
        Console.WriteLine(mi);
        Console.WriteLine(mf);
        Console.WriteLine();
        if (!(mf - mi).IsZero())
            throw new($"step[{i}] {mf - mi}");
    }
}

void secondBFV()
{
    // RngSeed(87456);
    Ring.DisplayPolynomial = MonomDisplay.StarCaret;
    var (n, q) = (16, 424242);
    var bfv = new BFV(n, q);
    bfv.Show();

    for (int k = 0; k < 1000; ++k)
    {
        var m1 = BFV.GenUnif(n, bfv.P);
        var m2 = BFV.GenUnif(n, bfv.P);
        var m1m2 = bfv.ModP((m1 + m2).Div(bfv.PM).rem);
        var ct1 = bfv.Encrypt(m1);
        var ct2 = bfv.Encrypt(m2);

        var decrypt = bfv.Decrypt(bfv.Add(ct1, ct2));
        Console.WriteLine($"m1 + m2:{m1m2}");
        Console.WriteLine($"decrypt:{decrypt}");
        Console.WriteLine();
        if (!(decrypt - m1m2).IsZero())
            throw new();
    }
}

// void thirdBFV()
{
    // RngSeed(87456);
    Ring.DisplayPolynomial = MonomDisplay.StarCaret;
    var (n, q) = (16, 424242);
    var p = n / 2 - 1;
    var bfv = new BFV(n, q);
    bfv.Show();

    for (int k = 0; k < 1000; ++k)
    {
        var m1 = BFV.GenUnif(n, bfv.P);
        var m2 = BFV.GenUnif(n, bfv.P);
        var m1m2 = bfv.ModP((m1 * m2).Div(bfv.PM).rem);
        var (ct1, ct2) = bfv.Encrypt(m1, m2);
        
        var multi_ct0 = bfv.ModQ((ct1.ct0 * ct2.ct0 * new Rational(p, q)).Div(bfv.PM).rem).RoundPoly();
        var multi_ct1 = bfv.ModQ(((ct1.ct0 * ct2.ct1 + ct1.ct1 * ct2.ct0) * new Rational(p, q)).Div(bfv.PM).rem)
            .RoundPoly();
        var multi_ct2 = bfv.ModQ((ct1.ct1 * ct2.ct1 * new Rational(p, q)).Div(bfv.PM).rem).RoundPoly();

        var s = bfv.SK;
        var decrypt = (multi_ct2 * s.Pow(2)) + multi_ct1 * s + multi_ct0;
        decrypt = bfv.ModQ(decrypt.Div(bfv.PM).rem);
        decrypt = bfv.ModP((new Rational(p, q) * decrypt).RoundPoly());
        Console.WriteLine($"m1 * m2:{m1m2}");
        Console.WriteLine($"decrypt:{decrypt}");
        Console.WriteLine();
        if (!(decrypt - m1m2).IsZero())
            throw new();
    }
}

// HomomorphicEncryption
// from
// https://CRAN.R-project.org/package=HomomorphicEncryption
// https://github.com/bquast/HomomorphicEncryption/
public class BFV
{
    public int N { get; }
    public int P { get; }
    public int Q { get; }
    public KPoly<Rational> PM { get; }
    public KPoly<Rational> SK { get; }
    public (KPoly<Rational> pk0, KPoly<Rational> pk1) PK { get; }
    
    public (KPoly<Rational> ek0, KPoly<Rational> ek1) GenEvalKey { get; }
    public Rational FloorQoverP { get; }
    public Rational PoverQ => new(P, Q);

    public BFV(int n, int q) : this(n, n / 2 - 1, q)
    {
    }
    
    public BFV(int n, int p, int q)
    {
        (N, Q, P) = (n, q, p);
        var x = FG.QPoly();
        PM = x.Pow(n) + 1;
        SK = DistributionExt.DiceSample(N, [-1, 0, 1]).ToKPoly(Rational.KZero());
        var a = DistributionExt.DiceSample(N, 0, Q - 1).ToKPoly(Rational.KZero());
        var e = DistributionExt.DiscreteGaussianSample(N, 0, 3.0).ToKPoly(Rational.KZero());
        Console.WriteLine(new { e });
        var pk = ModQ(-(a * SK + e).Div(PM).rem);
        var ek = ModQ(-(a * SK + e) + SK.Pow(2)).Div(PM).rem;
        PK = (pk, a);
        GenEvalKey = (ek, a);
        FloorQoverP = new Rational(Q, P).Floor;
    }

    public (KPoly<Rational>ct0, KPoly<Rational>ct1) Encrypt(KPoly<Rational> m1)
    {
        var e0 = GenDiscrGauss(N, 3.0);
        var u = GenTernary(N - 1);
        var tmp0 = (PK.pk0 * u + e0 + FloorQoverP * m1).Div(PM).rem;
        // var ct0 = tmp0.Coefs.Select(c => new Rational(AmodPbig(c.Num % Q, Q))).ToKPoly();
        var ct0 = ModQ(tmp0).RoundPoly();

        var e1 = GenDiscrGauss(N, 3.0);
        var tmp1 = (PK.pk1 * u + e1).Div(PM).rem;
        // var ct1 = tmp1.Coefs.Select(c => new Rational(AmodPbig(c.Num % Q, Q))).ToKPoly();
        var ct1 = ModQ(tmp1);
        return (ct0, ct1);
    }

    public ((KPoly<Rational>ct0, KPoly<Rational>ct1) c0, (KPoly<Rational>ct0, KPoly<Rational>ct1) c1) 
        Encrypt(KPoly<Rational> m1, KPoly<Rational> m2)
    {
        var e0 = GenDiscrGauss(N, 3.0);
        var e1 = GenDiscrGauss(N, 3.0);
        var u = GenTernary(N - 1);
        
        var tmp00 = (PK.pk0 * u + e0 + FloorQoverP * m1).Div(PM).rem;
        var ct00 = ModQ(tmp00).RoundPoly();

        var tmp01 = (PK.pk1 * u + e1).Div(PM).rem;
        var ct01 = ModQ(tmp01);
        
        var tmp10 = (PK.pk0 * u + e0 + FloorQoverP * m2).Div(PM).rem;
        var ct10 = ModQ(tmp10).RoundPoly();

        var tmp11 = (PK.pk1 * u + e1).Div(PM).rem;
        var ct11 = ModQ(tmp11);
        return ((ct00, ct01), (ct10, ct11));
    }

    public KPoly<Rational> Decrypt((KPoly<Rational>ct0, KPoly<Rational>ct1) cypher)
    {
        var tmp = cypher.ct1 * SK + cypher.ct0;
        return ModP((P * tmp.Div(PM).rem / Q).RoundPoly());
        // tmp = tmp.Div(PM).rem.Coefs.Select(c => new Rational(P * AmodPbig(c.Num % Q, Q), Q)).ToKPoly().RoundPoly();
        // return ModP(ModQ(tmp.Div(PM).rem * new Rational(P, Q)).RoundPoly());
    }

    public (KPoly<Rational>ct0, KPoly<Rational>ct1) Add((KPoly<Rational>ct0, KPoly<Rational>ct1) cipher0,
        (KPoly<Rational>ct0, KPoly<Rational>ct1) cipher1)
    {
        var ct0 = ModQ((cipher0.ct0 + cipher1.ct0).Div(PM).rem);
        var ct1 = ModQ((cipher0.ct1 + cipher1.ct1).Div(PM).rem);
        return (ct0, ct1);
    }

    public KPoly<Rational> ModQ(KPoly<Rational> poly)
        => poly.Coefs.Select(c => c - (c / Q).Floor * Q).ToKPoly();

    public KPoly<Rational> ModP(KPoly<Rational> poly)
        => poly.Coefs.Select(c => c - ((c / P).Floor * P)).ToKPoly();

    public void Show()
    {
        Console.WriteLine($"N = {N} Q = {Q} P = {P}");
        Console.WriteLine("Private Key");
        Console.WriteLine(SK);
        Console.WriteLine("Public Key");
        Console.WriteLine(PK.pk0);
        Console.WriteLine(PK.pk1);
        Console.WriteLine();
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
        return DistributionExt.DiceSample(n, 0, q - 1).ToKPoly(Rational.KZero());
    }
}
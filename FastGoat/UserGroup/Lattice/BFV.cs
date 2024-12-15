using FastGoat.Commons;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Integers;

namespace FastGoat.UserGroup.Lattice;

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
    public Rational PoverQ => new(P, Q);
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
        return BFVPublic.CoefsMod((PoverQ * tmp.Div(BFVpublic.PM).rem).RoundPoly(), P);
    }

    public void Show(bool showKeys = true)
    {
        Console.WriteLine($"N = {N} Q = {Q} T = {T} P = {P}");
        if (showKeys)
        {
            Console.WriteLine("Private Key");
            Console.WriteLine(SK);
            Console.WriteLine("Public Key");
            Console.WriteLine(BFVpublic.PK.pk0);
            Console.WriteLine(BFVpublic.PK.pk1);
            Console.WriteLine();
            Console.WriteLine("Evaluation Key");
            Console.WriteLine(BFVpublic.EK.ek0);
            Console.WriteLine(BFVpublic.EK.ek1);
        }
        Console.WriteLine();
    }
}
using FastGoat.Commons;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Integers;

namespace FastGoat.UserGroup.Lattice;

public class RLWE
{
    static EPoly<ZnInt> RandRq(EPoly<ZnInt> e0, int degree, int bound)
    {
        var o = e0.KOne;
        var x = e0.Poly.x;
        var coefs = degree.Range().Select(_ => IntExt.Rng.Next(bound) * o).ToArray();
        var r = new KPoly<ZnInt>(x, o, coefs);
        return new(e0.F, r);
    }

    static EPoly<ZnInt> RandNormRq(EPoly<ZnInt> e0, int degree, double sigma)
    {
        var o = e0.KOne;
        var x = e0.Poly.x;
        var coefs = DistributionExt.DiscreteGaussianSample(degree, sigma).Select(e => e * o).ToArray();
        var r = new KPoly<ZnInt>(x, o, coefs);
        return new(e0.F, r);
    }

    public int N { get; }
    public int Q { get; }
    public int Err { get; }
    public double Sigma { get; }
    public (EPoly<ZnInt> a, EPoly<ZnInt> p) PK { get; }
    public (EPoly<ZnInt> a, EPoly<ZnInt> e) EK { get; }
    private EPoly<ZnInt> SK { get; }
    private EPoly<ZnInt> A { get; }
    public RLWE(int n, int q, double sigma = 1.0)
    {
        if (!int.IsPow2(n))
            throw new($"N = {n} must be 2^k");

        if (q % (2 * n) != 1)
            throw new($"Q = {q} must be = 1 mod 2n");
        
        (N, Q, Sigma) = (n, q, sigma);
        Err = (int)(double.Sqrt(Q) / 2.0);

        var x = FG.ZPoly(q);
        A = FG.EPoly(x.Pow(n) + 1, 'A').X;

        SK = RandNormRq(A, N, Sigma);
        var a = RandRq(A, N, Q);
        var ra = RandNormRq(A, N, Sigma);
        var pa = ra - a * SK;
        PK = (a, pa);
        var e = RandRq(A, N, Q);
        var re = RandNormRq(A, N, Sigma);
        var pe = re - e * SK + SK.Pow(2);
        EK = (e, pe);
    }

    public EPoly<ZnInt> Encode(int[] seq)
    {
        if (seq.Length != N)
            throw new();

        var x = A.Poly.x;
        var coefs = seq.Take(N).Reverse().Select(e => e == 0 ? A.KZero : A.KOne * (Q - 1) / 2).ToArray();
        return new(A.F, new KPoly<ZnInt>(x, A.KOne, coefs));
    }

    public int[] Decode(EPoly<ZnInt> a)
    {
        return N.Range().Reverse().Select(i => a[i].K).Select(k => k * 2 > Q ? Q - k : k)
            .Select(i => i < Q / 4 ? 0 : 1).ToArray();
    }

    public (EPoly<ZnInt> c1, EPoly<ZnInt> c2) Encrypt(EPoly<ZnInt> m)
    {
        var e1 = RandNormRq(A, N, Sigma);
        var c1 = PK.a * e1;
        var c2 = PK.p * e1 + m;
        return (c1, c2);
    }

    public EPoly<ZnInt> Decrypt((EPoly<ZnInt> c1, EPoly<ZnInt> c2) d) => d.c2 + d.c1 * SK;

    KPoly<Rational> RewritePolynomial(EPoly<ZnInt> P)
        => new('A', Rational.KZero(),
            P.Poly.Coefs.Select(c => c.K > Q / 2 ? c.K - Q : c.K).Select(c => new Rational(c)).ToArray());
    public void Show(bool showKeys = true)
    {
        var k = int.Log2(N);
        Console.WriteLine($"N = {N} = 2^{k} Q = {Q} Sigma = {Sigma} Err = {Err}");
        if (showKeys)
        {
            Console.WriteLine($"Private Key[{RewritePolynomial(SK)}]");
            Console.WriteLine("Public Key (a,p)");
            Console.WriteLine($"    a = {RewritePolynomial(PK.a)}");
            Console.WriteLine($"    p = {RewritePolynomial(PK.p)}");
        }
        
        Console.WriteLine();
    }
}
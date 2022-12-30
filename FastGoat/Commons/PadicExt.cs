using System.Numerics;

namespace FastGoat.Commons;

public static class PadicExt
{
    public static (int val, BigInteger nb) GetValuation(int p, BigInteger n)
    {
        if (n.IsZero)
            return (0, 0);

        var n0 = n;
        var v = 0;
        while (!n0.IsOne)
        {
            var (q, r) = BigInteger.DivRem(n0, p);
            if (r != 0)
                break;

            ++v;
            n0 = q;
        }

        return (v, n0);
    }

    public static (int val, BigInteger nb) AddValuation(int p, BigInteger n, int val)
    {
        var (v0, n0) = GetValuation(p, n);
        return (v0 + val, n0);
    }
}
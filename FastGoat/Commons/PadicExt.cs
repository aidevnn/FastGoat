using System.Numerics;

namespace FastGoat.Commons;

/// <summary>
/// Extension class for p-adic integers. 
/// </summary>
public static class PadicExt
{
    /// <summary>
    /// Calculates the valuation of a p-adic integer.
    /// </summary>
    /// <param name="p">The prime p.</param>
    /// <param name="n">The p-adic integer to calculate the valuation for.</param>
    /// <returns>A tuple containing the valuation and the remaining p-adic.</returns>
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

    /// <summary>
    /// Adds a given valuation to a p-adic integer.
    /// </summary>
    /// <param name="p">The prime p.</param>
    /// <param name="n">The p-adic integer to add the value to.</param>
    /// <param name="val">The valuation to add.</param>
    /// <returns>A tuple containing the sum of the new valuation and the remaining p-adic integer.</returns>
    public static (int val, BigInteger nb) AddValuation(int p, BigInteger n, int val)
    {
        var (v0, n0) = GetValuation(p, n);
        return (v0 + val, n0);
    }
}
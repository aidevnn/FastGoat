using System.Numerics;
using FastGoat.UserGroup.Integers;

namespace FastGoat.UserGroup.Polynoms;

public readonly struct Digits
{
    public int O { get; }
    public BigInteger Max { get; }
    public BigInteger Max4 { get; }

    public Digits(int o)
    {
        if (o < 3)
            throw new($"Nb digits must be > 3");

        O = o;
        Max = BigInteger.Pow(10, O);
        Max4 = BigInteger.Pow(10, O + 4);
    }
    
    public override int GetHashCode() => O;
    
    public BigInteger Clamp(BigInteger k)
    {
        if (k >= Max)
        {
            var n1 = Length(k) - O;
            var k1 = k / BigInteger.Pow(10, n1);
            return k1;
        }
        
        if (k <= -Max)
        {
            var n1 = Length(-k) - O;
            var k1 = k / BigInteger.Pow(10, n1);
            return k1;
        }
        
        return k;
    }

    public static int Length(BigInteger k) => k.IsZero ? 1 : (int)Double.Round(BigInteger.Log10(BigInteger.Abs(k)), 13) + 1;

    // public static (int, BigInteger) TrailingZero(BigInteger k)
    // {
    //     var n0 = BigInteger.TrailingZeroCount(k);
    //     var k0 = k;
    //     int i = 0;
    //     for (; i < n0; ++i)
    //     {
    //         var k1 = BigInteger.DivRem(k0, 10);
    //         if(!k1.Remainder.IsZero)
    //             break;
    //
    //         k0 = k1.Quotient;
    //     }
    //
    //     return (i, k0);
    // }
}
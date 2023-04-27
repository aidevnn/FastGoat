using System.Numerics;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.UserGroup.Integers;

namespace FastGoat.UserGroup.Polynoms;

public readonly struct BigReal : IElt<BigReal>, IRingElt<BigReal>, IFieldElt<BigReal>
{
    public static bool Debug = false;
    public Digits Digit { get; }
    public BigInteger K { get; }
    public int NbDigits { get; }
    public int V { get; }
    public int O => Digit.O;

    public (BigInteger K, int V, int NbDigits, int O) Details => (K, V, NbDigits, O);

    private BigReal(BigInteger k, int v, Digits d)
    {
        K = k;
        V = v;
        if (k.IsZero)
            V = 0;
        
        NbDigits = Digits.Length(K);
        Digit = d;
        Hash = (K, V, Digits: Digit).GetHashCode();
    }
    
    public bool Equals(BigReal other)
    {
        if (V != other.V)
            return false;

        if (K == other.K)
            return true;

        if (NbDigits == other.NbDigits)
            return false;

        if (NbDigits > other.NbDigits)
        {
            var k0 = other.K * BigInteger.Pow(10, NbDigits - other.NbDigits);
            return K == k0;
        }
        else
        {
            var k0 = K * BigInteger.Pow(10, other.NbDigits - NbDigits);
            return other.K == k0;
        }
    }

    public int CompareTo(BigReal other)
    {
        var compS = K.Sign.CompareTo(other.K.Sign);
        if (compS != 0)
            return compS;

        var compV = V.CompareTo(other.V);
        if (compV != 0)
            return compV;

        if (NbDigits == other.NbDigits)
            return K.CompareTo(other.K);
        
        if (NbDigits > other.NbDigits)
        {
            var k0 = other.K * BigInteger.Pow(10, NbDigits - other.NbDigits);
            return K.CompareTo(k0);
        }
        else
        {
            var k0 = K * BigInteger.Pow(10, other.NbDigits - NbDigits);
            return k0.CompareTo(other.K);
        }
    }

    public int Hash { get; }
    public bool IsZero() => K == 0;

    public BigReal Zero => new(0, 0, Digit);
    public BigReal One => new(1, 0, Digit);
    public int P => 0;
    public bool Invertible() => !IsZero();

    public BigReal Pow10(int n) => new(1, n, Digit);
    public BigReal Mul10PowN(int N) => new(K, V + N, Digit);

    public BigReal Add(BigReal e)
    {
        if (O < e.O)
            return ToBigReal(e.O).Add(e);
        
        var d = V - e.V;
        var ad = Int32.Abs(d);
        if (ad > O)
        {
            if (Debug)
                Console.WriteLine(new { cs = 3, e0 = this, d0 = Details, e1 = e, d1 = e.Details });
            if (d < 0)
                return e;
            else
                return this;
        }
        else
        {
            if (d >= 0)
            {
                if (NbDigits >= e.NbDigits)
                {
                    var shift0 = BigInteger.Pow(10, d);
                    var shift1 = BigInteger.Pow(10, NbDigits - e.NbDigits);
                    var k0 = e.K * shift1 + shift0 * K;
                    var k1 = Digit.Clamp(k0);
                    var a = d + NbDigits - Digits.Length(k0);
                    if (Debug)
                        Console.WriteLine(
                            new { cs = 0, e0 = this, d0 = Details, e1 = e, d1 = e.Details, shift0, shift1, k0, k1, V, a });
                    return new(k1, V - a, Digit);
                }
                else
                {
                    var shift0 = BigInteger.Pow(10, d);
                    var shift1 = BigInteger.Pow(10, e.NbDigits - NbDigits);
                    var k0 = e.K + shift1 * shift0 * K;
                    var k1 = Digit.Clamp(k0);
                    var a = d + e.NbDigits - Digits.Length(k0);
                    if (Debug)
                        Console.WriteLine(
                            new { cs = 1, e0 = this, d0 = Details, e1 = e, d1 = e.Details, shift0, shift1, k0, k1, V, a });
                    return new(k1, V - a, Digit);
                }
            }
            else
            {
                if (Debug)
                    Console.WriteLine(new { cs = 2, e0 = this, d0 = Details, e1 = e, d1 = e.Details });
                
                return e.Add(this);
            }
        }
    }

    public BigReal Opp() => new(-K, V, Digit);

    public BigReal Sub(BigReal e) => Add(e.Opp());

    public BigReal Mul(BigReal e)
    {
        if (IsZero() || e.IsZero())
            return Zero;
        
        if (O < e.O)
            return ToBigReal(e.O).Mul(e);

        var k0 = K * e.K;
        var n0 = V + e.V;
        var a = Digits.Length(k0) - NbDigits - e.NbDigits + 1;
        var k1 = Digit.Clamp(k0);
        return new(k1, n0 + a, Digit);
    }

    public (BigReal quo, BigReal rem) Div(BigReal e)
    {
        if (e.IsZero())
            throw new DivideByZeroException();
        
        if (O < e.O)
            return ToBigReal(e.O).Div(e);

        var k0 = (K * Digit.Max4) / e.K;
        var d = Digits.Length(k0) - (O + NbDigits - e.NbDigits + 5);
        var k1 = Digit.Clamp(k0);
        var n1 = V - e.V;
        return (new(k1, n1 + d, Digit), Zero);
    }

    public BigReal Mul(int k)
    {
        var e = FromBigInteger(k, O);
        return Mul(e);
    }

    public BigReal Inv() => One.Div(this).quo;

    public BigReal Pow(int k)
    {
        if (k == 0)
            return One;

        if (k < 0)
            return Inv().Pow(-k);

        var br = this;
        return Enumerable.Repeat(br, k).Aggregate(br.One, (acc, b) => acc * b);
    }

    public override int GetHashCode() => Hash;

    public BigReal ToBigReal(int O)
    {
        var digit = new Digits(O);
        var k0 = digit.Clamp(K);
        return new(k0, V, Digit);
    }

    public double ToDouble
    {
        get
        {
            if (IsZero())
                return 0.0;

            var k0 = DigitsDouble.Clamp(K);
            var d0 = Digits.Length(k0);
            return (double)k0 * Double.Pow(10, V - d0 + 1);
        }
    }

    public Rational ToRational
    {
        get
        {
            if (V - NbDigits + 1 >= 0)
            {
                var num = K * BigInteger.Pow(10, V - NbDigits + 1);
                return new(num);
            }
            else
            {
                var num = K;
                var denom = BigInteger.Pow(10, NbDigits - V - 1);
                return new(num, denom);
                
            }
        }
    }

    public override string ToString()
    {
        if (IsZero())
            return Enumerable.Repeat("0", O - 1).Prepend("0.").Append("E+000").Glue();
        
        var v = V < 0 ? $"-{-V:000}" : $"+{V:000}";
        var s0 = K < 0 ? $"{-K}" : $"{K}";
        var s1 = $"{s0[0]}.{s0.Skip(1).Concat(Enumerable.Repeat('0', O)).Take(O - 1).Glue()}";
        return K < 0 ? $"-{s1}E{v}" : $"{s1}E{v}";
    }

    public static BigReal Pow10(int n, int digit) => new(1, n, new Digits(digit));

    public static BigReal FromDouble(double d, int o)
    {
        var s = $"{d:E}".Split('E');
        var s10 = s[0].Replace(".", "");
        var k0 = long.Parse(s10);
        var sgn = s[1][0] == '+' ? 1 : -1;
        var exp = int.Parse(s[1].Skip(1).Glue()) * sgn;
        var digit = new Digits(o);
        var k1 = digit.Clamp(k0);
        return new(k1, exp, digit);
    }

    public static BigReal FromString(string s, int o)
    {
        if (!s.Contains('E'))
            throw new("Invalid format, scientific format ex:1.2345E+004");

        var s1 = s.Split('E');

        if (!s1[0].Contains('.') || s1[0].ToList().FindIndex(c => c == '.') != 1)
            throw new("Invalid format, scientific format ex:1.2345E+004");

        var s2 = s1[0].Replace(".", "");
        var k = BigInteger.Parse(s2);
        var digit = new Digits(o);
        var k1 = digit.Clamp(k);
        var v = int.Parse(s1[1]);
        return new(k1, v, digit);
    }

    public static BigReal FromBigInteger(BigInteger r, int o)
    {
        var digit = new Digits(o);
        var d = Digits.Length(r) - 1;
        var k = digit.Clamp(r);
        return new(k, d, digit);
    }

    public static BigReal FromRational(Rational r, int o)
    {
        var num = FromBigInteger(r.Num, o);
        var denom = FromBigInteger(r.Denom, o);
        return num / denom;
    }

    public static BigReal operator +(BigReal a, BigReal b) => a.Add(b);

    public static BigReal operator +(int a, BigReal b) => b.One.Mul(a).Add(b);

    public static BigReal operator +(BigReal a, int b) => a.Add(a.One.Mul(b));

    public static BigReal operator -(BigReal a) => a.Opp();

    public static BigReal operator -(BigReal a, BigReal b) => a.Sub(b);

    public static BigReal operator -(int a, BigReal b) => b.One.Mul(a).Sub(b);

    public static BigReal operator -(BigReal a, int b) => a.Sub(a.One.Mul(b));

    public static BigReal operator *(BigReal a, BigReal b) => a.Mul(b);

    public static BigReal operator *(int a, BigReal b) => b.Mul(a);

    public static BigReal operator *(BigReal a, int b) => a.Mul(b);

    public static BigReal operator /(BigReal a, BigReal b) => a.Div(b).quo;

    public static BigReal operator /(BigReal a, int b) => a.Div(a.One.Mul(b)).quo;

    public static BigReal operator /(int a, BigReal b) => b.One.Mul(a).Div(b).quo;

    private static Digits DigitsDouble = new(16);

    public static double Abs(BigReal t) => Double.Abs(t.ToDouble);

    public static bool IsValuedField => true;
}
using System.Numerics;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;

namespace FastGoat.UserGroup.Integers;

public readonly struct Rational : IElt<Rational>, IRingElt<Rational>, IFieldElt<Rational>
{
    public int P => 0;
    public BigInteger Num { get; }
    public BigInteger Denom { get; }
    public static Rational KZero() => new(0, 1);

    public static Rational KOne() => new(1, 1);

    public static double Abs(Rational r) => Absolute(r);
    public static bool IsValuedField => true;
    public Rational()
    {
        Num = 1;
        Denom = 1;
        Hash = (Num, Denom).GetHashCode();
    }

    public Rational(BigInteger num)
    {
        Num = num;
        Denom = 1;
        Hash = (Num, Denom).GetHashCode();
    }

    public Rational(BigInteger num, BigInteger denom)
    {
        if (denom == 0)
            throw new GroupException(GroupExceptionType.GroupDef);

        if (denom.IsOne)
        {
            Num = num;
            Denom = denom;
            Hash = (Num, Denom).GetHashCode();
        }
        else
        {
            var gcd = BigInteger.GreatestCommonDivisor(num, denom);
            Num = denom > 0 ? num / gcd : -num / gcd;
            Denom = denom > 0 ? denom / gcd : -denom / gcd;
            Hash = (Num, Denom).GetHashCode();
        }
    }

    public void Deconstruct(out BigInteger num, out BigInteger denom)
    {
        num = Num;
        denom = Denom;
    }

    public double Abs() => Abs(this);

    public bool Equals(Rational other) => Hash == other.Hash;

    public int CompareTo(Rational other) => (Num * other.Denom).CompareTo(other.Num * Denom);

    public int Hash { get; }
    public bool IsZero() => Num == 0;
    public Rational Zero => new(0, 1);
    public Rational One => new(1, 1);
    public Rational Add(Rational e) => new(Num * e.Denom + e.Num * Denom, Denom * e.Denom);
    public Rational Sub(Rational e) => new(Num * e.Denom - e.Num * Denom, Denom * e.Denom);

    public Rational Opp() => new(-Num, Denom);

    public Rational Mul(Rational e) => new(Num * e.Num, Denom * e.Denom);
    public Rational Mul(int k) => new(Num * k, Denom);

    public Rational Pow(int k)
    {
        if (k == 0)
            return One;

        if (k < 0)
            return Inv().Pow(-k);

        var num = BigInteger.Pow(Num, k);
        var denom = BigInteger.Pow(Denom, k);
        return new(num, denom);
    }

    public (Rational quo, Rational rem) Div(Rational e)
    {
        if (e.IsZero())
            throw new DivideByZeroException();

        return (new(Num * e.Denom, Denom * e.Num), new(0));
    }

    public Rational Inv() => new(Denom, Num);
    public bool Invertible() => !IsZero();
    public bool IsInteger() => Denom.IsOne;
    public override int GetHashCode() => Hash;
    public override string ToString()
    {
        if (IsZero())
            return "0";
        else
        {
            if (Ring.DisplayPolynomial != MonomDisplay.StarPowFct)
                return Denom == 1 ? $"{Num}" : $"{Num}/{Denom}";
            else
            {
                if (BigInteger.Abs(Num) > int.MaxValue || BigInteger.Abs(Denom) > int.MaxValue)
                    return Denom == 1 ? $"{Num}*r1" : $"{Num}*r1/{Denom}";

                return Denom == 1 ? $"{Num}" : $"{Num}/{Denom}";
            }
        }
    }

    public Rational RoundEven
    {
        get
        {
            var (q, r) = BigInteger.DivRem(Num, Denom);
            var rs = r.Sign;
            var r0 = r * rs * 2;
            if (r0 < Denom || (r0 == Denom && BigInteger.IsEvenInteger(q)))
                return new(q, 1);

            return new(q + rs, 1);
        }
    }
    public static implicit operator double(Rational e) => (double)e.Num / (double)e.Denom;
    public static Rational Absolute(Rational e) => new(BigInteger.Abs(e.Num), e.Denom);
    public static Rational operator +(Rational a, Rational b) => a.Add(b);
    public static Rational operator +(int a, Rational b) => b.Add(b.One.Mul(a));
    public static Rational operator +(Rational a, int b) => a.Add(a.One.Mul(b));
    public static Rational operator -(Rational a) => a.Opp();
    public static Rational operator -(Rational a, Rational b) => a + (-b);
    public static Rational operator -(int a, Rational b) => a + (-b);
    public static Rational operator -(Rational a, int b) => a + (-b);
    public static Rational operator *(Rational a, Rational b) => a.Mul(b);
    public static Rational operator *(Rational a, int b) => a.Mul(b);
    public static Rational operator *(int a, Rational b) => b.Mul(a);
    public static Rational operator *(BigInteger a, Rational b) => b.Mul(new Rational(a));
    public static Rational operator /(Rational a, Rational b) => a.Div(b).quo;
    public static Rational operator /(Rational a, int b) => a.Div(a.One.Mul(b)).quo;
    public static Rational operator /(Rational a, BigInteger b) => a.Div(a.One.Mul(new Rational(b))).quo;
    public static Rational operator /(int a, Rational b) => b.Inv().Mul(a);
}
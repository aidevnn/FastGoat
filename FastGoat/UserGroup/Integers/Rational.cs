using System.Numerics;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Floats;

namespace FastGoat.UserGroup.Integers;

public readonly struct Rational : IElt<Rational>, IRingElt<Rational>, IFieldElt<Rational>, IFloatElt<Rational>
{
    public static bool IsComplex => false;
    public int P => 0;
    public BigInteger Num { get; }
    public BigInteger Denom { get; }
    public int Sign => Num.Sign;

    public static Rational KZero() => new(0, 1);

    public static Rational KOne() => new(1, 1);

    public static double Abs(Rational r) => r.Absolute;
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

    public bool IsSquare
    {
        get
        {
            if (!Denom.IsOne)
                return false;

            var dec = IntExt.PrimesDec(Num);
            return dec.Values.All(e => e % 2 == 0);
        }
    }

    public bool Equals(Rational other) => Num * other.Denom == Denom * other.Num;

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
    public double ToDouble => (double)Num / (double)Denom;
    public bool Invertible() => !IsZero();
    public bool IsInteger() => Denom.IsOne;

    public Rational Floor
    {
        get
        {
            if (Denom.IsOne)
                return this;
            
            var div = BigInteger.Divide(Num, Denom);
            var shift = Sign == 1 ? 0 : -1;
            return new(div + shift);
        }
    }

    public (Rational quo, Rational rem) QuoPlusFrac
    {
        get
        {
            var div = BigInteger.Divide(Num, Denom);
            var f = new Rational(Num - div * Denom, Denom);
            return (new(div), f);
        }
    }

    public (Rational quo, Rational rem) FloorPlusFrac
    {
        get
        {
            var i = Floor.Num;
            var f = new Rational(Num - i * Denom, Denom);
            return (new(i), f);
        }
    }

    public Rational Mod(Rational Q) => Sub(Mul(Q.Inv()).Floor * Q);
    public Rational Mod(int q) => Mod(new Rational(q)); 
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

    public Rational Trunc
    {
        get
        {
            var (q, _) = BigInteger.DivRem(Num.Sign * Num, Denom);
            return new(Num.Sign * q);
        }
    }

    public Rational Absolute => new(BigInteger.Abs(Num), Denom);
    public Rational Absolute2 => this * this;
    public Rational Conj => this;

    public static Rational Parse(string r)
    {
        if (r.Any(c => !"0123456789./-".Contains(c)))
            throw new($"Unable to parse r={r}, must contains only 0123456789./-");

        var p = r.Split('/');
        var num = BigInteger.Parse(p[0]);
        if (p.Length == 1)
            return new(num);

        var denom = BigInteger.Parse(p[1]);
        return new(num, denom);
    }

    public static Rational Pi(int o) => BigReal.Pi(o).ToRational;

    public static Rational Sqrt(Rational r) => NthRoot(r, 2);

    public static Rational Sqrt(Rational r, int nbDigits)
    {
        var nb = new Rational(BigInteger.Pow(10, nbDigits));
        var nb2 = nb * nb;
        return Sqrt(r * nb2) / nb;
    }
    

    public static Rational NthRoot(Rational r, int n)
    {
        if (n < 0)
            return NthRoot(r.Inv(), -n);
        
        var num = IntExt.NthRootBigint(r.Num, n);
        var denom = IntExt.NthRootBigint(r.Denom, n);
        return new(num, denom);
    }
    public static implicit operator double(Rational e) => e.ToDouble;
    public static implicit operator Rational(string e) => Rational.Parse(e);
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

    public static Rational Min(Rational a, Rational b) => a > b ? b : a;
    public static Rational Max(Rational a, Rational b) => a > b ? a : b;
}
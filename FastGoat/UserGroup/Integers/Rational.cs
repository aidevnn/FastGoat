using FastGoat.Commons;
using FastGoat.Structures;

namespace FastGoat.UserGroup.Integers;

public struct Rational : IElt<Rational>, IRingElt<Rational>, IFieldElt<Rational>
{
    public int P => 0;
    public int Num { get; }
    public int Denom { get; }

    public static Rational KZero() => new Rational(0, 1);

    public Rational()
    {
        Num = 1;
        Denom = 1;
        Hash = (Num, Denom).GetHashCode();
    }

    public Rational(int num)
    {
        Num = num;
        Denom = 1;
        Hash = (Num, Denom).GetHashCode();
    }

    public Rational(int num, int denom)
    {
        if (denom == 0)
            throw new GroupException(GroupExceptionType.GroupDef);

        var gcd = IntExt.Gcd(num, denom);
        Num = denom > 0 ? num / gcd : -num / gcd;
        Denom = denom > 0 ? denom / gcd : -denom / gcd;
        Hash = (Num, Denom).GetHashCode();
    }

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

        var num = (int)Math.Pow(Num, k);
        var denom = (int)Math.Pow(Denom, k);
        return new(num, denom);
    }

    public (Rational quo, Rational rem) Div(Rational e)
    {
        if (e.IsZero())
            throw new DivideByZeroException();

        return (new(Num * e.Denom, Denom * e.Num), new(0));
    }

    public Rational Inv() => new(Denom, Num);

    public override int GetHashCode() => Hash;
    public override string ToString() => IsZero() ? "0" : Denom == 1 ? $"{Num}" : $"{Num}/{Denom}";

    public static Rational operator +(Rational a, Rational b) => a.Add(b);
    public static Rational operator -(Rational a, Rational b) => a.Sub(b);
    public static Rational operator *(Rational a, Rational b) => a.Mul(b);
    public static Rational operator /(Rational a, Rational b) => a.Div(b).quo;
    public static Rational operator /(Rational a, int b) => new Rational(a.Num, a.Denom * b);
    public static Rational operator *(int a, Rational b) => new Rational(a * b.Num, b.Denom);
}
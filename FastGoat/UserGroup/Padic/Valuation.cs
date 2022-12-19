using System.Numerics;
using FastGoat.Commons;

namespace FastGoat.UserGroup.Padic;

public readonly struct Valuation : IEquatable<Valuation>, IComparable<Valuation>,
    IComparisonOperators<Valuation, Valuation, bool>
{
    private readonly int v;
    private readonly bool isInfinity;
    public static Valuation Infinity => new();

    public Valuation()
    {
        isInfinity = true;
        v = 0;
    }

    public Valuation(int v)
    {
        isInfinity = false;
        this.v = v;
    }

    public bool IsInfinity => isInfinity;
    public int V => !isInfinity ? v : throw new ArgumentException("Infinite valuation");

    public static Valuation operator +(Valuation a, Valuation b) =>
        a.IsInfinity || b.IsInfinity ? new() : new(a.v + b.v);

    public static Valuation operator -(Valuation a, Valuation b) =>
        a.IsInfinity || b.IsInfinity ? new() : new(a.v - b.v);

    public bool Equals(Valuation other) =>
        (IsInfinity && other.IsInfinity) || (!IsInfinity && !other.IsInfinity && v == other.v);

    public int CompareTo(Valuation other)
    {
        if (Equals(other))
            return 0;

        if (IsInfinity)
            return 1;

        if (other.IsInfinity)
            return -1;

        return v.CompareTo(other.v);
    }

    public override int GetHashCode() => (IsInfinity, v).GetHashCode();
    public override string ToString() => IsInfinity ? "∞" : $"{v}";

    public static (Valuation, BigInteger) Of(int p, BigInteger n)
    {
        if (n.IsZero)
            return (new(), n);

        var (v, n0) = PadicExt.GetValuation(p, n);
        return (new(v), n0);
    }

    public static bool operator ==(Valuation left, Valuation right) => left.Equals(right);

    public static bool operator !=(Valuation left, Valuation right) => !left.Equals(right);

    public static bool operator >(Valuation left, Valuation right) => left.CompareTo(right) == 1;

    public static bool operator >=(Valuation left, Valuation right) => left.CompareTo(right) >= 0;

    public static bool operator <(Valuation left, Valuation right) => left.CompareTo(right) == -1;

    public static bool operator <=(Valuation left, Valuation right) => left.CompareTo(right) <= 0;
}
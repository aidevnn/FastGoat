using System.Numerics;

namespace FastGoat.UserGroup.Padic;

public readonly struct Modulus : IEquatable<Modulus>
{
    public int P { get; }
    public int O { get; }
    public BigInteger Mod { get; }

    public Modulus()
    {
        P = 2;
        O = 1;
        Mod = 2;
    }

    public Modulus(int p, int o)
    {
        if (p < 2 || o < 1)
            throw new ArgumentException($"p={p} must be greater than 2 and o={o} must be greater than 1.");

        P = p;
        O = o;
        Mod = BigInteger.Pow(P, O);
    }

    public bool Equals(Modulus other) => P == other.P && O == other.O;

    public static Modulus operator ++(Modulus i) => new(i.P, i.O + 1);
    public static Modulus operator *(Modulus i, int k) => k > 0 ? new(i.P, i.O * k) : throw new();
    public ZnBInt Zero => new(this, 0);
    public override int GetHashCode() => (P, O).GetHashCode();

    public override string ToString() => $"{P}^{O} = {Mod}";
}
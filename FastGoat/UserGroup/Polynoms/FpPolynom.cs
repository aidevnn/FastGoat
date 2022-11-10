using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.UserGroup.Integers;

namespace FastGoat.UserGroup.Polynoms;

public struct FpPolynom : IElt<FpPolynom>
{
    public int[] Coefs { get; }
    public int P { get; }
    private char X { get; }

    public FpPolynom(int p)
    {
        if (!IntExt.Primes10000.Contains(p))
            throw new GroupException(GroupExceptionType.GroupDef);

        P = p;
        Coefs = new[] { 0, 1 };
        X = 'x';
        Hash = Coefs.Aggregate((P, X).GetHashCode(), (acc, a) => (acc, a).GetHashCode());
    }

    public FpPolynom(int p, int[] coefs)
    {
        if (!IntExt.Primes10000.Contains(p))
            throw new GroupException(GroupExceptionType.GroupDef);

        X = 'x';
        P = p;
        Coefs = coefs;
        Hash = Coefs.Aggregate((P, X).GetHashCode(), (acc, a) => (acc, a).GetHashCode());
    }

    public FpPolynom(int p, char x, int[] coefs)
    {
        if (!char.IsLetter(x))
            throw new GroupException(GroupExceptionType.GroupDef);

        X = x;
        P = p;
        Coefs = coefs;
        Hash = Coefs.Aggregate((P, X).GetHashCode(), (acc, a) => (acc, a).GetHashCode());
    }

    public static FpPolynom Fpx(int p) => new(p, 'x', new[] { 0, 1 });
    public static FpPolynom One(int p) => new(p, 'x', new[] { 1 });
    public static FpPolynom Zero(int p) => new(p, 'x', new[] { 0 });
    public static FpPolynom Fpx(int p, char x) => new(p, x, new[] { 0, 1 });
    public static FpPolynom One(int p, char x) => new(p, x, new[] { 1 });
    public static FpPolynom Zero(int p, char x) => new(p, x, new[] { 0 });

    public int CoefMax => Coefs.Last();
    public int Degree => Coefs.Length - 1;

    public FpPolynom Shift(int n) => new(P, X, new int[n].Concat(Coefs).ToArray());

    public bool Equals(FpPolynom other) => Hash == other.Hash;

    public int CompareTo(FpPolynom other) => X == other.X
        ? Coefs.SequenceCompareTo(other.Coefs)
        : throw new GroupException(GroupExceptionType.BaseGroup);

    public int Hash { get; }
    public override int GetHashCode() => Hash;

    public FpPolynom Pow(int n)
    {
        if (n <= 0)
            throw new GroupException(GroupExceptionType.GroupDef);

        var acc = new FpPolynom(P, X, new[] { 1 });
        for (int i = 0; i < n; i++)
            acc *= this;

        return acc;
    }

    public (FpPolynom quo, FpPolynom rem) DivideBy(FpPolynom b)
    {
        if (X != b.X)
            throw new GroupException(GroupExceptionType.BaseGroup);

        var qr = PolynomExt.DivPoly(P, Coefs, b.Coefs);
        var quo = new FpPolynom(P, X, qr.quo);
        var rem = new FpPolynom(P, X, qr.rem);
        return (quo, rem);
    }

    public override string ToString()
    {
        var str = "";
        for (int i = 0; i < Coefs.Length; i++)
            str = PolynomExt.PolyStr(X, str, i, Coefs[i]);

        return str == "" ? "0" : str;
    }

    public static FpPolynom operator *(FpPolynom a, FpPolynom b)
    {
        if (a.X != b.X)
            throw new GroupException(GroupExceptionType.BaseGroup);
        return new FpPolynom(a.P, a.X, PolynomExt.MulPoly(b.P, a.Coefs, b.Coefs));
    }

    public static FpPolynom operator *(int k, FpPolynom b)
    {
        return new FpPolynom(b.P, b.X, PolynomExt.MulKPoly(b.P, k, b.Coefs));
    }

    public static FpPolynom operator *(FpPolynom b, int k)
    {
        return new FpPolynom(b.P, b.X, PolynomExt.MulKPoly(b.P, k, b.Coefs));
    }

    public static FpPolynom operator /(FpPolynom b, int k)
    {
        var ki = IntExt.InvModP(k, b.P);
        return b * ki;
    }

    public static FpPolynom operator /(FpPolynom a, FpPolynom b) => a.DivideBy(b).quo;

    public static FpPolynom operator -(FpPolynom a) => -1 * a;

    public static FpPolynom operator +(FpPolynom a, FpPolynom b)
    {
        if (a.X != b.X)
            throw new GroupException(GroupExceptionType.BaseGroup);

        return new FpPolynom(a.P, a.X, PolynomExt.AddPoly(a.P, a.Coefs, b.Coefs));
    }

    public static FpPolynom operator +(FpPolynom a, int b)
    {
        return new FpPolynom(a.P, a.X, PolynomExt.AddPoly(a.P, a.Coefs, new[] { b }));
    }

    public static FpPolynom operator +(int a, FpPolynom b)
    {
        return new FpPolynom(b.P, b.X, PolynomExt.AddPoly(b.P, new[] { a }, b.Coefs));
    }

    public static FpPolynom operator -(FpPolynom a, int k) => a + (-k);

    public static FpPolynom operator -(int k, FpPolynom a) => k + (-a);

    public static FpPolynom operator -(FpPolynom a, FpPolynom b) => a + ((-1) * b);

    public static FpPolynom Gcd(FpPolynom a, FpPolynom b)
    {
        if (b.Degree == 0 && b.CoefMax == 0)
            return a;

        return Gcd(b, a.DivideBy(b).rem);
    }

    public static (FpPolynom x, FpPolynom y) Bezout(FpPolynom a, FpPolynom b)
    {
        if (b.Degree == 0 && b.CoefMax == 0)
            return (new(a.P, a.X, new[] { 1 }), b);


        var (quo, rem) = a.DivideBy(b);
        var (x0, y0) = Bezout(b, rem);
        return (y0, x0 - y0 * quo);
    }
}
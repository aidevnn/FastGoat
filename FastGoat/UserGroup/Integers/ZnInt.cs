using FastGoat.Commons;
using FastGoat.Structures;

namespace FastGoat.UserGroup.Integers;

public readonly struct ZnInt : IElt<ZnInt>, IRingElt<ZnInt>, IFieldElt<ZnInt>
{
    public int P { get; }
    public int K { get; }

    public ZnInt(int p, int k)
    {
        P = p;
        K = P == 0 ? k : IntExt.AmodP(k, P);
        Hash = (K, P).GetHashCode();
    }

    public bool Equals(ZnInt other)
    {
        return other.Hash == Hash;
    }

    public int CompareTo(ZnInt other)
    {
        if (P != other.P)
            throw new GroupException(GroupExceptionType.BaseGroup);

        return K.CompareTo(other.K);
    }

    public int Hash { get; }

    public override int GetHashCode()
    {
        return Hash;
    }

    public override string ToString()
    {
        var digits = $"{P}".Length;
        var fmt = $"{{0,{digits}}}";
        return string.Format(fmt, K);
    }

    public bool IsZero() => K == 0;
    public ZnInt Zero => new(P, 0);
    public ZnInt One => new(P, 1);

    public ZnInt Add(ZnInt e) => new(P, K + e.K);

    public ZnInt Sub(ZnInt e) => new(P, K - e.K);

    public ZnInt Opp() => new(P, -K);

    public ZnInt Mul(ZnInt e) => new(P, K * e.K);
    public ZnInt Mul(int k) => new(P, K * k);

    public ZnInt Pow(int k)
    {
        if (k == 0)
            return One;

        if (k < 0)
            return Inv().Pow(-k);

        var r = P != 0 ? IntExt.PowMod(K, k, P) : (int)Math.Pow(K, k);
        return new(P, r);
    }

    public (ZnInt quo, ZnInt rem) Div(ZnInt e)
    {
        var ek = P == 0 ? e.K : e.K % P;
        if (ek == 0)
            throw new DivideByZeroException();

        var q = P == 0 ? K / ek : K * IntExt.InvModP(ek, P);
        var r = K - q * ek;
        return (new(P, q), new(P, r));
    }

    public ZnInt Inv()
    {
        if (P == 0)
            throw new DivideByZeroException();

        var (x, y) = IntExt.Bezout(K, P);
        var gcd = K * x + P * y;
        return new(P, x / gcd);
    }

    public override bool Equals(object? obj)
    {
        return obj is ZnInt z && Equals(z);
    }

    public static bool operator ==(ZnInt a, ZnInt b) => a.Equals(b);
    public static bool operator !=(ZnInt a, ZnInt b) => !a.Equals(b);

    public static ZnInt operator +(ZnInt a, ZnInt b) => a.Add(b);
    public static ZnInt operator -(ZnInt a) => a.Opp();
    public static ZnInt operator -(ZnInt a, ZnInt b) => a + (-b);
    public static ZnInt operator *(ZnInt a, ZnInt b) => a.Mul(b);
    public static ZnInt operator /(ZnInt a, ZnInt b) => a.Mul(b.Inv());
    public static ZnInt operator *(ZnInt a, int k) => new(a.P, k * a.K);
    public static ZnInt operator *(int k, ZnInt a) => new(a.P, k * a.K);
    public static ZnInt operator /(ZnInt a, int k) => a.Mul(new ZnInt(a.P, k).Inv());
    
}
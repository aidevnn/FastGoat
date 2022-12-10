using FastGoat.Commons;
using FastGoat.Structures;

namespace FastGoat.UserGroup.Integers;

public readonly struct ZnInt : IElt<ZnInt>, IRingElt<ZnInt>, IFieldElt<ZnInt>
{
    public int P { get; }
    public int K { get; }

    public static ZnInt KZero(int p = 0) => new ZnInt(p, 0);

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

    public ZnInt LeadingCoeff => One;
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
        if (e.IsZero())
            throw new DivideByZeroException();

        if (IsZero())
            return (Zero, Zero);
        
        if (P == 0)
        {
            var (q, r) = Int32.DivRem(K, e.K);
            return (new(P, q), new(P, r));
        }
        else
        {
            var (x, y) = IntExt.Bezout(e.K, P);
            var gcd = e.K * x + P * y;
            if (x % gcd != 0)
            {
                var (q, r) = Int32.DivRem(K, e.K);
                return (new(P, q), new(P, r)); // not uniq result
            }
            else
            {
                var inv = (x / gcd) % P;
                var q = new ZnInt(P, inv * K);
                var r = new ZnInt(P, K - q.K * e.K);
                return (q, r); // r = 0 always
            }
        }
    }

    public ZnInt Inv()
    {
        if (P == 0)
        {
            if (K == 1)
                return new(0, 1);
            if (K == -1)
                return new(0, -1);
            
            throw new DivideByZeroException();
        }

        var (x, y) = IntExt.Bezout(K, P);
        var gcd = K * x + P * y;
        if (x % gcd != 0)
            throw new DivideByZeroException();
        
        return new(P, x / gcd);
    }

    public override bool Equals(object? obj)
    {
        return obj is ZnInt z && Equals(z);
    }

    public static bool operator ==(ZnInt a, ZnInt b) => a.Equals(b);
    public static bool operator !=(ZnInt a, ZnInt b) => !a.Equals(b);

    public static ZnInt operator +(ZnInt a, ZnInt b) => a.Add(b);
    public static ZnInt operator +(int a, ZnInt b) => b.Add(b.One.Mul(a));
    public static ZnInt operator +(ZnInt a, int b) => a.Add(a.One.Mul(b));
    public static ZnInt operator -(ZnInt a) => a.Opp();
    public static ZnInt operator -(ZnInt a, ZnInt b) => a + (-b);
    public static ZnInt operator -(int a, ZnInt b) => a + (-b);
    public static ZnInt operator -(ZnInt a, int b) => a + (-b);
    public static ZnInt operator *(ZnInt a, ZnInt b) => a.Mul(b);
    public static ZnInt operator *(ZnInt a, int b) => a.Mul(b);
    public static ZnInt operator *(int a, ZnInt b) => b.Mul(a);
    public static ZnInt operator /(ZnInt a, ZnInt b) => a.Div(b).quo;
    public static ZnInt operator /(ZnInt a, int b) => a.Div(a.One.Mul(b)).quo;
    public static ZnInt operator /(int a, ZnInt b) => b.Inv().Mul(a);
}
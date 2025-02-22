using FastGoat.Commons;
using FastGoat.Structures;

namespace FastGoat.UserGroup.Integers;

public enum ZnDisplay
{
    Unsigned,
    Signed
}

public readonly struct ZnInt : IElt<ZnInt>, IRingElt<ZnInt>, IFieldElt<ZnInt>
{
    public static ZnDisplay Display = ZnDisplay.Unsigned;
    public int Mod { get; }
    public int K { get; }
    public int P => Mod;

    public static double Abs(ZnInt z) => throw new(); // z.P == 0 ? double.Abs(z.K) : z.K;
    public static bool IsValuedField => false;
    public static ZnInt ZnZero(int m = 0) => new(m, 0);

    public static ZnInt ZpZero(int p = 2) =>
        IntExt.Primes10000.Contains(p) ? new(p, 0) : throw new GroupException(GroupExceptionType.GroupDef);

    public ZnInt(int mod, int k)
    {
        Mod = mod;
        K = Mod == 0 ? k : IntExt.AmodP(k, Mod);
        Hash = (K, Mod).GetHashCode();
    }

    public bool Equals(ZnInt other)
    {
        return other.Hash == Hash;
    }

    public int CompareTo(ZnInt other)
    {
        if (Mod != other.Mod)
            throw new GroupException(GroupExceptionType.BaseGroup);

        return K.CompareTo(other.K);
    }

    public int Hash { get; }

    public override int GetHashCode()
    {
        return Hash;
    }

    public int Signed => 2 * K > Mod ? K - Mod : K;

    public override string ToString()
    {
        var digits = $"{Mod}".Length;
        var fmt = $"{{0,{digits}}}";
        var k0 = Display == ZnDisplay.Unsigned ? K : Signed;
        return string.Format(fmt, k0);
    }

    public bool IsZero() => K == 0;
    public ZnInt Zero => new(Mod, 0);
    public ZnInt One => new(Mod, 1);

    public ZnInt Add(ZnInt e) => new(Mod, K + e.K);

    public ZnInt Sub(ZnInt e) => new(Mod, K - e.K);

    public ZnInt Opp() => new(Mod, -K);

    public ZnInt Mul(ZnInt e) => new(Mod, K * e.K);
    public ZnInt Mul(int k) => new(Mod, K * k);

    public ZnInt Pow(int k)
    {
        if (k == 0)
            return One;

        if (k < 0)
            return Inv().Pow(-k);

        var r = Mod != 0 ? IntExt.PowMod(K, k, Mod) : (int)Math.Pow(K, k);
        return new(Mod, r);
    }

    public (ZnInt quo, ZnInt rem) Div(ZnInt e)
    {
        if (e.IsZero())
            throw new DivideByZeroException();

        if (IsZero())
            return (Zero, Zero);

        if (Mod == 0)
        {
            var (q, r) = Int32.DivRem(K, e.K);
            return (new(Mod, q), new(Mod, r));
        }
        else
        {
            var (x, y) = IntExt.Bezout(e.K, Mod);
            var gcd = e.K * x + Mod * y;
            if (x % gcd != 0)
            {
                var (q, r) = Int32.DivRem(K, e.K);
                return (new(Mod, q), new(Mod, r)); // not uniq result
            }
            else
            {
                var inv = (x / gcd) % Mod;
                var q = new ZnInt(Mod, inv * K);
                var r = new ZnInt(Mod, K - q.K * e.K);
                return (q, r); // r = 0 always
            }
        }
    }

    public ZnInt Inv()
    {
        if (Mod == 0)
        {
            if (K == 1)
                return new(0, 1);
            if (K == -1)
                return new(0, -1);

            throw new DivideByZeroException();
        }

        var (x, y) = IntExt.Bezout(K, Mod);
        var gcd = K * x + Mod * y;
        if (x % gcd != 0)
            throw new DivideByZeroException();

        return new(Mod, x / gcd);
    }

    public bool Invertible() => IntExt.Gcd(K, P) == 1;

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
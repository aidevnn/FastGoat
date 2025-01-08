using FastGoat.Commons;
using FastGoat.Structures;

namespace FastGoat.UserGroup.Integers;

public readonly struct ZnInt64 : IElt<ZnInt64>, IRingElt<ZnInt64>, IFieldElt<ZnInt64>
{
    public static ZnDisplay Display = ZnDisplay.Unsigned;
    public long Mod { get; }
    public long K { get; }
    public int P => (int)Mod; // TODO: field characteristic

    public static double Abs(ZnInt64 z) => throw new(); // z.P == 0 ? double.Abs(z.K) : z.K;
    public static bool IsValuedField => false;
    public static ZnInt64 ZnZero(long m = 0) => new(m, 0);

    public static ZnInt64 ZpZero(int p = 2) =>
        IntExt.Primes10000.Contains(p) ? new(p, 0) : throw new GroupException(GroupExceptionType.GroupDef);

    public ZnInt64(long mod, long k)
    {
        Mod = mod;
        K = Mod == 0 ? k : IntExt.AmodPlong(k, Mod);
        Hash = (K, Mod).GetHashCode();
    }

    public bool Equals(ZnInt64 other)
    {
        return other.Hash == Hash;
    }

    public int CompareTo(ZnInt64 other)
    {
        if (Mod != other.Mod)
            throw new GroupException(GroupExceptionType.BaseGroup);

        return K.CompareTo(other.K);
    }

    public int Hash { get; }

    public override int GetHashCode() => Hash;

    public long Signed => 2 * K > Mod ? K - Mod : K;

    public override string ToString()
    {
        var digits = $"{Mod}".Length;
        var fmt = $"{{0,{digits}}}";
        var k0 = Display == ZnDisplay.Unsigned ? K : Signed;
        return string.Format(fmt, k0);
    }

    public bool IsZero() => K == 0;
    public ZnInt64 Zero => new(Mod, 0);
    public ZnInt64 One => new(Mod, 1);

    public ZnInt64 Add(ZnInt64 e) => new(Mod, K + e.K);

    public ZnInt64 Sub(ZnInt64 e) => new(Mod, K - e.K);

    public ZnInt64 Opp() => new(Mod, -K);

    public ZnInt64 Mul(ZnInt64 e) => new(Mod, K * e.K);
    public ZnInt64 Mul(int k) => new(Mod, K * k);

    public ZnInt64 Pow(int k)
    {
        if (k == 0)
            return One;

        if (k < 0)
            return Inv().Pow(-k);

        var r = Mod != 0 ? IntExt.PowModLong(K, k, Mod) : (long)Math.Pow(K, k);
        return new(Mod, r);
    }

    public (ZnInt64 quo, ZnInt64 rem) Div(ZnInt64 e)
    {
        if (e.IsZero())
            throw new DivideByZeroException();

        if (IsZero())
            return (Zero, Zero);

        if (Mod == 0)
        {
            var (q, r) = Int64.DivRem(K, e.K);
            return (new(Mod, q), new(Mod, r));
        }
        else
        {
            var (x, y) = IntExt.BezoutLong(e.K, Mod);
            var gcd = e.K * x + Mod * y;
            if (x % gcd != 0)
            {
                var (q, r) = Int64.DivRem(K, e.K);
                return (new(Mod, q), new(Mod, r)); // not uniq result
            }
            else
            {
                var inv = (x / gcd) % Mod;
                var q = new ZnInt64(Mod, inv * K);
                var r = new ZnInt64(Mod, K - q.K * e.K);
                return (q, r); // r = 0 always
            }
        }
    }

    public ZnInt64 Inv()
    {
        if (Mod == 0)
        {
            if (K == 1)
                return new(0, 1);
            if (K == -1)
                return new(0, -1);

            throw new DivideByZeroException();
        }

        var (x, y) = IntExt.BezoutLong(K, Mod);
        var gcd = K * x + Mod * y;
        if (x % gcd != 0)
            throw new DivideByZeroException();

        return new(Mod, x / gcd);
    }

    public bool Invertible() => IntExt.GcdLong(K, Mod) == 1;

    public override bool Equals(object? obj)
    {
        return obj is ZnInt64 z && Equals(z);
    }

    public static bool operator ==(ZnInt64 a, ZnInt64 b) => a.Equals(b);
    public static bool operator !=(ZnInt64 a, ZnInt64 b) => !a.Equals(b);

    public static ZnInt64 operator +(ZnInt64 a, ZnInt64 b) => a.Add(b);
    public static ZnInt64 operator +(int a, ZnInt64 b) => b + a;
    public static ZnInt64 operator +(ZnInt64 a, int b) => new(a.Mod, a.K + b);
    public static ZnInt64 operator +(ZnInt64 a, long b) => new(a.Mod, a.K + b);
    public static ZnInt64 operator +(long a, ZnInt64 b) => b + a;
    public static ZnInt64 operator -(ZnInt64 a) => a.Opp();
    public static ZnInt64 operator -(ZnInt64 a, ZnInt64 b) => a + (-b);
    public static ZnInt64 operator -(int a, ZnInt64 b) => a + (-b);
    public static ZnInt64 operator -(ZnInt64 a, int b) => a + (-b);
    public static ZnInt64 operator -(ZnInt64 a, long b) => a + (-b);
    public static ZnInt64 operator -(long a, ZnInt64 b) => a + (-b);
    public static ZnInt64 operator *(ZnInt64 a, ZnInt64 b) => a.Mul(b);
    public static ZnInt64 operator *(ZnInt64 a, int b) => a.Mul(b);
    public static ZnInt64 operator *(int a, ZnInt64 b) => b.Mul(a);
    public static ZnInt64 operator *(ZnInt64 a, long b) => new(a.Mod, a.K * b);
    public static ZnInt64 operator *(long a, ZnInt64 b) => b * a;
    public static ZnInt64 operator /(ZnInt64 a, ZnInt64 b) => a.Div(b).quo;
    public static ZnInt64 operator /(ZnInt64 a, int b) => a.Div(a.One.Mul(b)).quo;
    public static ZnInt64 operator /(int a, ZnInt64 b) => b.Inv() * a;
    public static ZnInt64 operator /(ZnInt64 a, long b) => a.Div(new(a.Mod, b)).quo;
    public static ZnInt64 operator /(long a, ZnInt64 b) => b.Inv() * a;
}
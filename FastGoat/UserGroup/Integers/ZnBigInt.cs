using System.Numerics;
using FastGoat.Commons;
using FastGoat.Structures;

namespace FastGoat.UserGroup.Integers;

public readonly struct ZnBigInt : IElt<ZnBigInt>, IRingElt<ZnBigInt>, IFieldElt<ZnBigInt>
{
    public static ZnDisplay Display = ZnDisplay.Unsigned;
    public BigInteger Mod { get; }
    public BigInteger K { get; }
    public int P => 0; // TODO: field characteristic

    public static double Abs(ZnBigInt z) => throw new(); // z.P == 0 ? double.Abs(z.K) : z.K;
    public static bool IsValuedField => false;
    public static ZnBigInt ZnZero(BigInteger m = default) => new(m, 0);

    public static ZnBigInt ZpZero(int p = 2) =>
        IntExt.Primes10000.Contains(p) ? new(p, 0) : throw new GroupException(GroupExceptionType.GroupDef);

    public ZnBigInt(BigInteger mod, BigInteger k)
    {
        Mod = mod;
        K = Mod == 0 ? k : IntExt.AmodPbigint(k % Mod, Mod);
        K = K * 2 > Mod ? K - Mod : K;
        Hash = (K, Mod).GetHashCode();
    }

    public bool Equals(ZnBigInt other) => (Mod, K).Equals((other.Mod, other.K));

    public int CompareTo(ZnBigInt other)
    {
        if (Mod != other.Mod)
            throw new GroupException(GroupExceptionType.BaseGroup);

        return K.CompareTo(other.K);
    }

    public int Hash { get; }

    public override int GetHashCode() => Hash;

    public override string ToString()
    {
        var digits = $"{Mod}".Length + 1;
        var fmt = $"{{0,{digits}}}";
        return string.Format(fmt, K);
    }

    public BigInteger Unsigned => K < 0 ? Mod + K : K;

    public bool IsZero() => K == 0;
    public ZnBigInt Zero => new(Mod, 0);
    public ZnBigInt One => new(Mod, 1);
    public ZnBigInt Rng => new(Mod, DistributionExt.Dice(0, Mod - 1));

    public ZnBigInt Add(ZnBigInt e) => new(Mod, K + e.K);

    public ZnBigInt Sub(ZnBigInt e) => new(Mod, K - e.K);

    public ZnBigInt Opp() => new(Mod, -K);

    public ZnBigInt Mul(ZnBigInt e) => new(Mod, K * e.K);
    public ZnBigInt Mul(int k) => new(Mod, K * k);

    public ZnBigInt Pow(int k)
    {
        if (k == 0)
            return One;

        if (k < 0)
            return Inv().Pow(-k);

        var r = Mod != 0 ? BigInteger.ModPow(K, k, Mod) : BigInteger.Pow(K, k);
        return new(Mod, r);
    }

    public (ZnBigInt quo, ZnBigInt rem) Div(ZnBigInt e)
    {
        if (e.IsZero())
            throw new DivideByZeroException();

        if (IsZero())
            return (Zero, Zero);

        if (Mod == 0)
        {
            var (q, r) = BigInteger.DivRem(K, e.K);
            return (new(Mod, q), new(Mod, r));
        }
        else
        {
            var (x, _, gcd) = IntExt.BezoutBigInt(e.K, Mod);
            if (x % gcd != 0)
            {
                var (q, r) = BigInteger.DivRem(K, e.K);
                return (new(Mod, q), new(Mod, r)); // not uniq result
            }
            else
            {
                var inv = (x / gcd) % Mod;
                var q = new ZnBigInt(Mod, inv * K);
                var r = new ZnBigInt(Mod, K - q.K * e.K);
                return (q, r); // r = 0 always
            }
        }
    }

    public ZnBigInt Inv()
    {
        if (Mod == 0)
        {
            if (K == 1)
                return new(0, 1);
            if (K == -1)
                return new(0, -1);

            throw new DivideByZeroException();
        }

        var (x, _, gcd) = IntExt.BezoutBigInt(K, Mod);
        if (x % gcd != 0)
            throw new DivideByZeroException();

        return new(Mod, x / gcd);
    }

    public bool Invertible() => IntExt.GcdBigInt(K, Mod) == 1;

    public override bool Equals(object? obj)
    {
        return obj is ZnBigInt z && Equals(z);
    }

    public static bool operator ==(ZnBigInt a, ZnBigInt b) => a.Equals(b);
    public static bool operator !=(ZnBigInt a, ZnBigInt b) => !a.Equals(b);

    public static ZnBigInt operator +(ZnBigInt a, ZnBigInt b) => a.Add(b);
    public static ZnBigInt operator +(int a, ZnBigInt b) => b + a;
    public static ZnBigInt operator +(ZnBigInt a, int b) => new(a.Mod, a.K + b);
    public static ZnBigInt operator +(ZnBigInt a, BigInteger b) => new(a.Mod, a.K + b);
    public static ZnBigInt operator +(BigInteger a, ZnBigInt b) => b + a;
    public static ZnBigInt operator -(ZnBigInt a) => a.Opp();
    public static ZnBigInt operator -(ZnBigInt a, ZnBigInt b) => a.Sub(b);
    public static ZnBigInt operator -(int a, ZnBigInt b) => a + (-b);
    public static ZnBigInt operator -(ZnBigInt a, int b) => a + (-b);
    public static ZnBigInt operator -(ZnBigInt a, BigInteger b) => a + (-b);
    public static ZnBigInt operator -(BigInteger a, ZnBigInt b) => a + (-b);
    public static ZnBigInt operator *(ZnBigInt a, ZnBigInt b) => a.Mul(b);
    public static ZnBigInt operator *(ZnBigInt a, int b) => a.Mul(b);
    public static ZnBigInt operator *(int a, ZnBigInt b) => b.Mul(a);
    public static ZnBigInt operator *(ZnBigInt a, BigInteger b) => new(a.Mod, a.K * b);
    public static ZnBigInt operator *(BigInteger a, ZnBigInt b) => b * a;
    public static ZnBigInt operator /(ZnBigInt a, ZnBigInt b) => a.Div(b).quo;
    public static ZnBigInt operator /(ZnBigInt a, int b) => a.Div(a.One.Mul(b)).quo;
    public static ZnBigInt operator /(int a, ZnBigInt b) => b.Inv() * a;
    public static ZnBigInt operator /(ZnBigInt a, BigInteger b) => a.Div(new(a.Mod, b)).quo;
    public static ZnBigInt operator /(BigInteger a, ZnBigInt b) => b.Inv() * a;
}
using System.Numerics;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.UserGroup.Integers;

namespace FastGoat.UserGroup.Padic;

public struct ZnBInt : IElt<ZnBInt>, IRingElt<ZnBInt>, IFieldElt<ZnBInt>
{
    public BigInteger Mod => Details.Mod;
    public BigInteger K { get; }
    public static ZnBInt KZero(int mod, int o = 1) => new(mod, 0, o);
    public Modulus Details { get; }

    public static double Abs(ZnBInt z)
    {
        if (z.IsZero())
            return 0;
        
        var a0 = z.K;
        for (int k = 0; k < z.Details.O; ++k)
        {
            a0 /= z.P;
            if (a0.IsZero)
                return k;
        }
        
        return z.Details.O;
    }
    public static bool IsValuedField => true;
    public int P => Details.P;

    public ZnBInt(int mod, BigInteger k, int o = 1)
    {
        Details = new(mod, o);
        var k0 = BigInteger.Remainder(k, Details.Mod);
        K = k0.Sign == -1 ? k0 + Details.Mod : k0;
        Hash = (K, Details).GetHashCode();
    }

    public ZnBInt(Modulus details, BigInteger k)
    {
        Details = details;
        var k0 = BigInteger.Remainder(k, details.Mod);
        K = k0.Sign == -1 ? k0 + details.Mod : k0;
        Hash = (K, Details).GetHashCode();
    }

    public bool Equals(ZnBInt other) => Hash == other.Hash;

    public int CompareTo(ZnBInt other) => K.CompareTo(other.K);

    public int Hash { get; }
    public bool IsZero() => K.IsZero;

    public ZnBInt Zero => new(Details, 0);
    public ZnBInt One => new(Details, 1);
    public ZnBInt Add(ZnBInt e) => new(Details, K + e.K);

    public ZnBInt Sub(ZnBInt e) => new(Details, K - e.K);

    public ZnBInt Opp() => new(Details, -K);

    public ZnBInt Mul(ZnBInt e) => new(Details, K * e.K);

    public (ZnBInt quo, ZnBInt rem) Div(ZnBInt e)
    {
        if (e.IsZero())
            throw new DivideByZeroException();

        var gcd0 = BigInteger.GreatestCommonDivisor(K, e.K);
        var e0 = K / gcd0;
        var e1 = e.K / gcd0;

        // var (x, _, gcd1) = IntExt.BezoutBigInt(e1, Mod);
        var (x, y) = BezoutBigInteger(e1, Mod);
        var gcd1 = e1 * x + Mod * y;
        if (x % gcd1 != 0)
        {
            var (q, r) = BigInteger.DivRem(K, e.K);
            return (new(Details, q), new(Details, r)); // not uniq result
        }
        else
        {
            var inv = (x / gcd1) % Mod;
            var q = new ZnBInt(Details, inv * e0);
            var r = new ZnBInt(Details, K - q.K * e.K);
            return (q, r); // r = 0 always
        }
    }

    public ZnBInt Mul(int k) => new(Details, K * k);

    public ZnBInt Pow(int k)
    {
        if (k == 0)
            return One;

        if (k < 0)
            return Inv().Pow(-k);

        var r = BigInteger.ModPow(K, k, Mod);
        return new(Details, r);
    }

    public ZnBInt Inv()
    {
        var (x, y) = BezoutBigInteger(K, Mod);
        var gcd = K * x + Mod * y;
        if (x % gcd != 0)
            throw new DivideByZeroException();

        return new(Details, x / gcd);
    }

    public bool Invertible() => !K.IsZero && BigInteger.GreatestCommonDivisor(K, P).IsOne;

    public int Sgn => K * 2 <= Mod ? 1 : -1;
    public BigInteger ToSignedBigInt => K * 2 <= Mod ? K : K - Mod;

    public string ToSignedBigIntString()
    {
        var digits = Mod.ToString().Length + 1;
        var fmt = $"{{0,{digits}}}";
        return string.Format(fmt, ToSignedBigInt);
    }

    public IEnumerable<int> PadicSequence()
    {
        var a0 = K;
        for (int k = 0; k < Details.O; ++k)
        {
            var (q, r) = BigInteger.DivRem(a0, P);
            var r0 = (int)r;
            yield return r0 < 0 ? P + r0 : r0;
            if (q.IsZero)
                break;

            a0 = q;
        }
    }

    public string PadicNumericString()
    {
        var pstr = Ring.Xi($"{P}", Details.O).ToString();
        var s = "";

        if (IsZero())
            return $"[0({pstr})]";

        var seq = PadicSequence().ToArray();
        for (int i = 0; i < seq.Length; i++)
        {
            var sep = P < 10 ? "" : "'";
            sep = i == 1 ? "." : sep;
            sep = i == 0 ? "" : sep;
            s = $"{seq[i]}" + sep + s;
        }

        s = s.Reverse().Glue();
        return $"[{s}({pstr})]";
    }

    public string LandauString()
    {
        var r1 = Rational.KOne();
        var xp = Ring.Polynomial($"{P}", r1);
        var pstr = $"O({xp.Pow(Details.O)})";
        var sgn = Sgn;
        var p = P;
        var seq = sgn == 1 ? PadicSequence().ToArray() : Opp().PadicSequence().ToArray();
        var sum = PadicSequence().Select((e, i) => (e, i))
            .Where(e => e.e != 0).Select(e => e.i == 0 ? $"{e.e}" : $"{e.e}·{xp.Pow(e.i)}").Append(pstr).ToArray();

        return sum.Glue(" + ");
    }

    public static IEnumerable<ZnBInt> Generate(int mod, int o = 1)
    {
        for (int k = 0; k < mod.Pow(o); ++k)
            yield return KZero(mod, o) * k;
    }

    public static (ZnBInt trunc, ZnBInt rem) Truncate(ZnBInt z, Modulus details)
    {
        if (details.O < 1)
            throw new($"tau:{details} sigma:{z.Details}");

        if (details.O > z.Details.O)
            return (KZero(z.P), new(details, z.K));

        var k0 = z.K % details.Mod;
        var mod1 = new Modulus(z.P, z.Details.O - details.O);
        if (k0 * 2 > details.Mod)
        {
            var r = k0 - details.Mod;
            var a = (z.K - r) / details.Mod;
            return (new ZnBInt(mod1, a), new ZnBInt(details, r));
        }
        else
        {
            var a = (z.K - k0) / details.Mod;
            return (new ZnBInt(mod1, a), new ZnBInt(details, k0));
        }
    }

    public static (ZnBInt trunc, ZnBInt rem) Truncate(ZnBInt z, int tau) => Truncate(z, new Modulus(z.P, tau));

    public static (BigInteger x, BigInteger y) BezoutBigInteger(BigInteger a, BigInteger b)
    {
        if (b == 0)
        {
            var x = a < 0 ? -1 : 1;
            return (x, 0);
        }

        var q = a / b;
        var (x0, y0) = BezoutBigInteger(b, a - b * q);
        return (y0, x0 - y0 * q);
    }

    public override int GetHashCode() => Hash;

    public override string ToString()
    {
        var digits = Mod.ToString().Length;
        var fmt = $"{{0,{digits}}}";
        return string.Format(fmt, K);
    }

    public static ZnBInt operator +(ZnBInt a, ZnBInt b) => a.Add(b);

    public static ZnBInt operator +(int a, ZnBInt b) => b.One.Mul(a).Add(b);

    public static ZnBInt operator +(ZnBInt a, int b) => a.Add(a.One.Mul(b));

    public static ZnBInt operator -(ZnBInt a) => a.Opp();

    public static ZnBInt operator -(ZnBInt a, ZnBInt b) => a.Sub(b);

    public static ZnBInt operator -(int a, ZnBInt b) => b.One.Mul(a).Sub(b);

    public static ZnBInt operator -(ZnBInt a, int b) => a.Sub(a.One.Mul(b));

    public static ZnBInt operator *(ZnBInt a, ZnBInt b) => a.Mul(b);

    public static ZnBInt operator *(int a, ZnBInt b) => b.Mul(a);

    public static ZnBInt operator *(ZnBInt a, int b) => a.Mul(b);

    public static ZnBInt operator *(ZnBInt a, BigInteger b) => new(a.Details, a.K * b);
    public static ZnBInt operator *(BigInteger a, ZnBInt b) => new(b.Details, b.K * a);

    public static ZnBInt operator /(ZnBInt a, ZnBInt b) => a.Div(b).quo;

    public static ZnBInt operator /(ZnBInt a, int b) => a.Div(a.One.Mul(b)).quo;

    public static ZnBInt operator /(int a, ZnBInt b) => b.Inv().Mul(a);
}
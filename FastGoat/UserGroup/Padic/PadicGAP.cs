using System.Numerics;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.UserGroup.Integers;

namespace FastGoat.UserGroup.Padic;

public readonly struct PadicGAP : IElt<PadicGAP>, IRingElt<PadicGAP>, IFieldElt<PadicGAP>
{
    public static PadicGAP KZero(int p, int o) => new(p, o);
    public Modulus Details { get; }
    public Valuation Val { get; }
    public BigInteger K { get; }
    public int P => Details.P;
    public int O => Details.O;
    public BigInteger Mod => Details.Mod;
    public (int, int, BigInteger, Valuation) POKV => (P, O, K, Val);
    public int Sign => K * 2 > Mod ? -1 : 1;

    public PadicGAP(int p, int o)
    {
        if (o < 1 || !IntExt.Primes10000.Contains(p))
            throw new ArgumentException($"p={p} must be prime and o={o} must be greater than 1.");

        Details = new(p, o);
        Val = new();
        K = BigInteger.Zero;
    }

    public PadicGAP(int p, int o, BigInteger k)
    {
        if (o < 1 || !IntExt.Primes10000.Contains(p))
            throw new ArgumentException($"p={p} must be prime and o={o} must be greater than 1.");

        Details = new(p, o);
        if (k.IsZero)
            Val = new();
        else
        {
            (Val, BigInteger k0) = Valuation.Of(p, k);
            var mod = BigInteger.Pow(p, o);
            var k1 = k0 % mod;
            K = k1 < 0 ? k1 + mod : k1;
        }
    }

    public PadicGAP(int p, int o, BigInteger k, int v)
    {
        if (o < 1 || !IntExt.Primes10000.Contains(p))
            throw new ArgumentException($"p={p} must be prime and o={o} must be greater than 1.");

        Details = new(p, o);
        if (k.IsZero)
            Val = new();
        else
        {
            var (v0, k0) = Valuation.Of(p, k);
            Val = v0 + new Valuation(v);
            var mod = BigInteger.Pow(p, o);
            var k1 = k0 % mod;
            K = k1 < 0 ? k1 + mod : k1;
        }
    }

    public PadicGAP(int p, int o, Rational r)
    {
        var qp = new PadicGAP(p, o, r.Num) / new PadicGAP(p, o, r.Denom);
        Details = qp.Details;
        K = qp.K;
        Val = qp.Val;
    }

    public PadicGAP(int p, int o, (BigInteger num, BigInteger denom) r)
    {
        var qp = new PadicGAP(p, o, new Rational(r.num, r.denom));
        Details = qp.Details;
        K = qp.K;
        Val = qp.Val;
    }

    private PadicGAP(Modulus details, BigInteger k, Valuation v)
    {
        Details = details;
        if (k.IsZero)
            Val = new();
        else
        {
            K = k;
            Val = v;
        }
    }

    public int Hash => POKV.GetHashCode();

    public bool Equals(PadicGAP other) => POKV.Equals(other.POKV);

    public int CompareTo(PadicGAP other)
    {
        var sign = Sign;
        var compSign = sign.CompareTo(other.Sign);
        if (compSign != 0)
            return compSign;

        var compV = Val.CompareTo(other.Val);
        if (compV != 0)
            return sign * compV;

        return sign * K.CompareTo(other.K);
    }

    public bool IsZero() => K.IsZero;

    public PadicGAP Zero => new PadicGAP(Details, 0, new());
    public PadicGAP One => new PadicGAP(Details, 1, new(0));
    public PadicGAP Ppow(int k) => new PadicGAP(P, O, 1, k);

    public PadicGAP Add(PadicGAP other)
    {
        if (!Details.Equals(other.Details))
            throw new ArgumentException($"Elements must have the same structure.");

        if (Val.IsInfinity)
            return other;

        if (other.Val.IsInfinity)
            return this;

        if (Val <= other.Val)
        {
            var o = other.Val - Val;
            if (o.V >= O)
                return new(Details, K, Val);

            var z1 = new ZnBInt(Details, K);
            var z2 = new ZnBInt(Details, other.K * BigInteger.Pow(P, o.V));
            var k = (z1 + z2).K;
            return new(P, O, k, Val.V);
        }
        else
        {
            var o = Val - other.Val;
            if (o.V >= O)
                return new(Details, other.K, other.Val);

            var z1 = new ZnBInt(Details, K * BigInteger.Pow(P, o.V));
            var z2 = new ZnBInt(Details, other.K);
            var k = (z1 + z2).K;
            return new(P, O, k, other.Val.V);
        }
    }

    public PadicGAP Opp() => new(P, O, new ZnBInt(Details, K).Opp().K, Val.V);
    public PadicGAP Sub(PadicGAP other) => Add(other.Opp());

    public PadicGAP Mul(PadicGAP other)
    {
        if (!Details.Equals(other.Details))
            throw new ArgumentException($"Elements must have the same structure.");

        if (K.IsZero || other.K.IsZero)
            return Zero;

        var v = Val + other.Val;
        var z1 = new ZnBInt(Details, K);
        var z2 = new ZnBInt(other.Details, other.K);
        var k = (z1 * z2).K;
        return new(P, O, k, v.V);
    }

    public PadicGAP Inv() => new(P, O, new ZnBInt(Details, K).Inv().K, -Val.V);

    public (PadicGAP quo, PadicGAP rem) Div(PadicGAP other) => (Mul(other.Inv()), Zero);

    public PadicGAP Mul(int k) => K.IsZero ? this : new(P, O, k * K, Val.V);

    public PadicGAP Pow(int k)
    {
        if (k == 0)
            return One;

        if (K.IsZero)
            return this;

        if (k < 0)
            return Inv().Pow(-k);

        return new(P, O, new ZnBInt(Details, K).Pow(k).K, k * Val.V);
    }

    public PadicGAP LeadingCoeff => One;
    public Rational Norm => Val.IsInfinity ? new(0) : new Rational(P).Pow(-Val.V);
    public PadicGAP Normalized => new(Details, K, new(0));

    public SortedList<int, int> Digits()
    {
        var table = new SortedList<int, int>() { [0] = 0 };
        if (K.IsZero)
            return table;

        var v = Val.V;
        for (int i = Int32.Min(0, v); i < Int32.Max(0, v + O); i++)
            table[i] = 0;

        var a0 = K;
        for (int i = v; i < v + O; ++i)
        {
            var (q, r) = BigInteger.DivRem(a0, P);
            var r0 = (int)r;
            table[i] = r0 < 0 ? P + r0 : r0;
            a0 = q;
        }

        return table;
    }

    public (int exp, BigInteger n) ToSignedBigInteger => K.IsZero ? (0, 0) : (Val.V, K * 2 > Mod ? K : K - Mod);

    public static PadicGAP Convert(int p, int v, params int[] coefs)
    {
        if (coefs.Length == 0 || !IntExt.Primes10000.Contains(p) || coefs[0] % p == 0)
            throw new ArgumentException();

        var o = coefs.Length;
        var k = coefs.Reverse().Aggregate(BigInteger.Zero, (acc, i) => acc * p + new ZnInt(p, i).K);
        return new(p, o, k, v);
    }

    public static PadicGAP operator +(PadicGAP a, PadicGAP b) => a.Add(b);

    public static PadicGAP operator +(int a, PadicGAP b) => b.One.Mul(a).Add(b);

    public static PadicGAP operator +(PadicGAP a, int b) => a.Add(a.One.Mul(b));

    public static PadicGAP operator -(PadicGAP a) => a.Opp();

    public static PadicGAP operator -(PadicGAP a, PadicGAP b) => a.Sub(b);

    public static PadicGAP operator -(int a, PadicGAP b) => b.One.Mul(a).Sub(b);

    public static PadicGAP operator -(PadicGAP a, int b) => a.Sub(a.One.Mul(b));

    public static PadicGAP operator *(PadicGAP a, PadicGAP b) => a.Mul(b);

    public static PadicGAP operator *(int a, PadicGAP b) => b.Mul(a);

    public static PadicGAP operator *(PadicGAP a, int b) => a.Mul(b);

    public static PadicGAP operator /(PadicGAP a, PadicGAP b) => a.Div(b).quo;

    public static PadicGAP operator /(PadicGAP a, int b) => a.Div(a.One.Mul(b)).quo;

    public static PadicGAP operator /(int a, PadicGAP b) => b.Inv().Mul(a);

    public static PadicGAP operator +(PadicGAP a, Rational b) => a.Add(new(a.P, a.O, b));
    public static PadicGAP operator +(Rational a, PadicGAP b) => new PadicGAP(b.P, b.O, a).Add(b);
    public static PadicGAP operator -(PadicGAP a, Rational b) => a.Sub(new(a.P, a.O, b));
    public static PadicGAP operator -(Rational a, PadicGAP b) => new PadicGAP(b.P, b.O, a).Sub(b);
    public static PadicGAP operator *(PadicGAP a, Rational b) => a.Mul(new PadicGAP(a.P, a.O, b));
    public static PadicGAP operator *(Rational a, PadicGAP b) => new PadicGAP(b.P, b.O, a).Mul(b);
    public static PadicGAP operator /(PadicGAP a, Rational b) => a.Div(new PadicGAP(a.P, a.O, b)).quo;
    public static PadicGAP operator /(Rational a, PadicGAP b) => new PadicGAP(b.P, b.O, a).Div(b).quo;

    public override int GetHashCode() => POKV.GetHashCode();

    public override string ToString()
    {
        var pstr = $"({P})~{O}";
        var s = "";

        if (K.IsZero)
            return $"[0{pstr}]";

        var Table = Digits();
        var Start = Table.Keys.Min();
        var emax = Table.Reverse().First(e => e.Value != 0).Key;
        for (int i = Start; i <= Int32.Max(0, emax); ++i)
        {
            var sep = P < 10 ? "" : "'";
            sep = i == 1 ? "." : sep;
            sep = i == Start ? "" : sep;
            s = s + sep + $"{Table[i]}";
        }

        s = s.Glue();
        return $"[{s}{pstr}]";
    }
}
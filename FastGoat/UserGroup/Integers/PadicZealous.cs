using System.Numerics;
using FastGoat.Commons;
using FastGoat.Structures;

namespace FastGoat.UserGroup.Integers;

public readonly struct PadicZealous : IElt<PadicZealous>, IRingElt<PadicZealous>, IFieldElt<PadicZealous>
{
    public static PadicZealous KZero(int p, int o) => new(p, o);
    public int N { get; }
    public int Val { get; }
    public BigInteger S { get; }
    public BigInteger Mod { get; }
    public int P { get; }
    public (int, int, int, BigInteger, BigInteger) PNVSM => (P, N, Val, S, Mod);

    public PadicZealous(int p, int n)
    {
        if (!IntExt.Primes10000.Contains(p))
            throw new ArgumentException($"p={p} must be prime.");

        P = p;
        N = n;
        Val = n;
        S = BigInteger.Zero;
        Mod = 1;
    }

    public PadicZealous(int p, int n, BigInteger k)
    {
        if (!IntExt.Primes10000.Contains(p))
            throw new ArgumentException($"p={p} must be prime.");

        var qp = new PadicZealous(p, n, k, 0);
        P = p;
        N = n;
        Val = qp.Val;
        S = qp.S;
        Mod = qp.Mod;
    }

    public PadicZealous(int p, int n, Rational r)
    {
        var qp = new PadicZealous(p, n, r.Num) / new PadicZealous(p, n, r.Denom);
        P = qp.P;
        N = qp.N;
        S = qp.S;
        Mod = qp.Mod;
        Val = qp.Val;
        if (S.IsZero)
        {
            Val = n;
            S = 0;
            Mod = 1;
        }
    }

    public PadicZealous(int p, int n, (BigInteger num, BigInteger denom) r)
    {
        var qp = new PadicZealous(p, n, new Rational(r.num, r.denom));
        P = qp.P;
        N = qp.N;
        S = qp.S;
        Mod = qp.Mod;
        Val = qp.Val;
        if (S.IsZero)
        {
            Val = n;
            S = 0;
            Mod = 1;
        }
    }

    private PadicZealous(int p, int n, int v, BigInteger s, BigInteger mod)
    {
        P = p;
        N = n;
        Val = v;
        S = s;
        Mod = mod;
        if (S.IsZero)
        {
            Val = n;
            S = 0;
            Mod = 1;
        }
    }

    private PadicZealous(int p, int n, BigInteger s, int v)
    {
        P = p;
        N = n;
        if (s.IsZero)
        {
            Val = n;
            S = 0;
            Mod = 1;
        }
        else
        {
            var (v0, k0) = PadicExt.GetValuation(p, s);
            if (v0 + v >= n)
            {
                Val = n;
                S = 0;
                Mod = 1;
            }
            else
            {
                Val = v0 + v;
                Mod = BigInteger.Pow(p, n - Val);
                var k1 = k0 % Mod;
                S = k1 < 0 ? k1 + Mod : k1;
                if (S.IsZero)
                {
                    Val = n;
                    S = 0;
                    Mod = 1;
                }
            }
        }
    }

    public BigInteger SignedK => S * 2 > Mod ? S - Mod : S;

    public int Sign => SignedK.Sign;

    public int Hash => PNVSM.GetHashCode();

    public bool Equals(PadicZealous other) => Sub(other).Norm.IsZero();

    public int CompareTo(PadicZealous other)
    {
        var sign = Sign;
        var compSign = sign.CompareTo(other.Sign);
        if (compSign != 0)
            return compSign;

        var compV = Val.CompareTo(other.Val);
        if (compV != 0)
            return sign * compV;

        return sign * S.CompareTo(other.S);
    }

    public bool IsZero() => S.IsZero;

    public PadicZealous Zero => new(P, N);
    public PadicZealous One => new(P, N, 1);
    public PadicZealous Ppow(int k) => new PadicZealous(P, N + k, 1, k);

    public PadicZealous Add(PadicZealous other)
    {
        if (P != other.P)
            throw new ArgumentException($"Elements must have the same prime p={P}.");

        var n0 = Int32.Min(N, other.N);

        if (Val <= other.Val)
        {
            var o = other.Val - Val;
            var k0 = S + other.S * BigInteger.Pow(P, o);
            return new(P, n0, k0, Val);
        }
        else
        {
            var o = Val - other.Val;
            var k0 = other.S + S * BigInteger.Pow(P, o);
            return new(P, n0, k0, other.Val);
        }
    }

    public PadicZealous Opp() => new(P, N, Val, Mod - S, Mod);
    public PadicZealous Sub(PadicZealous other) => Add(other.Opp());

    public PadicZealous Mul(PadicZealous other)
    {
        if (P != other.P)
            throw new ArgumentException($"Elements must have the same prime p={P}.");

        var n0 = Int32.Min(Val + other.N, N + other.Val);
        var v0 = Val + other.Val;
        var k0 = S * other.S;
        return new(P, n0, k0, v0);
    }

    public PadicZealous Inv()
    {
        var (x, y) = ZnBInt.BezoutBigInteger(S, Mod);
        var k0 = x / (x * S + y * Mod);
        var k1 = k0 < 0 ? k0 + Mod : k0;
        var n0 = N - 2 * Val;
        return new(P, n0, -Val, k1, Mod);
    }

    public (PadicZealous quo, PadicZealous rem) Div(PadicZealous other) => (Mul(other.Inv()), Zero);

    public PadicZealous Mul(int k) => Mul(new PadicZealous(P, N, k));

    public PadicZealous Pow(int k)
    {
        if (k == 0)
            return One;

        if (S.IsZero)
            return this;

        if (k < 0)
            return Inv().Pow(-k);

        var a0 = this;
        return Enumerable.Repeat(a0, k).Aggregate(One, (acc, e) => acc.Mul(e));
    }

    public PadicZealous LeadingCoeff => One;
    public Rational Norm => Val==N ? new(0) : new Rational(P).Pow(-Val);
    public PadicZealous Normalized => new(P, N, 0, S, Mod);

    public SortedList<int, int> Digits()
    {
        var table = new SortedList<int, int>() { [0] = 0 };
        if (S.IsZero)
            return table;

        var v = Val;
        for (int i = Int32.Min(0, v); i < Int32.Max(0, N); i++)
            table[i] = 0;

        var a0 = S;
        for (int i = v; i < N; ++i)
        {
            var (q, r) = BigInteger.DivRem(a0, P);
            var r0 = (int)r;
            table[i] = r0 < 0 ? P + r0 : r0;
            a0 = q;
        }

        return table;
    }

    public static PadicZealous Convert(int p, int v, params int[] coefs)
    {
        if (coefs.Length == 0 || !IntExt.Primes10000.Contains(p) || coefs[0] % p == 0)
            throw new ArgumentException();

        var o = coefs.Length;
        var k = coefs.Reverse().Aggregate(BigInteger.Zero, (acc, i) => acc * p + new ZnInt(p, i).K);
        return new(p, o, k, v);
    }

    public static PadicZealous operator +(PadicZealous a, PadicZealous b) => a.Add(b);

    public static PadicZealous operator +(int a, PadicZealous b) => b.One.Mul(a).Add(b);

    public static PadicZealous operator +(PadicZealous a, int b) => a.Add(a.One.Mul(b));

    public static PadicZealous operator -(PadicZealous a) => a.Opp();

    public static PadicZealous operator -(PadicZealous a, PadicZealous b) => a.Sub(b);

    public static PadicZealous operator -(int a, PadicZealous b) => b.One.Mul(a).Sub(b);

    public static PadicZealous operator -(PadicZealous a, int b) => a.Sub(a.One.Mul(b));

    public static PadicZealous operator *(PadicZealous a, PadicZealous b) => a.Mul(b);

    public static PadicZealous operator *(int a, PadicZealous b) => b.Mul(a);

    public static PadicZealous operator *(PadicZealous a, int b) => a.Mul(b);

    public static PadicZealous operator /(PadicZealous a, PadicZealous b) => a.Div(b).quo;

    public static PadicZealous operator /(PadicZealous a, int b) => a.Div(a.One.Mul(b)).quo;

    public static PadicZealous operator /(int a, PadicZealous b) => b.Inv().Mul(a);

    public static PadicZealous operator +(PadicZealous a, Rational b) => a.Add(new(a.P, a.N, b));
    public static PadicZealous operator +(Rational a, PadicZealous b) => new PadicZealous(b.P, b.N, a).Add(b);
    public static PadicZealous operator -(PadicZealous a, Rational b) => a.Sub(new(a.P, a.N, b));
    public static PadicZealous operator -(Rational a, PadicZealous b) => new PadicZealous(b.P, b.N, a).Sub(b);
    public static PadicZealous operator *(PadicZealous a, Rational b) => a.Mul(new PadicZealous(a.P, a.N, b));
    public static PadicZealous operator *(Rational a, PadicZealous b) => new PadicZealous(b.P, b.N, a).Mul(b);
    public static PadicZealous operator /(PadicZealous a, Rational b) => a.Div(new PadicZealous(a.P, a.N, b)).quo;
    public static PadicZealous operator /(Rational a, PadicZealous b) => new PadicZealous(b.P, b.N, a).Div(b).quo;

    public override int GetHashCode() => PNVSM.GetHashCode();

    public override string ToString()
    {
        var pstr = $"({P})~{N}";
        var s = "";

        if (S.IsZero)
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
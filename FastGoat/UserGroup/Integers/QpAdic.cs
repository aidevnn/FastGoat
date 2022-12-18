using System.Numerics;
using FastGoat.Commons;
using FastGoat.Structures;

namespace FastGoat.UserGroup.Integers;

public readonly struct QpAdic : IElt<QpAdic>, IRingElt<QpAdic>, IFieldElt<QpAdic>
{
    public static QpAdic KZero(int p, int o) => new(p, o);
    public Modulus Details { get; }
    public int V { get; } // TODO when x=0 val(x) = infinity
    public BigInteger K { get; }
    public int P => Details.P;
    public int O => Details.O;
    public BigInteger Mod => Details.Mod;
    public (int, int, BigInteger, int) POKV => (P, O, K, V);
    public int Sign => K * 2 > Mod ? -1 : 1;

    public QpAdic(int p, int o)
    {
        if (o < 1 || !IntExt.Primes10000.Contains(p))
            throw new ArgumentException($"p={p} must be prime and o={o} must be greater than 1.");

        Details = new(p, o);
        V = 0;
        K = BigInteger.Zero;
    }

    public QpAdic(int p, int o, BigInteger k)
    {
        if (o < 1 || !IntExt.Primes10000.Contains(p))
            throw new ArgumentException($"p={p} must be prime and o={o} must be greater than 1.");

        Details = new(p, o);
        var (v0, k0) = PadicExt.GetValuation(p, k);
        var mod = BigInteger.Pow(p, o);
        var k1 = k0 % mod;
        K = k1 < 0 ? k1 + mod : k1;
        V = v0;

        if (K.IsZero)
            V = 0;
    }

    public QpAdic(int p, int o, BigInteger k, int v)
    {
        if (o < 1 || !IntExt.Primes10000.Contains(p))
            throw new ArgumentException($"p={p} must be prime and o={o} must be greater than 1.");

        Details = new(p, o);
        var (v0, k0) = PadicExt.AddValuation(p, k, v);
        var mod = BigInteger.Pow(p, o);
        var k1 = k0 % mod;
        K = k1 < 0 ? k1 + mod : k1;
        V = v0;

        if (K.IsZero)
            V = 0;
    }

    public QpAdic(int p, int o, Rational r)
    {
        var qp = new QpAdic(p, o, r.Num) / new QpAdic(p, o, r.Denom);
        Details = qp.Details;
        K = qp.K;
        V = qp.V;
    }

    public QpAdic(int p, int o, (BigInteger num, BigInteger denom) r)
    {
        var qp = new QpAdic(p, o, new Rational(r.num, r.denom));
        Details = qp.Details;
        K = qp.K;
        V = qp.V;
    }
    
    private QpAdic(Modulus details, BigInteger k, int v)
    {
        V = v;
        K = k;
        Details = details;
        if (K.IsZero)
            V = 0;
    }

    public int Hash => POKV.GetHashCode();
    
    public bool Equals(QpAdic other) => POKV.Equals(other.POKV);

    public int CompareTo(QpAdic other)
    {
        var sign = Sign;
        var compSign = sign.CompareTo(other.Sign);
        if (compSign != 0)
            return compSign;

        var compV = V.CompareTo(other.V);
        if (compV != 0)
            return sign * compV;

        return sign * K.CompareTo(other.K);
    }

    public bool IsZero() => K.IsZero;

    public QpAdic Zero => new QpAdic(Details, 0, 0);
    public QpAdic One => new QpAdic(Details, 1, 0);

    public QpAdic Add(QpAdic other)
    {
        if (!Details.Equals(other.Details))
            throw new ArgumentException($"Elements must have the same structure.");

        if (V <= other.V)
        {
            var o = other.V - V;
            if (o >= O)
                return new(Details, K, V);

            var z1 = new ZnBInt(Details, K);
            var z2 = new ZnBInt(Details, other.K * BigInteger.Pow(P, o));
            var k = (z1 + z2).K;
            return new(Details, k, V);
        }
        else
        {
            var o = V - other.V;
            if (o >= O)
                return new(Details, other.K, other.V);

            var z1 = new ZnBInt(Details, K * BigInteger.Pow(P, o));
            var z2 = new ZnBInt(Details, other.K);
            var k = (z1 + z2).K;
            return new(Details, k, other.V);
        }
    }

    public QpAdic Opp() => new(Details, new ZnBInt(Details, K).Opp().K, V);
    public QpAdic Sub(QpAdic other) => Add(other.Opp());

    public QpAdic Mul(QpAdic other)
    {
        if (!Details.Equals(other.Details))
            throw new ArgumentException($"Elements must have the same structure.");

        var v = V + other.V;
        var z1 = new ZnBInt(Details, K);
        var z2 = new ZnBInt(other.Details, other.K);
        var k = (z1 * z2).K;
        return new(Details, k, v);
    }

    public QpAdic Inv() => new(Details, new ZnBInt(Details, K).Inv().K, -V);

    public (QpAdic quo, QpAdic rem) Div(QpAdic other) => (Mul(other.Inv()), Zero);

    public QpAdic Mul(int k) => new(P, O, k * K, V);

    public QpAdic Pow(int k)
    {
        if (k == 0)
            return One;

        if (k < 0)
            return Inv().Pow(-k);

        return new(Details, new ZnBInt(Details, K).Pow(k).K, k * V);
    }

    public QpAdic LeadingCoeff => One;
    public Rational Norm => new Rational(P).Pow(-V);
    public QpAdic Normalized => new(Details, K, 0);

    public SortedList<int, int> Digits()
    {
        var table = new SortedList<int, int>() { [0] = 0 };
        for (int i = Int32.Min(0, V); i < Int32.Max(0, V + O); i++)
            table[i] = 0;
        
        var a0 = K;
        for (int i = V; i < V + O; ++i)
        {
            var (q, r) = BigInteger.DivRem(a0, P);
            var r0 = (int)r;
            table[i] = r0 < 0 ? P + r0 : r0;
            a0 = q;
        }

        return table;
    }

    public (int exp, BigInteger n) ToSignedBigInteger => (V, K * 2 > Mod ? K : K - Mod);

    public static QpAdic operator +(QpAdic a, QpAdic b) => a.Add(b);

    public static QpAdic operator +(int a, QpAdic b) => b.One.Mul(a).Add(b);

    public static QpAdic operator +(QpAdic a, int b) => a.Add(a.One.Mul(b));

    public static QpAdic operator -(QpAdic a) => a.Opp();

    public static QpAdic operator -(QpAdic a, QpAdic b) => a.Sub(b);

    public static QpAdic operator -(int a, QpAdic b) => b.One.Mul(a).Sub(b);

    public static QpAdic operator -(QpAdic a, int b) => a.Sub(a.One.Mul(b));

    public static QpAdic operator *(QpAdic a, QpAdic b) => a.Mul(b);

    public static QpAdic operator *(int a, QpAdic b) => b.Mul(a);

    public static QpAdic operator *(QpAdic a, int b) => a.Mul(b);

    public static QpAdic operator /(QpAdic a, QpAdic b) => a.Div(b).quo;

    public static QpAdic operator /(QpAdic a, int b) => a.Div(a.One.Mul(b)).quo;

    public static QpAdic operator /(int a, QpAdic b) => b.Inv().Mul(a);
    
    public static QpAdic operator +(QpAdic a, Rational b) => a.Add(new(a.P, a.O, b));
    public static QpAdic operator +(Rational a, QpAdic b) => new QpAdic(b.P, b.O, a).Add(b);
    public static QpAdic operator -(QpAdic a, Rational b) => a.Sub(new(a.P, a.O, b));
    public static QpAdic operator -(Rational a, QpAdic b) => new QpAdic(b.P, b.O, a).Sub(b);
    public static QpAdic operator *(QpAdic a, Rational b) => a.Mul(new QpAdic(a.P, a.O, b));
    public static QpAdic operator *(Rational a, QpAdic b) => new QpAdic(b.P, b.O, a).Mul(b);
    public static QpAdic operator /(QpAdic a, Rational b) => a.Div(new QpAdic(a.P, a.O, b)).quo;
    public static QpAdic operator /(Rational a, QpAdic b) => new QpAdic(b.P, b.O, a).Div(b).quo;
    
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
            s = $"{Table[i]}" + sep + s;
        }

        s = s.Reverse().Glue();
        return $"[{s}{pstr}]";
    }
}

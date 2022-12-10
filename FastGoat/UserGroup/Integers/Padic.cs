using System.Numerics;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;

namespace FastGoat.UserGroup.Integers;

public struct Padic : IVsElt<ZnInt, Padic>, IElt<Padic>, IRingElt<Padic>, IFieldElt<Padic>
{
    public ZnInt[] Coefs { get; }
    public int O { get; }
    public int P { get; }

    public Padic(int p, int o)
    {
        O = o;
        P = p;
        Coefs = Enumerable.Repeat(ZnInt.KZero(p), O).ToArray();
        Hash = Coefs.Aggregate(ZnInt.KZero(P).Hash, (acc, a) => (acc, a.Hash).GetHashCode());
        Valuation = 1;
    }

    public Padic(ZnInt[] coefs)
    {
        O = coefs.Length;
        P = coefs[0].P;
        Coefs = coefs;
        Valuation = Coefs.Select((z, i) => (z, i)).FirstOrDefault(e => !e.Item1.IsZero(), (ZnInt.KZero(P), 0)).Item2;
        Hash = Coefs.Aggregate(ZnInt.KZero(P).Hash, (acc, a) => (acc, a.Hash).GetHashCode());
    }

    public bool Equals(Padic other) => P == other.P && Coefs.SequenceEqual(other.Coefs);

    public bool EqualsApproxOne()
    {
        var e0 = KOne;
        return Coefs[0].Equals(e0) && Coefs.Skip(1).SkipLast(1).All(e => e.IsZero());
    }

    public int CompareTo(Padic other) => Coefs.Reverse().SequenceCompareTo(other.Coefs.Reverse());

    public Padic Resize(int o)
    {
        var coefs = Enumerable.Repeat(ZnInt.KZero(P), o).ToArray();
        for (int i = 0; i < Int32.Min(o, O); i++)
            coefs[i] = Coefs[i];

        return new(coefs);
    }

    public ZnInt KZero => new(P, 0);
    public ZnInt KOne => new(P, 1);
    public int Hash { get; }
    public bool IsZero() => Coefs.All(e => e.IsZero());

    public Padic Zero => new(P, O);

    public Padic One
    {
        get
        {
            var coefs = Enumerable.Repeat(ZnInt.KZero(P), O).ToArray();
            coefs[0] = new ZnInt(P, 1);
            return new(coefs);
        }
    }

    public ZnInt this[int index] => index >= 0 && index < O ? Coefs[index] : KZero;

    public Padic Add(Padic e)
    {
        var coefs = Enumerable.Repeat(ZnInt.KZero(P), O).ToArray();
        var carry = 0;
        for (int i = 0; i < O; i++)
        {
            var k = Coefs[i].K + e.Coefs[i].K + carry;
            var z = coefs[i] = new ZnInt(P, k);
            carry = (k - z.K) / P;
        }

        return new(coefs);
    }

    public Padic Sub(Padic e)
    {
        var coefs = Enumerable.Repeat(ZnInt.KZero(P), O).ToArray();
        var carry = 0;
        for (int i = 0; i < O; i++)
        {
            var k0 = Coefs[i].K;
            var k1 = e.Coefs[i].K + carry;
            if (k0 >= k1)
            {
                var k = k0 - k1;
                var z = coefs[i] = new ZnInt(P, k);
                carry = (k - z.K) / P;
            }
            else
            {
                var k = P + k0 - k1;
                var z = coefs[i] = new ZnInt(P, k);
                carry = (P + k - z.K) / P;
            }
        }

        return new(coefs);
    }

    public Padic Opp() => Zero - this;

    public Padic Mul(Padic e)
    {
        var coefs = Enumerable.Repeat(ZnInt.KZero(P), O).ToArray();
        var carry = 0;
        for (int i = 0; i < O; i++)
        {
            var sum = carry;
            for (int k0 = 0; k0 <= i; k0++)
            {
                var k1 = i - k0;
                sum += Coefs[k0].K * e.Coefs[k1].K;
            }

            var e0 = coefs[i] = new ZnInt(P, sum);
            carry = (sum - e0.K) / P;
        }

        return new(coefs);
    }

    public Padic Mul(int k)
    {
        var coefs = Enumerable.Repeat(ZnInt.KZero(P), O).ToArray();
        var carry = 0;
        for (int i = 0; i < O; i++)
        {
            var k0 = Coefs[i].K * k + carry;
            var z = coefs[i] = new ZnInt(P, k0);
            carry = (k0 - z.K) / P;
        }

        return new(coefs);
    }

    public (Padic quo, Padic rem) Div(Padic e)
    {
        if (e.IsZero())
            throw new DivideByZeroException();
        
        var e0 = e.Trimed;
        var coefs = Enumerable.Repeat(ZnInt.KZero(P), O).ToArray();
        var rem = new Padic(Coefs);
        var b0 = e0[0].Inv();
        for (int i = 0; i < O && !rem.IsZero(); i++)
        {
            var a0 = rem[0] * b0;
            coefs[i] = a0;
            rem = (rem - a0 * e0).Shift();
        }

        var e1 = new Padic(coefs).Shift(e.Valuation);
        var r = Sub(e * e1);
        return (e1, r);
    }

    public Padic Pow(int k)
    {
        if (k == 0)
            return One;

        if (k < 0)
            throw new GroupException(GroupExceptionType.GroupDef);

        var pi = this;
        return Enumerable.Repeat(pi, k).Aggregate((a, b) => a.Mul(b));
    }

    public Padic Inv()
    {
        var qr = One.Div(this);
        if (qr.rem.IsZero())
            return qr.quo;

        throw new DivideByZeroException();
    }

    public Padic KMul(ZnInt k) => Mul(k.K);

    public int Valuation { get; }

    public Padic Trimed
    {
        get
        {
            var coefs = Enumerable.Repeat(ZnInt.KZero(P), O).ToArray();
            for (int i = Valuation, j = 0; i < O; i++, j++)
                coefs[j] = Coefs[i];

            return new(coefs);
        }
    }

    public Padic Shift(int s = 1)
    {
        if (s == 0)
            return new(Coefs);
        else if (s > 0)
        {
            var coefs = Enumerable.Repeat(ZnInt.KZero(P), O).ToArray();
            for (int i = s; i < O; i++)
                coefs[i - s] = Coefs[i];

            return new(coefs);
        }
        else
        {
            var coefs = Enumerable.Repeat(ZnInt.KZero(P), O).ToArray();
            for (int i = -s; i < O; i++)
                coefs[i] = Coefs[i + s];

            return new(coefs);
        }
    }

    public Padic PpowV
    {
        get
        {
            if (IsZero())
                return Zero;

            var v = Valuation;
            var z0 = KZero;
            var coefs = Coefs.Select((e, i) => i == v ? e.One : z0).ToArray();
            return new(coefs);
        }
    }

    public Padic LeadingCoeff
    {
        get
        {
            if (IsZero())
                return Zero;

            var v = Valuation;
            var z0 = KZero;
            var coefs = Coefs.Select((e, i) => i == v ? e : z0).ToArray();
            return new(coefs);
        }
    }

    public BigInteger ToBigInt
    {
        get
        {
            var acc = BigInteger.Zero;
            var ppow = BigInteger.One;
            foreach (var z in Coefs)
            {
                acc += z.K*ppow;
                ppow *= P;
            }

            return acc;
        }
    }

    public static IEnumerable<Padic> Generate(int p, int o)
    {
        var e0 = new Padic(p, o);
        var e1 = e0.One;
        do
        {
            yield return e0;
            e0 += e1;
        } while (!e0.IsZero());
    }

    public static Padic Convert(int p, int o, BigInteger a)
    {
        var coefs = Enumerable.Repeat(ZnInt.KZero(p), o).ToArray();
        var a0 = a;
        var i = 0;
        while (i < o && !a0.IsZero)
        {
            var r = (int)BigInteger.Remainder(a0, p);
            var z = coefs[i] = new ZnInt(p, r);
            a0 = (a0 - z.K) / p;
            ++i;
        }

        return new(coefs);
    }

    public static Padic operator +(Padic a, Padic b) => a.Add(b);

    public static Padic operator +(int a, Padic b) => b.Add(b.One.Mul(a));

    public static Padic operator +(Padic a, int b) => a.Add(a.One.Mul(b));

    public static Padic operator -(Padic a) => a.Opp();

    public static Padic operator -(Padic a, Padic b) => a.Sub(b);

    public static Padic operator -(int a, Padic b) => b.Sub(b.One.Mul(a));

    public static Padic operator -(Padic a, int b) => a.Sub(a.One.Mul(b));

    public static Padic operator *(Padic a, Padic b) => a.Mul(b);

    public static Padic operator *(int a, Padic b) => b.Mul(a);

    public static Padic operator *(Padic a, int b) => a.Mul(b);

    public static Padic operator /(Padic a, Padic b) => a.Div(b).quo;

    public static Padic operator /(Padic a, int b) => a / (a.One * b);

    public static Padic operator /(int a, Padic b) => (a * b.One) / b;

    public static Padic operator +(Padic a, ZnInt b) => a + b.K;

    public static Padic operator +(ZnInt a, Padic b) => a.K + b;

    public static Padic operator -(Padic a, ZnInt b) => a - b.K;

    public static Padic operator -(ZnInt a, Padic b) => a.K - b;

    public static Padic operator *(Padic a, ZnInt b) => a * b.K;

    public static Padic operator *(ZnInt a, Padic b) => a.K * b;

    public static Padic operator /(Padic a, ZnInt b) => a / b.K;

    public override int GetHashCode() => Hash;

    public override string ToString()
    {
        var xi = Ring.Polynomial('p', KZero);
        var pO = Ring.Xi('p', O);
        var fx = Coefs.Select((k, i) => k * xi.Pow(i)).Aggregate((a, b) => a + b).GetString(reverse: true);
        return $"{fx} + O({pO})";
    }
}
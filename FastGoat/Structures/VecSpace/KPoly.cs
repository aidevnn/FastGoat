using System.Numerics;
using FastGoat.Commons;

namespace FastGoat.Structures.VecSpace;

public readonly struct KPoly<K> : IVsElt<K, KPoly<K>>, IElt<KPoly<K>>, IRingElt<KPoly<K>>
    where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
{
    public K[] Coefs { get; }
    public int Degree { get; }
    public K KZero { get; }
    public K KOne => KZero.One;
    public int P { get; }
    public char x { get; }

    public KPoly(char x0)
    {
        KZero = new K().Zero;
        P = KZero.P;
        x = x0;
        Coefs = new[] { KZero, KOne };
        Hash = Coefs.Aggregate(0, (acc, a) => (acc, a.Hash).GetHashCode());
        Degree = Coefs.Length - 1;
    }

    public KPoly(char x0, K k0)
    {
        KZero = k0.Zero;
        P = k0.P;
        x = x0;
        Coefs = new[] { k0 };
        Hash = Coefs.Aggregate(0, (acc, a) => (acc, a.Hash).GetHashCode());
        Degree = Coefs.Length - 1;
    }

    public KPoly(char x0, K kZero, K[] coefs)
    {
        KZero = kZero.Zero;
        P = kZero.P;
        x = x0;

        Coefs = coefs.Length != 0 ? coefs : new[] { KZero };
        Hash = Coefs.Aggregate(0, (acc, a) => (acc, a.Hash).GetHashCode());
        Degree = Coefs.Length - 1;
    }

    public K this[int idx]
    {
        get
        {
            if (idx < 0 || idx > Degree)
                return KZero;

            return Coefs[idx];
        }
    }

    public bool Equals(KPoly<K> other) => Coefs.SequenceEqual(other.Coefs); // Avoid collisions

    public int CompareTo(KPoly<K> other)
    {
        var compDegree = Degree.CompareTo(other.Degree);
        if (compDegree != 0)
            return compDegree;

        for (int i = Degree; i >= 0; i--)
        {
            var s = (Degree - i) % 2 == 0 ? 1 : -1;
            var comp = s * this[i].CompareTo(other[i]);
            if (comp != 0)
                return comp;
        }

        return 0;
    }

    public int Hash { get; }

    public bool IsZero() => Degree == 0 && Coefs[0].IsZero();

    public KPoly<K> Zero => new(x, KZero, new[] { KZero });
    public KPoly<K> ZeroExtended(int degree) => new(x, KZero, Enumerable.Repeat(KZero, degree + 1).ToArray());
    public KPoly<K> One => new(x, KZero, new[] { KOne });
    public KPoly<K> X => new(x, KZero, new[] { KZero, KOne });

    public KPoly<K> Derivative => new(x, KZero, Coefs.Select((e, i) => e.Mul(i)).Skip(1).TrimSeq().ToArray());
    public K Substitute(K f) => Coefs.Select((k, i) => k * f.Pow(i)).Aggregate((a, b) => a + b);
    public KPoly<K> Substitute(KPoly<K> f) => Coefs.Select((k, i) => k * f.Pow(i)).Aggregate((a, b) => a + b);
    public EPoly<K> Substitute(EPoly<K> f) => Coefs.Select((k, i) => k * f.Pow(i)).Aggregate((a, b) => a + b);

    public KPoly<EPoly<K>> Substitute(KPoly<EPoly<K>> f)
    {
        var poly = new KPoly<EPoly<K>>(f.x, f.KZero, Coefs.Select(k => k * f.KOne).ToArray());
        return poly.Substitute(f);
    }

    public KPoly<FracPoly<K>> Substitute(KPoly<FracPoly<K>> f)
    {
        var poly = new KPoly<FracPoly<K>>(f.x, f.KZero, Coefs.Select(k => k * f.KOne).ToArray());
        return poly.Substitute(f);
    }

    public KPoly<K> Add(KPoly<K> e)
    {
        var maxDegree = Math.Max(Degree, e.Degree);
        var coefs = new K[maxDegree + 1];
        for (int i = 0; i <= maxDegree; i++)
            coefs[i] = this[i].Add(e[i]);

        return new(x, KZero, coefs.TrimSeq().ToArray());
    }

    public KPoly<K> Sub(KPoly<K> e)
    {
        var maxDegree = Math.Max(Degree, e.Degree);
        var coefs = new K[maxDegree + 1];
        for (int i = 0; i <= maxDegree; i++)
            coefs[i] = this[i].Sub(e[i]);

        return new(x, KZero, coefs.TrimSeq().ToArray());
    }

    public KPoly<K> Opp() => Zero.Sub(this);

    public KPoly<K> Mul(KPoly<K> e)
    {
        var deg = Degree + e.Degree;
        var coefs = Enumerable.Repeat(KZero, deg + 1).ToArray();
        for (int i = 0; i <= Degree; i++)
        {
            var ai = this[i];
            for (int j = 0; j <= e.Degree; j++)
                coefs[i + j] = coefs[i + j].Add(ai.Mul(e[j]));
        }

        return new(x, KZero, coefs.TrimSeq().ToArray());
    }

    public void InPlaceAdd(KPoly<K> a)
    {
        if (a.Degree > Degree)
            throw new GroupException(GroupExceptionType.GroupDef);

        for (int i = 0; i <= a.Degree; i++)
            Coefs[i] += a[i];
    }

    public void InPlaceAddProd(KPoly<K> a, KPoly<K> b)
    {
        if (a.Degree + b.Degree > Degree)
            throw new GroupException(GroupExceptionType.GroupDef);

        for (int i = 0; i <= a.Degree; i++)
        {
            var ai = a[i];
            for (int j = 0; j <= b.Degree; j++)
                Coefs[i + j] += ai * b[j];
        }
    }

    public KPoly<K> KMul(K k)
    {
        var coefs = Coefs.Select(e => e.Mul(k)).TrimSeq().ToArray();
        return new(x, KZero, coefs);
    }

    public KPoly<K> Mul(int k)
    {
        var k0 = KOne.Mul(k);
        return KMul(k0);
    }

    public KPoly<K> Monic => IsZero() ? this : KMul(Coefs.Last().Inv());

    public KPoly<K> Pow(int k)
    {
        if (k == 0)
            return One;

        if (k < 0)
            throw new GroupException(GroupExceptionType.GroupDef);

        var pi = this;
        return Enumerable.Repeat(pi, k).Aggregate((a, b) => a.Mul(b));
    }

    public (KPoly<K> quo, KPoly<K> rem) Div(KPoly<K> e)
    {
        if (e.IsZero())
            throw new DivideByZeroException();

        if (Degree < e.Degree)
            return (Zero, new(x, KZero, Coefs));

        var em = e.Coefs.Last();
        var quo = Enumerable.Repeat(KZero, Degree - e.Degree + 1).ToArray();
        var rem = Coefs.ToArray();
        for (int i = Degree; i >= e.Degree; i--)
        {
            var ai = rem[i];
            var qr = ai.Div(em);
            if (!qr.rem.IsZero())
                throw new GroupException(GroupExceptionType.GroupDef);

            quo[i - e.Degree] = qr.quo;
            for (int j = 0; j <= i; j++)
                rem[j] = rem[j].Sub(e[e.Degree - i + j].Mul(qr.quo));
        }

        return (new(x, KZero, quo.TrimSeq().ToArray()), new(x, KZero, rem.TrimSeq().ToArray()));
    }

    public override int GetHashCode() => Hash;

    public override string ToString()
    {
        var xi = Ring.Polynomial(KZero, $"{x}")[0];
        var fx = Coefs.Select((k, i) => k * xi.Pow(i)).Aggregate((a, b) => a + b);
        return $"{fx}";
    }

    public EPoly<K> ToEPoly(KPoly<K> f) => new(f, this);

    public static KPoly<K> operator +(KPoly<K> a, KPoly<K> b) => a.Add(b);
    public static KPoly<K> operator +(int a, KPoly<K> b) => b.Add(b.One.Mul(a));
    public static KPoly<K> operator +(KPoly<K> a, int b) => a.Add(a.One.Mul(b));
    public static KPoly<K> operator -(KPoly<K> a) => a.Opp();
    public static KPoly<K> operator -(KPoly<K> a, KPoly<K> b) => a + (-b);
    public static KPoly<K> operator -(int a, KPoly<K> b) => a + (-b);
    public static KPoly<K> operator -(KPoly<K> a, int b) => a + (-b);
    public static KPoly<K> operator *(KPoly<K> a, KPoly<K> b) => a.Mul(b);
    public static KPoly<K> operator *(KPoly<K> a, int b) => a.Mul(b);
    public static KPoly<K> operator *(int a, KPoly<K> b) => b.Mul(a);
    public static KPoly<K> operator /(KPoly<K> a, KPoly<K> b) => a.Div(b).quo;
    public static KPoly<K> operator /(KPoly<K> a, int b) => a.Div(a.One.Mul(b)).quo;

    public static KPoly<K> operator +(KPoly<K> a, K b) => a + a.One.KMul(b);
    public static KPoly<K> operator +(K a, KPoly<K> b) => b.One.KMul(a) + b;
    public static KPoly<K> operator -(KPoly<K> a, K b) => a - a.One.KMul(b);
    public static KPoly<K> operator -(K a, KPoly<K> b) => b.One.KMul(a) - b;
    public static KPoly<K> operator *(KPoly<K> a, K b) => a.KMul(b);
    public static KPoly<K> operator *(K a, KPoly<K> b) => b.KMul(a);
    public static KPoly<K> operator /(KPoly<K> a, K b) => a.KMul(b.Inv());
}
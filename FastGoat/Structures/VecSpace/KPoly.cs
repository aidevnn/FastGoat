using FastGoat.Commons;
using FastGoat.UserGroup.Integers;

namespace FastGoat.Structures.VecSpace;

public readonly struct KPoly<K> : IVsElt<K, KPoly<K>>, IElt<KPoly<K>>, IRingElt<KPoly<K>>
    where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
{
    public K[] Coefs { get; }
    public int Degree { get; }
    public K KZero { get; }
    public K KOne { get; }
    public int P { get; }
    public char x { get; }

    /// <summary>
    /// Return P(X) = X in Ring Z
    /// </summary>
    public KPoly(char x0)
    {
        KZero = new K().Zero;
        KOne = new K().One;
        P = 0;
        x = x0;
        Coefs = new[] { KZero, KOne };
        Hash = Coefs.Aggregate(0, (acc, a) => (acc, a.Hash).GetHashCode());
        Degree = Coefs.Length - 1;
    }

    public KPoly(char x0, K k0)
    {
        KZero = k0.Zero;
        KOne = k0.One;
        P = k0.P;
        x = x0;
        Coefs = new[] { k0 };
        Hash = Coefs.Aggregate(0, (acc, a) => (acc, a.Hash).GetHashCode());
        Degree = Coefs.Length - 1;
    }

    public KPoly(char x0, K kZero, K[] coefs)
    {
        KZero = kZero;
        KOne = kZero.One;
        P = kZero.P;
        x = x0;

        Coefs = coefs.Length != 0 ? coefs : new[] { kZero };
        Hash = coefs.Aggregate(0, (acc, a) => (acc, a.Hash).GetHashCode());
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

    public bool Equals(KPoly<K> other) => Hash == other.Hash;

    public int CompareTo(KPoly<K> other)
    {
        var compDegree = Degree.CompareTo(other.Degree);
        if (compDegree != 0)
            return compDegree;

        for (int i = Degree; i >= 0; i--)
        {
            var comp = this[i].CompareTo(other[i]);
            if (comp != 0)
                return comp;
        }

        return 0;
    }

    public int Hash { get; }

    public bool IsZero() => Degree == 0 && Coefs[0].IsZero();

    public KPoly<K> Zero => new(x, KZero, new[] { KZero });
    public KPoly<K> One => new(x, KZero, new[] { KOne });
    public KPoly<K> X => new(x, KZero, new[] { KZero, KOne });

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

    public KPoly<K> Monic => KMul(Coefs.Last().Inv());

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
            return (Zero, this);

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
        string Str(K e, int i, char x0)
        {
            if (i == 0) return $"{e}";

            if (e.Equals(e.One))
            {
                if (i == 1) return $"{x0}";
                return $"{x0}^{i}";
            }

            if (e.P == 0 && e.Equals(e.One.Opp()))
            {
                if (i == 1) return $"-{x0}";
                return $"-{x0}^{i}";
            }

            if (i == 1) return $"{e}*{x0}";
            return $"{e}*{x0}^{i}";
        }

        if (IsZero())
            return "0";

        var x0 = x;
        return Coefs.Select((e, i) => (e, i)).Where(e => !e.e.IsZero())
            .Select(e => Str(e.e, e.i, x0)).Glue(" + ");
    }

    public static KPoly<K> operator +(KPoly<K> a, KPoly<K> b) => a.Add(b);
    public static KPoly<K> operator -(KPoly<K> a, KPoly<K> b) => a.Sub(b);
    public static KPoly<K> operator -(KPoly<K> a) => a.Opp();
    public static KPoly<K> operator *(KPoly<K> a, KPoly<K> b) => a.Mul(b);
    public static KPoly<K> operator /(KPoly<K> a, KPoly<K> b) => a.Div(b).quo;

    public static KPoly<K> operator +(KPoly<K> a, K b) => a.Add(new(a.x, b));
    public static KPoly<K> operator +(K b, KPoly<K> a) => a.Add(new(a.x, b));
    public static KPoly<K> operator +(KPoly<K> a, int b) => a.Add(new(a.x, a.KOne.Mul(b)));
    public static KPoly<K> operator +(int b, KPoly<K> a) => a.Add(new(a.x, a.KOne.Mul(b)));

    public static KPoly<K> operator -(KPoly<K> a, K b) => a.Sub(new(a.x, b));
    public static KPoly<K> operator -(K b, KPoly<K> a) => new KPoly<K>(a.x, b).Sub(a);
    public static KPoly<K> operator -(KPoly<K> a, int b) => a.Sub(new(a.x, a.KOne.Mul(b)));
    public static KPoly<K> operator -(int b, KPoly<K> a) => new KPoly<K>(a.x, a.KOne.Mul(b)).Sub(a);

    public static KPoly<K> operator *(KPoly<K> a, K b) => a.KMul(b);
    public static KPoly<K> operator *(K b, KPoly<K> a) => a.KMul(b);
    public static KPoly<K> operator *(KPoly<K> a, int b) => a.Mul(b);
    public static KPoly<K> operator *(int b, KPoly<K> a) => a.Mul(b);

    public static KPoly<K> operator /(KPoly<K> a, K b) => a.KMul(b.Inv());
    public static KPoly<K> operator /(KPoly<K> a, int b) => a.KMul(a.KOne.Mul(b).Inv());
}
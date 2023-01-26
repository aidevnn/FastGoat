using FastGoat.Commons;
using Microsoft.VisualBasic;

namespace FastGoat.Structures.VecSpace;

public readonly struct Polynomial<K, T> : IVsElt<K, Polynomial<K, T>>, IElt<Polynomial<K, T>>,
    IRingElt<Polynomial<K, T>>
    where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    where T : struct, IElt<T>
{
    public SortedList<Monom<T>, K> Coefs { get; }

    public Indeterminates<T> Indeterminates { get; }

    public K KZero { get; }
    public K KOne => KZero.One;

    public Polynomial(Indeterminates<T> indeterminates, K scalar)
    {
        Indeterminates = indeterminates;
        KZero = scalar.Zero;
        Coefs = new() { [new(indeterminates)] = scalar };
        Hash = Coefs.Aggregate(0, (acc, a) => (acc, (a.Key.Hash, a.Value.Hash)).GetHashCode());
    }

    public Polynomial(T x, K scalar)
    {
        Indeterminates = new Indeterminates<T>(x);
        KZero = scalar.Zero;
        Coefs = new() { [new(Indeterminates, x)] = scalar };
        Hash = Coefs.Aggregate(0, (acc, a) => (acc, (a.Key.Hash, a.Value.Hash)).GetHashCode());
    }

    public Polynomial(Monom<T> m, K scalar)
    {
        Indeterminates = m.Indeterminates;
        KZero = scalar.Zero;
        Coefs = new() { [m] = scalar };
        Hash = Coefs.Aggregate(0, (acc, a) => (acc, (a.Key.Hash, a.Value.Hash)).GetHashCode());
    }

    public Polynomial(Indeterminates<T> indeterminates, K zero, SortedList<Monom<T>, K> coefs)
    {
        if (!zero.IsZero())
            throw new GroupException(GroupExceptionType.GroupDef);

        Indeterminates = indeterminates;
        KZero = zero;
        Coefs = coefs.Count != 0 ? coefs : new() { [new(Indeterminates)] = zero };
        Hash = Coefs.Aggregate(0, (acc, a) => (acc, (a.Key.Hash, a.Value.Hash)).GetHashCode());
    }

    public bool Equals(Polynomial<K, T> other) => Hash == other.Hash;

    public int CompareTo(Polynomial<K, T> other)
    {
        return Coefs.OrderByDescending(kp => kp.Key).Select(kp => (kp.Key, kp.Value))
            .SequenceCompareTo(other.Coefs.OrderByDescending(kp => kp.Key).Select(kp => (kp.Key, kp.Value)));
    }

    public int P => KZero.P;

    public int Hash { get; }

    public bool IsZero()
    {
        var zero = new Monom<T>(Indeterminates);
        foreach (var (key, value) in Coefs)
        {
            if (!key.Equals(zero) || !value.IsZero())
                return false;
        }

        return true;
    }

    public Polynomial<K, T> Zero => new(Indeterminates, KZero);
    public Polynomial<K, T> One => new(Indeterminates, KZero.One);

    public T ExtractIndeterminate
    {
        get
        {
            if (NbIndeterminates != 1)
                throw new ArgumentException("Not univariate polynomial");

            var coefs = Coefs.Keys.First(m => m.Degree > 0);
            return Indeterminates.First(xi => coefs[xi] > 0);
        }
    }

    public Polynomial<K, T> Substitute(Polynomial<K, T> f, T xi)
    {
        var poly = Zero;
        foreach (var (m, c) in Coefs)
        {
            var (n, m1) = m.Remove(xi);
            var p0 = new Polynomial<K, T>(m1, c);
            poly += p0 * f.Pow(n);
        }

        return poly;
    }

    public Polynomial<K, T> Substitute(Polynomial<K, T> f, Polynomial<K, T> xi) => Substitute(f, xi.ExtractIndeterminate);

    public Polynomial<K, T> Substitute(K k, T xi) => Substitute(k * One, xi);
    public Polynomial<K, T> Substitute(K k, Polynomial<K, T> xi) => Substitute(k * One, xi);

    public Polynomial<K, T> Substitute(int k, T xi) => Substitute(k * One, xi);
    public Polynomial<K, T> Substitute(int k, Polynomial<K, T> xi) => Substitute(k * One, xi);
    public Indeterminates<T> Xi => Indeterminates;

    public (K lc, Monom<T> lm, Polynomial<K, T> lt) LeadingDetails
    {
        get
        {
            if (IsZero())
                throw new ArgumentException();

            var kp = Coefs.Last();
            var lt = new SortedList<Monom<T>, K>() { [kp.Key] = kp.Value };
            return (kp.Value, kp.Key, new(Indeterminates, KZero, lt));
        }
    }

    public Polynomial<K, T> Add(Polynomial<K, T> e)
    {
        SortedList<Monom<T>, K> coefs = new(Coefs);
        foreach (var kp in e.Coefs)
        {
            coefs[kp.Key] = coefs.ContainsKey(kp.Key)
                ? coefs[kp.Key].Add(kp.Value)
                : kp.Value;
        }

        foreach (var a in coefs.Where(a => a.Value.IsZero()).ToArray())
            coefs.Remove(a.Key);

        return new(Indeterminates, KZero, coefs);
    }

    public Polynomial<K, T> D(T t)
    {
        var derivate = Zero;
        foreach (var coef in Coefs)
        {
            var (n, m) = coef.Key.D(t);
            if (n != 0)
            {
                derivate = derivate.Add(new(Indeterminates, KZero, new() { [m] = coef.Value.Mul(n) }));
            }
        }

        return derivate;
    }

    public Polynomial<K, T> D(Polynomial<K, T> m)
    {
        if (m.Coefs.Count != 1 || m.Coefs.First().Key.Degree != 1)
            throw new GroupException(GroupExceptionType.GroupDef);

        return D(m.Coefs.First().Key.Indeterminates.First());
    }

    public int DegreeOf(T t) => Coefs.Max(e => e.Key[t]);

    public Polynomial<K, T> Monic() => IsZero() ? this : this.KMul(LeadingDetails.lc.Inv());

    public Polynomial<K, T> CoefMax(T t)
    {
        var n = DegreeOf(t);
        var xi = new Monom<T>(Indeterminates, t, n);
        var k0 = KZero;
        var coefMax = Coefs.Select(e => (e, Monom<T>.Reduce(e.Key, xi))).Where(e => e.Item2.pa.IsOne)
            .Select(e => new Polynomial<K, T>(xi.Indeterminates, k0, new() { [e.Item2.pb] = e.e.Value }))
            .Aggregate(Zero, (acc, a) => acc + a);

        return coefMax;
    }

    public Polynomial<K, T> Sub(Polynomial<K, T> e)
    {
        SortedList<Monom<T>, K> coefs = new(Coefs);
        foreach (var kp in e.Coefs)
        {
            coefs[kp.Key] = coefs.ContainsKey(kp.Key)
                ? coefs[kp.Key].Sub(kp.Value)
                : kp.Value.Opp();
        }

        foreach (var a in coefs.Where(a => a.Value.IsZero()).ToArray())
            coefs.Remove(a.Key);

        return new(Indeterminates, KZero, coefs);
    }

    public Polynomial<K, T> Opp()
    {
        SortedList<Monom<T>, K> coefs = new(Coefs.ToDictionary(a => a.Key, a => a.Value.Opp()));
        return new(Indeterminates, KZero, coefs);
    }

    public Polynomial<K, T> Mul(Polynomial<K, T> e)
    {
        SortedList<Monom<T>, K> coefs = new();
        foreach (var m0 in Coefs)
        {
            foreach (var m1 in e.Coefs)
            {
                var m2 = m0.Key.Mul(m1.Key);
                var a = m0.Value.Mul(m1.Value);

                coefs[m2] = coefs.ContainsKey(m2)
                    ? coefs[m2].Add(a)
                    : a;
            }
        }

        foreach (var a in coefs.Where(a => a.Value.IsZero()).ToArray())
            coefs.Remove(a.Key);

        return new(Indeterminates, KZero, coefs);
    }

    private Polynomial<K, T> InPlaceSubMul(Polynomial<K, T> e, K c, Monom<T> m)
    {
        SortedList<Monom<T>, K> coefs = new(Coefs);
        foreach (var kp in e.Coefs)
        {
            var key = kp.Key.Mul(m);
            coefs[key] = coefs.ContainsKey(key)
                ? coefs[key].Sub(kp.Value * c)
                : -kp.Value * c;
        }

        foreach (var a in coefs.Where(a => a.Value.IsZero()).ToArray())
            coefs.Remove(a.Key);

        return new(Indeterminates, KZero, coefs);
    }

    public (Polynomial<K, T> quo, Polynomial<K, T> rem) Div(Polynomial<K, T> e)
    {
        if (e.IsZero())
            throw new DivideByZeroException();

        // if (e.Coefs.Keys.SelectMany(a => a.ContentIndeterminates).Distinct().Count() > 1)
        //     throw new GroupException(GroupExceptionType.GroupDef);

        var rem = new Polynomial<K, T>(Indeterminates, KZero, new(Coefs));
        var (elm, elc) = e.Coefs.Last();
        int k = rem.Coefs.Count - 1;
        SortedList<Monom<T>, K> quo = new();
        while (k >= 0)
        {
            var (alm, alc) = (rem.Coefs.Keys[k], rem.Coefs.Values[k]);
            // var mnm = alm.Div(elm);
            var (mnm0, mnm1) = Monom<T>.Reduce(alm, elm);
            if (!mnm0.IsOne)
            {
                --k;
                continue;
            }

            var m = mnm1; //mnm.Value;
            var (q, r) = alc.Div(elc);
            if (!r.IsZero())
                break;

            rem = rem.InPlaceSubMul(e, q, m);
            k = rem.Coefs.Count - 1;
            quo[m] = q;
            if (rem.IsZero())
                break;
        }

        return (new(Indeterminates, KZero, quo), rem);
    }

    public Polynomial<K, T> KMul(K k)
    {
        SortedList<Monom<T>, K> coefs = new(Coefs.ToDictionary(a => a.Key, a => a.Value.Mul(k)));
        foreach (var a in coefs.Where(a => a.Value.IsZero()).ToArray())
            coefs.Remove(a.Key);

        return new(Indeterminates, KZero, coefs);
    }

    public Polynomial<K, T> Mul(int k)
    {
        if (k == 0)
            return new(Indeterminates, KZero);

        SortedList<Monom<T>, K> coefs = new(Coefs.ToDictionary(a => a.Key, a => a.Value.Mul(k)));
        return new(Indeterminates, KZero, coefs);
    }

    public Polynomial<K, T> Pow(int k)
    {
        if (k < 0)
            throw new GroupException(GroupExceptionType.GroupDef);

        var pi = this;
        return Enumerable.Repeat(pi, k).Aggregate(One, (a, b) => a.Mul(b));
    }

    public override int GetHashCode() => Hash;

    public string GetString(bool reverse = false)
    {
        var one = KZero.One;
        var sep = (Monom.Display & MonomDisplay.Star) == MonomDisplay.Star ? "*" : "·";

        string Str(Monom<T> m, K k)
        {
            var sm = $"{m}";
            if (k.Equals(one))
                return string.IsNullOrEmpty(sm) ? "1" : sm;

            if (k.Equals(one.Opp()) && one.P == 0)
                return string.IsNullOrEmpty(sm) ? "-1" : $"-{sm}";

            var k0 = $"{k}";
            k0 = k0.Contains('+') ? $"({k0})" : k0;
            return string.IsNullOrEmpty(sm) ? $"{k}" : $"{k0}{sep}{sm}";
        }

        if (reverse)
            return Coefs.Reverse().Select(kp => Str(kp.Key, kp.Value)).Glue(" + ");

        return Coefs.Select(kp => Str(kp.Key, kp.Value)).Glue(" + ");
    }

    public int NbIndeterminates
    {
        get
        {
            var nb = 0;
            foreach (var x in Indeterminates)
            {
                if (Coefs.Keys.Any(m => m[x] != 0))
                    ++nb;
            }

            return nb;
        }
    }

    public K this[Monom<T> idx]
    {
        get
        {
            if (!Coefs.ContainsKey(idx))
                return KZero;

            return Coefs[idx];
        }
    }

    public override string ToString() => GetString(true);

    public static Polynomial<K, T> operator +(Polynomial<K, T> a, Polynomial<K, T> b) => a.Add(b);
    public static Polynomial<K, T> operator +(int a, Polynomial<K, T> b) => b.Add(b.One.Mul(a));
    public static Polynomial<K, T> operator +(Polynomial<K, T> a, int b) => a.Add(a.One.Mul(b));
    public static Polynomial<K, T> operator -(Polynomial<K, T> a) => a.Opp();
    public static Polynomial<K, T> operator -(Polynomial<K, T> a, Polynomial<K, T> b) => a + (-b);
    public static Polynomial<K, T> operator -(int a, Polynomial<K, T> b) => a + (-b);
    public static Polynomial<K, T> operator -(Polynomial<K, T> a, int b) => a + (-b);
    public static Polynomial<K, T> operator *(Polynomial<K, T> a, Polynomial<K, T> b) => a.Mul(b);
    public static Polynomial<K, T> operator *(Polynomial<K, T> a, int b) => a.Mul(b);
    public static Polynomial<K, T> operator *(int a, Polynomial<K, T> b) => b.Mul(a);
    public static Polynomial<K, T> operator /(Polynomial<K, T> a, Polynomial<K, T> b) => a.Div(b).quo;
    public static Polynomial<K, T> operator /(Polynomial<K, T> a, int b) => a.Div(a.One.Mul(b)).quo;

    public static Polynomial<K, T> operator +(Polynomial<K, T> a, K b) => a + a.One.KMul(b);
    public static Polynomial<K, T> operator +(K a, Polynomial<K, T> b) => b.One.KMul(a) + b;
    public static Polynomial<K, T> operator -(Polynomial<K, T> a, K b) => a - a.One.KMul(b);
    public static Polynomial<K, T> operator -(K a, Polynomial<K, T> b) => b.One.KMul(a) - b;
    public static Polynomial<K, T> operator *(Polynomial<K, T> a, K b) => a.KMul(b);
    public static Polynomial<K, T> operator *(K a, Polynomial<K, T> b) => b.KMul(a);
    public static Polynomial<K, T> operator /(Polynomial<K, T> a, K b) => a.KMul(b.Inv());
}
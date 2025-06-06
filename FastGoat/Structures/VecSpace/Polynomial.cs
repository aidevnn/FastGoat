using FastGoat.Commons;
using FastGoat.UserGroup.Polynoms;

namespace FastGoat.Structures.VecSpace;

public readonly struct Polynomial<K, T> : IElt<Polynomial<K, T>>, IRingElt<Polynomial<K, T>>,
    IFieldElt<Polynomial<K, T>>, IModuleElt<K, Polynomial<K, T>>, IVsElt<K, Polynomial<K, T>>
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

    public bool Equals(Polynomial<K, T> other)
    {
        // if (Hash != other.Hash)
        //     return false;

        return Coefs.Count == other.Coefs.Count &&
               Coefs.All(e => other.Coefs.ContainsKey(e.Key) && other.Coefs[e.Key].Equals(e.Value));
    }

    public int CompareTo(Polynomial<K, T> other)
    {
        return Coefs.OrderByDescending(kp => kp.Key).Select(kp => (kp.Key, kp.Value))
            .SequenceCompareTo(other.Coefs.OrderByDescending(kp => kp.Key).Select(kp => (kp.Key, kp.Value)));
    }

    public int P => KZero.P;

    public Polynomial<K, T> Recreate(Indeterminates<T> ind)
    {
        var coefs = Coefs.Where(e => !e.Value.IsZero()).ToDictionary(e => new Monom<T>(ind, e.Key), e => e.Value);
        return new(ind, KZero, new(coefs));
    }

    public Polynomial<K, T> Recreate() => Recreate(Indeterminates);

    public Polynomial<K, T>[] Variables
    {
        get
        {
            var list = new List<Polynomial<K, T>>();
            foreach (var xi in ExtractAllIndeterminates)
                list.Add(X(xi));

            return list.ToArray();
        }
    }

    public Polynomial<K, T>[] AllVariables
    {
        get
        {
            var list = new List<Polynomial<K, T>>();
            foreach (var xi in Indeterminates)
                list.Add(X(xi));

            return list.ToArray();
        }
    }

    public (T, Polynomial<K, T>)[] IndeterminatesAndVariables
    {
        get
        {
            var list = new List<(T, Polynomial<K, T>)>();
            foreach (var xi in ExtractAllIndeterminates)
                list.Add((xi, X(xi)));

            return list.ToArray();
        }
    }

    public (T, Polynomial<K, T>)[] AllIndeterminatesAndVariables
    {
        get
        {
            var list = new List<(T, Polynomial<K, T>)>();
            foreach (var xi in Indeterminates)
                list.Add((xi, X(xi)));

            return list.ToArray();
        }
    }

    public Polynomial<K, T> X(T t)
    {
        if (!Indeterminates.Contains(t))
            throw new();

        return new Polynomial<K, T>(new Monom<T>(Indeterminates, t), KOne);
    }

    public Polynomial<K, T> XY(Monom<T> mn)
    {
        var ind = Indeterminates;
        if (mn.ContentIndeterminates.Any(t => !ind.Contains(t)))
            throw new();

        return new Polynomial<K, T>(mn, KOne);
    }

    public Polynomial<K, T> Inv()
    {
        if (Coefs.Count == 1)
            return new(Indeterminates, ConstTerm.Inv());

        throw new();
    }

    public bool Invertible() => false;

    public static Polynomial<K, T> operator /(int a, Polynomial<K, T> b) => b.Inv() * a;

    public static double Abs(Polynomial<K, T> t) => throw new();

    public static bool IsValuedField => false;

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
                throw new ArgumentException($"Not univariate polynomial P={this}");

            var coefs = Coefs.Keys.First(m => m.Degree > 0);
            return Indeterminates.First(xi => coefs[xi] > 0);
        }
    }

    public T[] ExtractAllIndeterminates
    {
        get
        {
            return Coefs.Keys.SelectMany(m => m.ContentIndeterminatesMonoms).Distinct().Order()
                .SelectMany(e => e.ContentIndeterminates).ToArray();
        }
    }

    public Monom<T> ExtractMonom => new Monom<T>(Indeterminates, ExtractIndeterminate);

    public Polynomial<K, T> Substitute(Polynomial<K, T> f, T xi)
    {
        var poly = Zero;
        foreach (var (m, c) in Coefs)
        {
            var (n, m1) = m.Remove(xi);
            var p0 = new Polynomial<K, T>(m1, c);
            var fn = f.Pow(n);
            var p0fn = p0 * fn;
            poly = poly + p0fn;
        }

        return poly;
    }

    public Polynomial<K, T> Substitute(params (Polynomial<K, T> f, T xi)[] subs)
    {
        return subs.Aggregate(this, (acc, e) => acc.Substitute(e.f, e.xi));
    }

    public Polynomial<K, T> Substitute(List<(T, Polynomial<K, T>)> dico)
    {
        var poly = Zero;
        foreach (var (m, c) in Coefs)
        {
            var poly0 = One;
            var m0 = m;
            foreach (var (xi, f) in dico)
            {
                (int n, m0) = m0.Remove(xi);
                poly0 *= f.Pow(n);
            }

            var p0 = new Polynomial<K, T>(m0, c);
            poly = poly + poly0 * p0;
        }

        return poly;
    }

    public Polynomial<U, T> Substitute<U>(U f, T xi)
        where U : struct, IElt<U>, IRingElt<U>, IFieldElt<U>, IModuleElt<K, U>, IVsElt<K, U>
    {
        var poly = new Polynomial<U, T>(Indeterminates, f.Zero);
        foreach (var (m, c) in Coefs)
        {
            var (n, m1) = m.Remove(xi);
            var p0 = new Polynomial<U, T>(m1, c * f.One);
            poly += p0 * f.Pow(n);
        }

        return poly;
    }

    public EPolynomial<K> Substitute(EPolynomial<K> f, T xi)
    {
        var poly = f.Zero;
        foreach (var (m, c) in Coefs)
        {
            var (n, m1) = m.Remove(xi);

            if (m1.IsOne)
                poly += c * f.Pow(n);
            else if (m1 is Monom<Xi> m2)
            {
                var p0 = new EPolynomial<K>(new Polynomial<K, Xi>(m2, c), f.Basis);
                poly += p0 * f.Pow(n);
            }
            else
                throw new Exception();
        }

        return poly;
    }

    public Polynomial<K, T> Substitute(Polynomial<K, T> f, Polynomial<K, T> xi) =>
        Substitute(f, xi.ExtractIndeterminate);

    public Polynomial<K, T> Substitute(K k, T xi) => Substitute(k * One, xi);

    public Polynomial<K, T> Substitute(K k, Polynomial<K, T> xi) =>
        Substitute(new Polynomial<K, T>(Indeterminates, k), xi);

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

    public K ConstTerm => this[new(Indeterminates)];

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

        return D(m.ExtractIndeterminate);
    }

    public int DegreeOf(T t) => Coefs.Max(e => e.Key[t]);
    public int Degree => Coefs.Max(e => e.Key.Degree);

    public Polynomial<K, T> Monic() => IsZero() ? this : KMul(LeadingDetails.lc.Inv());

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

    private void InPlaceSubMul(SortedList<Monom<T>, K> A, SortedList<Monom<T>, K> B, K c, Monom<T> m)
    {
        foreach (var kp in B)
        {
            var key = kp.Key.Mul(m);
            A[key] = A.ContainsKey(key)
                ? A[key].Sub(kp.Value * c)
                : -kp.Value * c;
        }

        foreach (var a in A.Where(a => a.Value.IsZero()).ToArray())
            A.Remove(a.Key);

        if (A.Count == 0)
            A[m.One] = c.Zero;
    }

    public (Polynomial<K, T> quo, Polynomial<K, T> rem) Div(Polynomial<K, T> e)
    {
        if (e.IsZero())
            throw new DivideByZeroException();

        var rem = new SortedList<Monom<T>, K>(Coefs);
        var (elm, elc) = e.Coefs.Last();
        int k = rem.Count - 1;
        SortedList<Monom<T>, K> quo = new();
        var mn0 = new Monom<T>(Indeterminates);
        while (k >= 0)
        {
            var (alm, alc) = (rem.Keys[k], rem.Values[k]);
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

            InPlaceSubMul(rem, e.Coefs, q, m);
            k = rem.Count - 1;
            quo[m] = q;
            if (rem.Count == 1 && rem.ContainsKey(mn0) && rem[mn0].IsZero())
                break;
        }

        return (new(Indeterminates, KZero, quo), new(Indeterminates, KZero, rem));
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

        return this.FastPow(k);
    }

    public override int GetHashCode() => Hash;

    public string GetString()
    {
        var one = KZero.One;
        var sep = (Ring.DisplayPolynomial & MonomDisplay.Star) == MonomDisplay.Star ? "*" : "·";

        string Str(Monom<T> m, K k)
        {
            var sm = $"{m}";
            if (k.Equals(one))
                return string.IsNullOrEmpty(sm) ? "1" : sm;

            if (k.Equals(one.Opp()) && (one.P == 0 || $"{k}".Equals("-1")))
                return string.IsNullOrEmpty(sm) ? "-1" : $"-{sm}";

            var k0 = $"{k}";
            k0 = k0.Contains('+') || (k0.Any(char.IsLetter) && k0.Contains('-')) ? $"({k0})" : k0;
            return string.IsNullOrEmpty(sm) ? $"{k}" : $"{k0}{sep}{sm}";
        }

        return Coefs.OrderByDescending(e => e.Key).Select(kp => Str(kp.Key, kp.Value)).Glue(" + ");
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

    public override string ToString() => GetString().Replace("+ -", "- ");

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
using FastGoat.Commons;

namespace FastGoat.Structures.VecSpace;

public readonly struct Polynomial<K, T> : IVsElt<K, Polynomial<K, T>>, IElt<Polynomial<K, T>>, IRingElt<Polynomial<K, T>>
    where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    where T : struct, IElt<T>
{
    private SortedList<Monom<T>, K> Coefs { get; }
    public IEnumerable<T> Indeterminates { get; }
    public K KZero { get; }

    public Polynomial(IEnumerable<T> indeterminates, K zero)
    {
        if (!zero.IsZero() || !indeterminates.Any())
            throw new GroupException(GroupExceptionType.GroupDef);

        KZero = zero;
        Indeterminates = indeterminates;
        Coefs = new() { [new()] = zero };
        Hash = Coefs.Aggregate(0, (acc, a) => (acc, (a.Key.Hash, a.Value.Hash)).GetHashCode());
    }

    private Polynomial(IEnumerable<T> indeterminates, K zero, SortedList<Monom<T>, K> coefs)
    {
        if (!zero.IsZero())
            throw new GroupException(GroupExceptionType.GroupDef);

        KZero = zero;
        Indeterminates = indeterminates;
        Coefs = coefs.Count != 0 ? coefs : new() { [new()] = zero };
        Hash = Coefs.Aggregate(0, (acc, a) => (acc, (a.Key.Hash, a.Value.Hash)).GetHashCode());
    }

    public Polynomial<K, T>[] Xi()
    {
        List<Polynomial<K, T>> list = new();
        foreach (var xi in Indeterminates)
        {
            var m = new Monom<T>(xi);
            var poly = new Polynomial<K, T>(Indeterminates, KZero, new() { [m] = KZero.One });
            list.Add(poly);
        }
        return list.ToArray();
    }

    public bool Equals(Polynomial<K, T> other) => Hash == other.Hash;

    public int CompareTo(Polynomial<K, T> other)
    {
        return Coefs.Select(kp => (kp.Key, kp.Value)).SequenceCompareTo(other.Coefs.Select(kp => (kp.Key, kp.Value)));
    }

    public int P => 0;

    public int Hash { get; }

    public bool IsZero() => Coefs.All(a=>a.Key.Equals(new()) && a.Value.IsZero());

    public Polynomial<K, T> Zero => new(Indeterminates, KZero);
    public Polynomial<K, T> One => new(Indeterminates, KZero, new() { [new()] = KZero.One });

    public Polynomial<K, T> Add(Polynomial<K, T> e)
    {
        if (!Indeterminates.Equals(e.Indeterminates))
            throw new GroupException(GroupExceptionType.GroupDef);
        
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

    public Polynomial<K, T> Sub(Polynomial<K, T> e)
    {
        if (!Indeterminates.Equals(e.Indeterminates))
            throw new GroupException(GroupExceptionType.GroupDef);

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
        SortedList<Monom<T>, K> coefs = new(Coefs.ToDictionary(a=>a.Key,a=>a.Value.Opp()));
        return new(Indeterminates, KZero, coefs);
    }

    public Polynomial<K, T> Mul(Polynomial<K, T> e)
    {
        if (!Indeterminates.Equals(e.Indeterminates))
            throw new GroupException(GroupExceptionType.GroupDef);

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

    public (Polynomial<K, T> quo, Polynomial<K, T> rem) Div(Polynomial<K, T> e)
    {
        if (!Indeterminates.Equals(e.Indeterminates))
            throw new GroupException(GroupExceptionType.GroupDef);

        if (e.IsZero())
            throw new DivideByZeroException();

        if (e.Coefs.Keys.SelectMany(a => a.Indeterminates).Distinct().Count() > 1)
            throw new GroupException(GroupExceptionType.GroupDef);
        
        var rem = new Polynomial<K, T>(Indeterminates, KZero, new(Coefs));
        var coefs = new Stack<KeyValuePair<Monom<T>, K>>(rem.Coefs.Reverse());
        SortedList<Monom<T>, K> quo = new();
        var em = e.Coefs.First();
        while (coefs.Count != 0)
        {
            var am = coefs.Pop();
            if (am.Key.Degree < em.Key.Degree)
                break;

            var mnm = am.Key.Div(em.Key);
            if(!mnm.HasValue)
                continue;

            var m = mnm.Value;
            var qr = am.Value.Div(em.Value);
            if (!qr.rem.IsZero())
                break;

            var p = new Polynomial<K, T>(Indeterminates, KZero, new() { [m] = qr.quo });
            rem = rem.Sub(e.Mul(p));
            quo[m] = qr.quo;
            coefs = new Stack<KeyValuePair<Monom<T>, K>>(rem.Coefs.Reverse());
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
        SortedList<Monom<T>, K> coefs = new(Coefs.ToDictionary(a => a.Key, a => a.Value.Mul(k)));
        foreach (var a in coefs.Where(a => a.Value.IsZero()).ToArray())
            coefs.Remove(a.Key);
        
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

    public override string ToString()
    {
        var one = KZero.One;
        string Str(Monom<T> m, K k)
        {
            var sm = $"{m}";
            if (k.Equals(one))
                return string.IsNullOrEmpty(sm) ? "1" : sm;
            
            if (k.Equals(one.Opp()))
                return string.IsNullOrEmpty(sm) ? "-1" : $"-{sm}";

            return string.IsNullOrEmpty(sm) ?$"{k}" : $"{k}·{sm}";
        }

        return Coefs.Select(kp => Str(kp.Key, kp.Value)).Glue(" + ");
    }

    public static Polynomial<K, T> operator +(Polynomial<K, T> a, Polynomial<K, T> b) => a.Add(b);
    public static Polynomial<K, T> operator -(Polynomial<K, T> a, Polynomial<K, T> b) => a.Sub(b);
    public static Polynomial<K, T> operator *(Polynomial<K, T> a, Polynomial<K, T> b) => a.Mul(b);
    public static Polynomial<K, T> operator /(Polynomial<K, T> a, Polynomial<K, T> b) => a.Div(b).quo;
    
    public static Polynomial<K, T> operator -(Polynomial<K, T> a) => a.Opp();
    public static Polynomial<K, T> operator *(int a, Polynomial<K, T> b) => b.Mul(a);
    public static Polynomial<K, T> operator /(Polynomial<K, T> a, int k) => a.KMul(a.KZero.One.Mul(k).Inv());
    
    public static Polynomial<K, T> operator +(Polynomial<K, T> a, int b) => a.Add(a.One.Mul(b));
    public static Polynomial<K, T> operator +(int b, Polynomial<K, T> a) => a.Add(a.One.Mul(b));
    public static Polynomial<K, T> operator -(Polynomial<K, T> a, int b) => a.Sub(a.One.Mul(b));
    public static Polynomial<K, T> operator -(int b, Polynomial<K, T> a) => a.One.Mul(b).Sub(a);

}
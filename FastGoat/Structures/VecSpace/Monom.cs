using FastGoat.Commons;

namespace FastGoat.Structures.VecSpace;

[Flags]
public enum MonomOrder
{
    Lex = 0b001,
    Graded = 0b010,
    Reverse = 0b100,
    GrLex = Graded | Lex,
    RevLex = Reverse | Lex,
    GrevLex = Graded | Reverse
}

[Flags]
public enum MonomDisplay
{
    Caret = 0b00001,
    Superscript = 0b00010,
    PowFct = 0b00100,
    Star = 0b01000,
    Dot = 0b10000,
    Default = Superscript,
    StarSuperscript = Star | Superscript,
    StarCaret = Star | Caret,
    StarPowFct = Star | PowFct,
    DotSuperscript = Dot | Superscript,
    DotCaret = Dot | Caret,
    DotPowFct = Dot | PowFct
}

public readonly struct Monom<T> : IElt<Monom<T>> where T : struct, IElt<T>
{
    private Dictionary<T, int> Content { get; }
    public Indeterminates<T> Indeterminates { get; }
    public IEnumerable<T> ContentIndeterminates => Content.Keys;
    public IEnumerable<Monom<T>> ContentIndeterminatesMonoms
    {
        get
        {
            var monom = this;
            return Content.Keys.Select(t => new Monom<T>(monom.Indeterminates, t, 1));
        }
    }

    public int Degree { get; }

    public Monom()
    {
        throw new ArgumentException();
    }

    public Monom(Indeterminates<T> indeterminates)
    {
        Content = new();
        Degree = 0;
        Indeterminates = indeterminates;
        Hash = Indeterminates.Hash;
    }

    public Monom(Indeterminates<T> indeterminates, T s)
    {
        if (!indeterminates.Contains(s))
            throw new ArgumentException();

        Indeterminates = indeterminates;
        Content = new() { [s] = 1 };
        Degree = 1;
        Hash = (Indeterminates.Hash, (s.Hash, 1)).GetHashCode();
    }

    public Monom(Indeterminates<T> indeterminates, T s, int n)
    {
        if (!indeterminates.Contains(s))
            throw new ArgumentException();

        if (n < 0)
            throw new GroupException(GroupExceptionType.BaseGroup);
        else if (n == 0)
        {
            Content = new();
            Degree = 0;
            Indeterminates = indeterminates;
            Hash = Indeterminates.Hash;
        }
        else
        {
            Indeterminates = indeterminates;
            Content = new() { [s] = n };
            Degree = n;
            Hash = (Indeterminates.Hash, (s.Hash, n)).GetHashCode();
        }
    }

    public Monom(T s, int n)
    {
        if (n < 1)
            throw new GroupException(GroupExceptionType.BaseGroup);

        Indeterminates = new(s);
        Content = new() { [s] = n };
        Degree = n;
        Hash = (Indeterminates.Hash, (s.Hash, n)).GetHashCode();
    }

    private Monom(Indeterminates<T> indeterminates, Dictionary<T, int> content)
    {
        Indeterminates = indeterminates;
        Content = content;
        Degree = content.Sum(e => e.Value);
        Hash = Indeterminates.Hash;
        foreach (var t in Indeterminates)
        {
            if (!content.ContainsKey(t))
                continue;

            Hash = (Hash, (t.Hash, content[t])).GetHashCode();
        }
    }

    public bool IsOne => Content.Count == 0;

    public Monom<T> Mul(Monom<T> g)
    {
        var content = new Dictionary<T, int>();
        foreach (var kp in Indeterminates)
        {
            var m = g[kp] + this[kp];
            if (m != 0)
                content[kp] = m;
        }

        return new(Indeterminates, content);
    }

    public Monom<T> Pow(int k)
    {
        if (k < 0) throw new Exception();

        var content = new Dictionary<T, int>();
        foreach (var kp in Indeterminates)
        {
            var m = k * this[kp];
            if (m != 0)
                content[kp] = m;
        }

        return new(Indeterminates, content);
    }

    public (bool, Monom<T>) Div(Monom<T> g)
    {
        var content = new Dictionary<T, int>();
        foreach (var t in Indeterminates)
        {
            var m = this[t] - g[t];
            if (m < 0)
                return (false, new(Indeterminates));

            if (m > 0)
                content[t] = m;
        }

        return (true, new(Indeterminates, content));
    }

    public (int n, Monom<T> m) D(T t)
    {
        if (Content.TryGetValue(t, out var n))
        {
            var (b, m) = Div(new(Indeterminates, t));
            if (b)
                return (n, m);
        }

        return (0, new(Indeterminates));
    }

    public (int n, Monom<T> m) Remove(T t)
    {
        var content = new Dictionary<T, int>(Content);
        if (content.ContainsKey(t))
        {
            var i = content[t];
            content.Remove(t);
            return (i, new(Indeterminates, content));
        }

        return (0, new(Indeterminates, content));
    }

    public bool Equals(Monom<T> other) => Hash == other.Hash && Content.Count == other.Content.Count &&
                                          Content.All(e => e.Value.Equals(other.Content[e.Key]));

    public int this[T t] => Content.ContainsKey(t) ? Content[t] : 0;

    public IEnumerable<(T, int)> ToTuples()
    {
        foreach (var c in Indeterminates)
            yield return (c, this[c]);
    }

    public int CompareTo(Monom<T> other)
    {
        if (Indeterminates.Graded)
        {
            var comp = Degree.CompareTo(other.Degree);
            if (comp != 0)
                return comp;
        }

        var sgn = Indeterminates.Reverse ? -1 : 1;
        return sgn * ToTuples().SequenceCompareTo(other.ToTuples());
    }

    public int Hash { get; }

    public override int GetHashCode() => Hash;

    public override string ToString()
    {
        var sep = (Ring.DisplayPolynomial & MonomDisplay.Star) == MonomDisplay.Star ? "*" :
            (Ring.DisplayPolynomial & MonomDisplay.Dot) == MonomDisplay.Dot ? "." : "";

        if ((Ring.DisplayPolynomial & MonomDisplay.Caret) == MonomDisplay.Caret)
            return Content.OrderBy(e => e.Key).Select(kp => kp.Value == 1 ? $"{kp.Key}" : $"{kp.Key}^{kp.Value}")
                .Glue(sep);

        if ((Ring.DisplayPolynomial & MonomDisplay.PowFct) == MonomDisplay.PowFct)
            return Content.OrderBy(e => e.Key).Select(kp => kp.Value == 1 ? $"{kp.Key}" : $"{kp.Key}.Pow({kp.Value})")
                .Glue(sep);

        string Sup(int v)
        {
            var v0 = $"{v}";
            for (int i = 0; i < 10; i++)
                v0 = v0.Replace($"{i}", $"{superscripts[i]}");
            return v0;
        }

        var s = Content.OrderBy(e => e.Key)
            .Select(kp => kp.Value == 1 ? $"{kp.Key}" : $"{kp.Key}{Sup(kp.Value)}")
            .Glue(sep);

        return s;
    }
    
    private static string superscripts = "⁰¹²³⁴⁵⁶⁷⁸⁹";
    
    public static Monom<T> Gcd(Monom<T> a, Monom<T> b)
    {
        if (!a.Indeterminates.Equals(b.Indeterminates))
            throw new ArgumentException();

        var content = a.Indeterminates.Select(e => (e, k: Int32.Min(a[e], b[e]))).Where(e => e.k != 0)
            .ToDictionary(e => e.e, e => e.k);

        return new(a.Indeterminates, content);
    }

    public static (Monom<T> pa, Monom<T>pb) Reduce(Monom<T> a, Monom<T> b)
    {
        if (!a.Indeterminates.Equals(b.Indeterminates))
            throw new ArgumentException();

        var da = new Dictionary<T, int>();
        var db = new Dictionary<T, int>();
        foreach (var t in a.Indeterminates)
        {
            var at = a[t];
            var bt = b[t];
            var mx = Int32.Max(at, bt);
            if (mx == 0 || at == bt)
                continue;

            if (mx == bt)
                da[t] = mx - at;
            else
                db[t] = mx - bt;
        }

        return (new(a.Indeterminates, da), new(a.Indeterminates, db));
    }
}

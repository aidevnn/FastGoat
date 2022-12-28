using FastGoat.Commons;

namespace FastGoat.Structures.VecSpace;

[Flags]
public enum MonomOrder
{
    Lex = 0b00,
    Graded = 0b01,
    Reverse = 0b10,
    GrLex = Graded | Lex,
    RevLex = Reverse | Lex,
    GrevLex = Graded | Reverse
}

public readonly struct Monom<T> : IElt<Monom<T>> where T : struct, IElt<T>
{
    private Dictionary<T, int> Content { get; }
    public Indeterminates<T> Indeterminates { get; }
    public IEnumerable<T> ContentIndeterminates => Content.Keys;
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

        if (n < 1)
            throw new GroupException(GroupExceptionType.BaseGroup);

        Indeterminates = indeterminates;
        Content = new() { [s] = n };
        Degree = n;
        Hash = (Indeterminates.Hash, (s.Hash, n)).GetHashCode();
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
        Hash = content.Aggregate(Indeterminates.Hash, (acc, a) => (acc, (a.Key, a.Value)).GetHashCode());
    }

    public Monom<T> Mul(Monom<T> g)
    {
        var content = new Dictionary<T, int>(Content);
        foreach (var kp in g.Content)
        {
            if (content.ContainsKey(kp.Key))
                content[kp.Key] += kp.Value;
            else
                content[kp.Key] = kp.Value;
        }

        return new(Indeterminates, content);
    }

    public Monom<T>? Div(Monom<T> g)
    {
        var content = new Dictionary<T, int>(Content);
        foreach (var kp in g.Content)
        {
            if (content.ContainsKey(kp.Key))
            {
                var e = content[kp.Key] -= kp.Value;
                if (e < 0)
                    return null;

                if (e == 0)
                    content.Remove(kp.Key);
            }
            else
            {
                return null;
            }
        }

        return new(Indeterminates, content);
    }

    public (int n, Monom<T> m) D(T t)
    {
        if (Content.ContainsKey(t))
        {
            var m = Div(new(Indeterminates, t));
            if (m.HasValue)
            {
                var n = Content[t];
                return (n, m.Value);
            }
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

    public int DegreeOf(T t) => Content.ContainsKey(t) ? Content[t] : 0;

    public bool Equals(Monom<T> other) => Hash == other.Hash;

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
        var sep = (Monom.Display & MonomDisplay.Star) == MonomDisplay.Star ? "*" :
            (Monom.Display & MonomDisplay.Dot) == MonomDisplay.Dot ? "." : "";
        
        if ((Monom.Display & MonomDisplay.Caret) == MonomDisplay.Caret)
            return Content.OrderBy(e => e.Key).Select(kp => kp.Value == 1 ? $"{kp.Key}" : $"{kp.Key}^{kp.Value}")
                .Glue(sep);

        if ((Monom.Display & MonomDisplay.PowFct) == MonomDisplay.PowFct)
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
}

[Flags]
public enum MonomDisplay
{
    Caret = 1,
    Superscript = 2,
    PowFct = 4,
    Star = 8,
    Dot = 16,
    Default = Superscript,
    StarSuperscript = Star | Superscript,
    StarCaret = Star | Caret,
    StarPowFct = Star | PowFct,
    DotSuperscript = Dot | Superscript,
    DotCaret = Dot | Caret,
    DotPowFct = Dot | PowFct
}

public static class Monom
{
    public static MonomDisplay Display = MonomDisplay.Default;
}
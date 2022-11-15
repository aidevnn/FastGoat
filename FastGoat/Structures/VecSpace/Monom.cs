using FastGoat.Commons;

namespace FastGoat.Structures.VecSpace;

public readonly struct Monom<T> : IElt<Monom<T>> where T : struct, IElt<T>
{
    private SortedList<T, int> Content { get; }
    public int Degree { get; }
    public IEnumerable<T> Indeterminates => Content.Keys;

    public Monom()
    {
        Content = new();
        Degree = 0;
        Hash = 0;
    }

    public Monom(T s)
    {
        Content = new() { [s] = 1 };
        Degree = 1;
        Hash = (0, (s.Hash, 1)).GetHashCode();
    }

    public Monom(T s, int n)
    {
        if (n < 1)
            throw new GroupException(GroupExceptionType.BaseGroup);

        Content = new() { [s] = n };
        Degree = n;
        Hash = (0, (s.Hash, n)).GetHashCode();
    }

    private Monom(SortedList<T, int> content)
    {
        Content = content;
        Degree = content.Sum(e => e.Value);
        Hash = content.Aggregate(0, (acc, a) => (acc, (a.Key, a.Value)).GetHashCode());
    }

    public Monom<T> Mul(Monom<T> g)
    {
        var content = new SortedList<T, int>(Content);
        foreach (var kp in g.Content)
        {
            if (content.ContainsKey(kp.Key))
                content[kp.Key] += kp.Value;
            else
                content[kp.Key] = kp.Value;
        }

        return new(content);
    }

    public Monom<T>? Div(Monom<T> g)
    {
        var content = new SortedList<T, int>(Content);
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

        return new(content);
    }

    public (int n, Monom<T> m) D(T t)
    {
        if (Content.ContainsKey(t))
        {
            var m = Div(new(t));
            if (m.HasValue)
            {
                var n = Content[t];
                return (n, m.Value);
            }
        }

        return (0, new());
    }

    public (int n, Monom<T> m) Remove(T t)
    {
        var content = new SortedList<T, int>(Content);
        if (content.ContainsKey(t))
        {
            var i = content[t];
            content.Remove(t);
            return (i, new(content));
        }

        return (0, new(content));
    }

    public int DegreeOf(T t) => Content.ContainsKey(t) ? Content[t] : 0;

    public bool Equals(Monom<T> other) => Hash == other.Hash;

    public int CompareTo(Monom<T> other)
    {
        var compD = Degree.CompareTo(other.Degree);
        if (compD != 0)
            return -compD;

        return Content.Select(kp => (kp.Key, -kp.Value))
            .SequenceCompareTo(other.Content.Select(kp => (kp.Key, -kp.Value)));
    }

    public int Hash { get; }

    public override int GetHashCode() => Hash;

    public override string ToString()
    {
        var sep = (Monom.Display & MonomDisplay.Star) == MonomDisplay.Star ? "*" : "";
        if ((Monom.Display & MonomDisplay.Caret) == MonomDisplay.Caret)
            return Content.Select(kp => kp.Value == 1 ? $"{kp.Key}" : $"{kp.Key}^{kp.Value}").Glue(sep);
        
        if ((Monom.Display & MonomDisplay.PowFct) == MonomDisplay.PowFct)
            return Content.Select(kp => kp.Value == 1 ? $"{kp.Key}" : $"{kp.Key}.Pow({kp.Value})").Glue(sep);
        
        var s = Content.Select(kp => kp.Value == 1 ? $"{kp.Key}" : $"{kp.Key}{kp.Value}").Glue(sep);
        for (int i = 0; i < 10; i++)
            s = s.Replace($"{i}", $"{superscripts[i]}");

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
    Default = Superscript,
    StarSuperscript = Star|Superscript,
    StarCaret=Star|Caret,
    StarPowFct = Star|PowFct
}
public static class Monom
{
    public static MonomDisplay Display = MonomDisplay.Default;
}
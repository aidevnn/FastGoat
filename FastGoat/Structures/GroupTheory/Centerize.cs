using FastGoat.Structures.SetTheory;

namespace FastGoat.Structures.GroupTheory;

public class Centerize<U> : SubGroup<U> where U : struct, IElt<U>
{
    public Centerize(SubGroup<U> h, SubSet<U> s) : base(h.UpperGroup)
    {
        if (!h.UpperGroup.Equals(s.UpperSet))
            return;

        Generate(h.AllElements.ToHashSet(new EltEquality<U>()), s.AllElements.ToHashSet(new EltEquality<U>()));
        SetName($"Z[{h.Infos.Name}]({s.Infos.Name})");
    }

    public override U Neutral => UpperGroup.Neutral;
    public override U Invert(U a) => UpperGroup.Invert(a);
    public override U Op(U a, U b) => UpperGroup.Op(a, b);

    void Generate(HashSet<U> h, HashSet<U> s)
    {
        foreach (var x in h)
        {
            if (s.All(e1 => e1.Equals(Op(Op(x, e1), Invert(x)))))
                AddElement(x);
        }
    }
}

public static partial class GroupExt
{
    public static Centerize<U> CenterizeIn<U>(this SubSet<U> s, SubGroup<U> h) where U : struct, IElt<U>
    {
        return new Centerize<U>(h, s);
    }

    public static Centerize<U> Centerize<U>(this SubGroup<U> h, SubSet<U> s) where U : struct, IElt<U>
    {
        return new Centerize<U>(h, s);
    }
}
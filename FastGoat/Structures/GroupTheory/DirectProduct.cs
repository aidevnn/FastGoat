using FastGoat.Structures.SetTheory;

namespace FastGoat.Structures.GroupTheory;


public class DirectProduct<U> : SubGroup<U> where U : struct, IElt<U>
{
    public DirectProduct(SubGroup<U> g, SubGroup<U> h) : base(g.UpperGroup)
    {
        if (!g.UpperGroup.Equals(h.UpperGroup))
            return;

        Generate(g, h);
        SetName($"{g.Infos.Name}.{h.Infos.Name}");
    }

    public override U Neutral => UpperGroup.Neutral;
    public override U Invert(U a) => UpperGroup.Invert(a);
    public override U Op(U a, U b) => UpperGroup.Op(a, b);

    void Generate(SubGroup<U> g, SubGroup<U> h)
    {
        foreach (var e0 in g.AllElements())
        {
            foreach (var e1 in h.AllElements())
            {
                var e2 = Op(e0, e1);
                AddElement(e2);
            }
        }
    }
}

public static partial class GroupExt
{
    public static SubGroup<U> DirectProduct<U>(this SubGroup<U> g, params SubGroup<U>[] subGroups) where U : struct, IElt<U>
    {
        var acc = new DirectProduct<U>(g, g.UpperGroup.Singleton());
        return subGroups.Aggregate(acc, (a, h) => new DirectProduct<U>(a, h));
    }

    public static SubGroup<U> Generate<U>(this SubGroup<U> g) where U : struct, IElt<U>
    {
        SubGroup<U> h = g.UpperGroup.GroupUnion(g, g.UpperGroup.Singleton());
        int sz = 0;
        do
        {
            sz = h.Count;
            h = h.DirectProduct(h);
        } while (sz != h.Count);

        return h;
    }

    public static SubGroup<U> Generate<U>(this SubSet<U> g, IGroup<U> group) where U : struct, IElt<U>
    {
        return group.GroupElement(g).Generate();
    }
}
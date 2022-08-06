using FastGoat.Structures.SetTheory;

namespace FastGoat.Structures.GroupTheory;

public class GroupElement<U> : SubGroup<U> where U : struct, IElt<U>
{
    public GroupElement(ISubGroup<U> group, IEnumerable<U> us) : base(group)
    {
        if (us.Any(e => !group.Ancestor.Equals(e.FSet)))
            return;

        foreach (var e in us)
            AddElement(Op(Neutral, e));

        SetName("H");
    }

    public override U Neutral => UpperGroup.Neutral;
    public override U Invert(U a) => UpperGroup.Invert(a);
    public override U Op(U a, U b) => UpperGroup.Op(a, b);
}

public static partial class GroupExt
{
    public static GroupElement<U> GroupElement<U>(this Group<U> group, params U[] us) where U : struct, IElt<U>
    {
        return us.Length == 0 ? new GroupElement<U>(group.Singleton(), group.AllElements) : new GroupElement<U>(group.Singleton(), us);
    }

    public static GroupElement<U> GroupElement<U>(this Group<U> group, IEnumerable<U> us) where U : struct, IElt<U>
    {
        return group.GroupElement(us.ToArray());
    }

    public static GroupElement<U> GroupElement<U>(this ISubGroup<U> group, params U[] us) where U : struct, IElt<U>
    {
        return new GroupElement<U>(group, us);
    }

    public static GroupElement<U> GroupUnion<U>(this ISubGroup<U> group, params SubSet<U>[] subSets) where U : struct, IElt<U>
    {
        var Ugi = new GroupElement<U>(group, group.AllElements.Union(subSets.SelectMany(h => h.AllElements)));
        Ugi.Infos.Name = subSets.Select(h => h.Infos.Name).Glue("{0}", "u");
        return Ugi;
    }

}
using FastGoat.Structures.SetTheory;

namespace FastGoat.Structures.GroupTheory;


public class GroupElement<U> : SubGroup<U> where U : struct, IElt<U>
{
    public GroupElement(IGroup<U> group, IEnumerable<U> us) : base(group)
    {
        if (us.Any(e => !group.Equals(e.FSet)))
            return;

        foreach (var e in us)
            AddElement(e);

        SetName("H");
    }

    public override U Neutral => UpperGroup.Neutral;
    public override U Invert(U a) => UpperGroup.Invert(a);
    public override U Op(U a, U b) => UpperGroup.Op(a, b);
}

public static partial class GroupExt
{
    public static GroupElement<U> GroupElement<U>(this IGroup<U> group, params U[] us) where U : struct, IElt<U>
    {
        return new GroupElement<U>(group, us);
    }

    public static GroupElement<U> GroupElement<U>(this IGroup<U> group, SubSet<U> subSet) where U : struct, IElt<U>
    {
        return new GroupElement<U>(group, subSet.AllElements());
    }

    public static GroupElement<U> Singleton<U>(this IGroup<U> group) where U : struct, IElt<U>
    {
        return group.GroupElement(group.Neutral);
    }

    public static GroupElement<U> GroupUnion<U>(this IGroup<U> group, params SubSet<U>[] subSets) where U : struct, IElt<U>
    {
        var Ugi = new GroupElement<U>(group, subSets.SelectMany(h => h.AllElements()).ToArray());
        Ugi.Infos.Name = subSets.Select(h => h.Infos.Name).Glue("{0}", "u");
        return Ugi;
    }

}
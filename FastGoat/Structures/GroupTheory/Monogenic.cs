using FastGoat.Structures.SetTheory;

namespace FastGoat.Structures.GroupTheory;

public class Monogenic<U> : SubGroup<U> where U : struct, IElt<U>
{
    public Monogenic(IGroup<U> group, U e) : base(group)
    {
        if (group.Contains(e))
            Generate(e);

        SetName($"<{e}>".Replace(" ", ""));
    }

    public override U Neutral => UpperGroup.Neutral;
    public override U Invert(U a) => UpperGroup.Invert(a);
    public override U Op(U a, U b) => UpperGroup.Op(a, b);

    void Generate(U e)
    {
        AddElement(UpperGroup.Neutral);
        var acc = e;
        while (!acc.Equals(Neutral))
        {
            AddElement(acc);
            acc = Op(e, acc);
        }
    }
}

public static partial class GroupExt
{
    public static SubGroup<U> Monogenic<U>(this IGroup<U> group, U e) where U : struct, IElt<U>
    {
        return new Monogenic<U>(group, e);
    }
}
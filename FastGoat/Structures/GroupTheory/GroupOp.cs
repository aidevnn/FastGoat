using FastGoat.Structures.SetTheory;

namespace FastGoat.Structures.GroupTheory;

public class GroupOp<U> : SubGroup<U> where U : struct, IElt<U>
{
    public GroupOp(SubGroup<U> sub, U e) : base(sub.UpperGroup)
    {
        if (!sub.UpperSet.Equals(e.FSet))
            return;

        Generate(sub.AllElements, e);
        SetName("Hx");
    }

    public GroupOp(U e, SubGroup<U> sub) : base(sub.UpperGroup)
    {
        if (!sub.UpperSet.Equals(e.FSet))
            return;

        Generate(e, sub.AllElements);
        SetName("xH");
    }

    void Generate(U e, IEnumerable<U> sub)
    {
        foreach (var e0 in sub)
            AddElement(Op(e, e0));
    }

    void Generate(IEnumerable<U> sub, U e)
    {
        foreach (var e0 in sub)
            AddElement(Op(e0, e));
    }

    public override U Neutral => UpperGroup.Neutral;
    public override U Invert(U a) => UpperGroup.Invert(a);
    public override U Op(U a, U b) => UpperGroup.Op(a, b);
}


public static partial class GroupExt
{
    public static GroupOp<U> Gop<U>(this SubGroup<U> h, U x) where U : struct, IElt<U>
    {
        return new GroupOp<U>(h, x);
    }

    public static GroupOp<U> Gop<U>(this U x, SubGroup<U> h) where U : struct, IElt<U>
    {
        return new GroupOp<U>(x, h);
    }
}
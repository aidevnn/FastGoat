using FastGoat.Structures.SetTheory;

namespace FastGoat.Structures.GroupTheory;

public class Singleton<U> : SubGroup<U> where U : struct, IElt<U>
{
    public Singleton(IGroup<U> group) : base(group)
    {
        AddElement(group.Neutral);
        SetName("H");
    }

    public override void AddElement(U e)
    {
        base.AddElement(e);
    }

    public override U Neutral => UpperGroup.Neutral;
    public override U Invert(U a) => UpperGroup.Invert(a);
    public override U Op(U a, U b) => UpperGroup.Op(a, b);
}

public static partial class GroupExt
{
    public static Singleton<U> Singleton<U>(this IGroup<U> group) where U : struct, IElt<U>
    {
        return new Singleton<U>(group);
    }

}
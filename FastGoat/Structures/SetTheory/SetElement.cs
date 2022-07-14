namespace FastGoat.Structures.SetTheory;

public class SetElement<U> : SubSet<U> where U : struct, IElt<U>
{
    public SetElement(IFSet<U> fSet, params U[] us) : base(fSet)
    {
        if (us.Any(e => !UpperSet.Equals(e.FSet)))
            return;

        foreach (var e in us)
            AddElement(e);
    }
}

public static partial class SetExt
{
    public static SubSet<U> SnapShot<U>(this FSet<U> fSet) where U : struct, IElt<U>
    {
        return new SetElement<U>(fSet, fSet.AllElements.ToArray());
    }
}
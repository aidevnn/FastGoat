namespace FastGoat.Structures.SetTheory;
public class EmptySet<U> : SubSet<U> where U : struct, IElt<U>
{
    public EmptySet(IFSet<U> fSet) : base(fSet) { }
    public override void AddElement(U e) { }
}

public static partial class SetExt
{
    public static SubSet<U> EmptySet<U>(this IFSet<U> fSet) where U : struct, IElt<U>
    {
        return new EmptySet<U>(fSet);
    }
}
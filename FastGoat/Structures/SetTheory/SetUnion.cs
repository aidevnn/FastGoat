namespace FastGoat.Structures.SetTheory;

public class SetUnion<U> : SetElement<U> where U : struct, IElt<U>
{
    public SetUnion(SubSet<U> g, SubSet<U> h) : base(g.UpperSet, g.AllElements.Union(h.AllElements).ToArray())
    {

    }

    public SetUnion(SubSet<U> g, params U[] us) : base(g.UpperSet, g.AllElements.Union(us).ToArray())
    {

    }
}

public static partial class SetExt
{
    public static SubSet<U> Union<U>(this SubSet<U> subSet, params U[] us) where U : struct, IElt<U>
    {
        return new SetUnion<U>(subSet, us);
    }

    public static SubSet<U> Union<U>(this SubSet<U> subSet, SubSet<U> h) where U : struct, IElt<U>
    {
        return new SetUnion<U>(subSet, h.AllElements.ToArray());
    }

    public static SubSet<U> Union<U>(this IFSet<U> fSet, params U[] us) where U : struct, IElt<U>
    {
        return new SetUnion<U>(fSet.EmptySet(), us);
    }
}
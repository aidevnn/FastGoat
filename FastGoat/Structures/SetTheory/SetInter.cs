namespace FastGoat.Structures.SetTheory;

public class SetInter<U> : SetElement<U> where U : struct, IElt<U>
{
    public SetInter(SubSet<U> g, SubSet<U> h) : base(g.UpperSet, h.AllElements().Intersect(g.AllElements()).ToArray())
    {

    }
}

public static partial class SetExt
{
    public static SubSet<U> SetInter<U>(this SubSet<U> subSet, SubSet<U> h) where U : struct, IElt<U>
    {
        return new SetInter<U>(subSet, h);
    }
}
namespace FastGoat.Structures.SetTheory;

public class SetExcept<U> : SetElement<U> where U : struct, IElt<U>
{
    public SetExcept(SubSet<U> g, SubSet<U> h) : base(g.UpperSet, g.AllElements().Except(h.AllElements()).ToArray())
    {

    }

    public SetExcept(SubSet<U> g, params U[] us) : base(g.UpperSet, g.AllElements().Except(us).ToArray())
    {

    }
}

public static partial class SetExt
{
    public static SubSet<U> Except<U>(this SubSet<U> subSet, params U[] us) where U : struct, IElt<U>
    {
        return new SetExcept<U>(subSet, us);
    }
}
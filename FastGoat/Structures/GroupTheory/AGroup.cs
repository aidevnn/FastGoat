using FastGoat.Structures.SetTheory;

namespace FastGoat.Structures.GroupTheory;

public abstract class Group<U> : FSet<U>, IGroup<U> where U : struct, IElt<U>
{
    public abstract U Neutral { get; }
    public abstract U Invert(U a);
    public abstract U Op(U a, U b);

}

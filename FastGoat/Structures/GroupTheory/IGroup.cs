using FastGoat.Structures.SetTheory;

namespace FastGoat.Structures.GroupTheory;

public interface IGroup<U> : IFSet<U> where U : struct, IElt<U>
{
    U Neutral { get; }
    U Invert(U a);
    U Op(U a, U b);
}

public interface ISubGroup<U> : IGroup<U> where U : struct, IElt<U>
{
    IGroup<U> UpperGroup { get; }
    IEnumerable<IGroup<U>> UpperGroupChain();
    IGroup<U> Ancestor { get; }
}
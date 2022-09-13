namespace FastGoat;

public interface IConcreteGroup<U> : IEnumerable<U>, IEquatable<IConcreteGroup<U>> where U : struct, IElt<U>
{
    IGroup<U> BaseGroup { get; }
    IConcreteGroup<U> ControlGroup { get; }
    int Hash { get; }
    U Neutral();
    U Invert(U a);
    U Op(U a, U b);
    U this[int k] { get; }
}

namespace FastGoat;

public interface IGroup<U> : IEquatable<IGroup<U>> where U : struct, IElt<U>
{
    int Hash { get; }
    U Neutral();
    U Invert(U a);
    U Op(U a, U b);
    U this[int k] { get; }
    string ToString();
}

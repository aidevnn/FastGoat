namespace FastGoat;

public interface IElt<U> : IEquatable<U>, IComparable<U> where U : struct, IElt<U>
{
    IGroup<U> Group { get; }
    int Hash { get; }
    int GetHashCode();
    string ToString();
}

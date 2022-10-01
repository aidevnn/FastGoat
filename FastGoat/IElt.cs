namespace FastGoat;

public interface IElt<T> : IEquatable<T>, IComparable<T> where T : IElt<T>
{
    int Hash { get; }
    IGroup<T> BaseGroup { get; }
}
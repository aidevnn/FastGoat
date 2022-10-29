namespace FastGoat.Theory;

public interface IGroup<T> : IEnumerable<T>, IEquatable<IGroup<T>> where T : IElt<T>
{
    int Hash { get; }
    string Name { get; }
    T this[params ValueType[] us] { get; }
    IEnumerable<T> GetElements();
    IEnumerable<T> GetGenerators();
    T Neutral();
    T Invert(T e);
    T Op(T e1, T e2);
}
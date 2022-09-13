namespace FastGoat;

public class EqualityHashSet<T> : EqualityComparer<HashSet<T>> where T : struct, IEquatable<T>
{
    public override bool Equals(HashSet<T>? x, HashSet<T>? y) => x is not null && y is not null && x.SetEquals(y);
    public override int GetHashCode(HashSet<T> obj) => obj.Count;
}
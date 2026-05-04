using System.Collections;
using FastGoat.Commons;
using FastGoat.Structures;

namespace Craft;

public readonly struct XSet<T> : IElt<XSet<T>>, IEnumerable<T> where T:struct,IEquatable<T>,IComparable<T>
{
    public HashSet<T> X { get; }
    public int Count => X.Count;

    public XSet()
    {
        X = new();
        Hash = X.Count;
    }
    
    public XSet(IEnumerable<T> set)
    {
        X = new(set);
        Hash = X.Count;
    }

    public XSet(params T[] set) : this(set.AsEnumerable())
    {
        
    }

    public bool Equals(XSet<T> other) => X.SetEquals(other.X);
    public bool Equals(HashSet<T> other) => X.SetEquals(other);

    public int CompareTo(XSet<T> other)
    {
        var comp = X.Count.CompareTo(other.X.Count);
        if (comp != 0)
            return comp;

        return X.Order().SequenceCompareTo(other.X.Order());
    }

    public bool Overlaps(IEnumerable<T> other) => X.Overlaps(other);
    public XSet<T> Append(T e) => new(X.Append(e));
    public XSet<T> Append(IEnumerable<T> e) => new(X.Concat(e));

    public int Hash { get; }
    public IEnumerator<T> GetEnumerator() => X.GetEnumerator();

    IEnumerator IEnumerable.GetEnumerator() => GetEnumerator();

    public override int GetHashCode() => Hash;
    public override string ToString() => $"{{ {X.Order().Glue(", ")} }}";
}

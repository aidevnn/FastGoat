using System.Collections;
using FastGoat.Commons;

namespace FastGoat.Structures.Subgroups;

public readonly struct XArray<T> : IElt<XArray<T>>, IEnumerable<T> where T : struct, IEquatable<T>, IComparable<T>
{
    public T[] X { get; }
    public int Count => X.Length;
    
    public XArray()
    {
        X = new T[0];
        Hash = 0;
    }

    public XArray(IEnumerable<T> set)
    {
        X = set.ToArray();
        Hash = X.Length;
    }

    public XArray(params T[] set) : this(set.AsEnumerable())
    {
    }

    public bool Equals(XArray<T> other) => X.SequenceEqual(other.X);

    public int CompareTo(XArray<T> other) => X.SequenceCompareTo(other.X);

    public int Hash { get; }
    public IEnumerator<T> GetEnumerator() => X.AsEnumerable().GetEnumerator();

    IEnumerator IEnumerable.GetEnumerator() => GetEnumerator();

    public T this[int idx] => X[idx];

    public override int GetHashCode() => Hash;
    public override string ToString() => $"[{X.Glue(", ")}]";
}


public static class XArrayExt
{
    extension<T>(IEnumerable<T> e) where T : struct, IEquatable<T>, IComparable<T>
    {
        public XArray<T> ToXArray() => new(e);
    }
}
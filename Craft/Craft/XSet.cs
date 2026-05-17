using System.Collections;
using FastGoat.Commons;
using FastGoat.Structures;

namespace Craft;

public readonly struct XSet<T> : IElt<XSet<T>>, IEnumerable<T> where T : struct, IEquatable<T>, IComparable<T>
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

    public int Hash { get; }
    public IEnumerator<T> GetEnumerator() => X.GetEnumerator();

    IEnumerator IEnumerable.GetEnumerator() => GetEnumerator();

    public override int GetHashCode() => Hash;
    public override string ToString() => $"{{ {X.Order().Glue(", ")} }}";

    public void ExceptWith(HashSet<T> orbx)
    {
        X.ExceptWith(orbx);
    }
}

public static class XSetExt
{
    extension<T>(IEnumerable<T> e) where T : struct, IEquatable<T>, IComparable<T>
    {
        public XSet<T> ToXSet() => new(e);
    }

    extension<T>(XSet<T> a) where T : struct, IEquatable<T>, IComparable<T>
    {
        public bool Overlaps(IEnumerable<T> other) => a.X.Overlaps(other);
        public XSet<T> Append(T e) => new(a.X.Append(e));
        public XSet<T> Concat(IEnumerable<T> e) => new(a.X.Concat(e));
    }

    extension<T>(IEnumerable<XSet<T>> ts) where T : struct, IElt<T>
    {
        public XSet<T> Intersect() => ts.Aggregate((e, f) => new(f.X.Intersect(e)));
        public XSet<T> Union() => ts.Aggregate((e, f) => new(e.X.Concat(f)));
    }
}
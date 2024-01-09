using System.Collections;

namespace FastGoat.Structures.GenericGroup;

public readonly struct GroupSubset<T>(HashSet<T> gens, HashSet<T> content) : IEnumerable<T>, IElt<GroupSubset<T>> where T : struct, IElt<T>
{
    public HashSet<T> Generators => gens;
    public HashSet<T> Elements => content;
    public int Count => content.Count;
    public bool Contains(T e) => content.Contains(e);
    public bool SuperSetOf(IEnumerable<T> other) => Elements.IsSupersetOf(other);
    public bool SubSetOf(IEnumerable<T> other) => Elements.IsSubsetOf(other);
    public bool SetEquals(IEnumerable<T> other) => Elements.SetEquals(other);
    public bool Equals(GroupSubset<T> other) => Elements.SetEquals(other.Elements);

    public int CompareTo(GroupSubset<T> other) => Count.CompareTo(other.Count);

    private int IdHash => Elements.Order().Aggregate(1, (acc, e) => (acc, e.GetHashCode()).GetHashCode());

    public int Hash => Count;
    public IEnumerator<T> GetEnumerator() => content.GetEnumerator();
    IEnumerator IEnumerable.GetEnumerator()
    {
        return GetEnumerator();
    }

    public override int GetHashCode() => Hash;
    public override string ToString() => $"Set[{Count}]{IdHash}";
}
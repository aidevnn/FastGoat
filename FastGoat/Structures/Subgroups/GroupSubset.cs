using System.Collections;
using FastGoat.Structures.GenericGroup;

namespace FastGoat.Structures.Subgroups;

public readonly struct GroupSubset<T>(HashSet<T> gens, HashSet<T> content)
    : IEnumerable<T>, IElt<GroupSubset<T>> where T : struct, IElt<T>
{
    public HashSet<T> Generators => gens;
    public HashSet<T> Elements => content;
    public int Count => content.Count;
    public bool Contains(T e) => content.Contains(e);
    public bool SuperSetOf(IEnumerable<T> other) => content.IsSupersetOf(other);
    public bool SubSetOf(IEnumerable<T> other) => content.IsSubsetOf(other);
    public bool SetEquals(IEnumerable<T> other) => content.SetEquals(other);
    public bool Equals(GroupSubset<T> other) => content.Count == other.Count && content.SetEquals(other.Elements);

    public int CompareTo(GroupSubset<T> other) => Count.CompareTo(other.Count);

    private int IdHash => content.Order().Aggregate(1, (acc, e) => (acc, e.GetHashCode()).GetHashCode());

    public int Hash => Count;
    public IEnumerator<T> GetEnumerator() => content.GetEnumerator();

    IEnumerator IEnumerable.GetEnumerator()
    {
        return GetEnumerator();
    }

    public GroupSubset<WElt> ToWElt() =>
        new(Generators.Select(e => new WElt(e)).ToHashSet(), Elements.Select(e => new WElt(e)).ToHashSet());

    public override int GetHashCode() => Hash;
    public override string ToString() => $"Set[{Count}]{IdHash}";
}
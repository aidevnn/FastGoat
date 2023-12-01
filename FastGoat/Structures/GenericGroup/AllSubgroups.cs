using System.Collections;

namespace FastGoat.Structures.GenericGroup;

public readonly struct AllSubgroups<T> : IEnumerable<SubgroupConjugates<T>> where T : struct, IElt<T>
{
    public SubgroupConjugates<T>[] AllSubgroupConjugates { get; }
    public ConcreteGroup<T> Parent => AllSubgroupConjugates.First().Parent;

    public AllSubgroups(ConcreteGroup<T> g)
    {
        var allSubs = Group.AllSubGroups(g);
        AllSubgroupConjugates = allSubs.Values.Select(l => new SubgroupConjugates<T>(g, l)).Order().ToArray();
    }

    public AllSubgroups(Dictionary<ConcreteGroup<T>, List<ConcreteGroup<T>>> allSubs)
    {
        var g = allSubs.MaxBy(sg => sg.Key.Count()).Key;
        AllSubgroupConjugates = allSubs.Values.Select(l => new SubgroupConjugates<T>(g, l)).Order().ToArray();
    }

    private AllSubgroups(HashSet<SubgroupConjugates<T>> all)
    {
        AllSubgroupConjugates = all.Order().ToArray();
    }

    public AllSubgroups<T> Restriction(ConcreteGroup<T> g) => new(AllSubgroupConjugates.Select(sc => sc.Restriction(g)).ToHashSet());

    public IEnumerator<SubgroupConjugates<T>> GetEnumerator() => AllSubgroupConjugates.AsEnumerable().GetEnumerator();

    IEnumerator IEnumerable.GetEnumerator()
    {
        return GetEnumerator();
    }

    public string Name => $"SubGroups of {Parent}";

    public IEnumerable<ConcreteGroup<T>> All => AllSubgroupConjugates.SelectMany(sc => sc.Conjugates);
    public IEnumerable<ConcreteGroup<T>> AllRepresentatives => AllSubgroupConjugates.Select(sc => sc.Representative);
    public override string ToString() => Name;
}
using System.Collections;
using FastGoat.Commons;

namespace FastGoat.Structures.GenericGroup;

public record SubGroupsInfos(int AllSubGr, int AllConjsCl, int AllNorms)
{
    public static explicit operator (int, int, int)(SubGroupsInfos infos) => (infos.AllSubGr, infos.AllConjsCl, infos.AllNorms);
}

public readonly struct AllSubgroups<T> : IEnumerable<SubgroupConjugates<T>> where T : struct, IElt<T>
{
    public SubgroupConjugates<T>[] AllSubgroupConjugates { get; }
    public ConcreteGroup<T> Parent => AllSubgroupConjugates.First().Parent;

    public AllSubgroups(ConcreteGroup<T> g)
    {
        var allSubs = Group.AllSubGroups(g);
        Infos = new(allSubs.Values.Sum(s => s.Count), allSubs.Count, allSubs.Count(s => s.Value.Count == 1));
        AllSubgroupConjugates = allSubs.Values.Select(l => new SubgroupConjugates<T>(g, l)).Order().ToArray();
    }

    public AllSubgroups(Dictionary<ConcreteGroup<T>, List<ConcreteGroup<T>>> allSubs)
    {
        var g = allSubs.MaxBy(sg => sg.Key.Count()).Key;
        Infos = new(allSubs.Values.Sum(s => s.Count), allSubs.Count, allSubs.Count(s => s.Value.Count == 1));
        AllSubgroupConjugates = allSubs.Values.Select(l => new SubgroupConjugates<T>(g, l)).Order().ToArray();
    }

    private AllSubgroups(HashSet<SubgroupConjugates<T>> all)
    {
        AllSubgroupConjugates = all.Order().ToArray();
        var conjs = AllSubgroupConjugates.Length;
        var subs = AllSubgroupConjugates.Sum(sc => sc.Size);
        var norms = AllSubgroupConjugates.Count(sc => sc.IsNormal);
        Infos = new(subs, conjs, norms);
    }

    public bool IsSimple()
    {
        if (Parent.GroupType == GroupType.AbelianGroup && IntExt.Primes10000.Contains(Parent.Count()))
            return true;

        return AllSubgroupConjugates.Count(sc => sc.IsNormal) == 2;
    }

    public AllSubgroups<T> Restriction(ConcreteGroup<T> g) => new(AllSubgroupConjugates.Select(sc => sc.Restriction(g)).ToHashSet());

    public IEnumerator<SubgroupConjugates<T>> GetEnumerator() => AllSubgroupConjugates.AsEnumerable().GetEnumerator();

    IEnumerator IEnumerable.GetEnumerator()
    {
        return GetEnumerator();
    }

    public SubGroupsInfos Infos { get; }

    public string Name => $"SubGroups of {Parent}";

    public IEnumerable<ConcreteGroup<T>> All => AllSubgroupConjugates.SelectMany(sc => sc.Conjugates);
    public IEnumerable<ConcreteGroup<T>> AllRepresentatives => AllSubgroupConjugates.Select(sc => sc.Representative);
    public override string ToString() => Name;
}
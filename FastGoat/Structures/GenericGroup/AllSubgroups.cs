using System.Collections;
using FastGoat.Commons;

namespace FastGoat.Structures.GenericGroup;

public record SubGroupsInfos(int AllSubGr, int AllConjsCl, int AllNorms) : IComparable<SubGroupsInfos>
{
    public (int, int, int) ToTuples() => (AllSubGr, AllConjsCl, AllNorms);
    public int CompareTo(SubGroupsInfos? other)
    {
        if (other is null)
            return 1;

        return ToTuples().CompareTo(other.ToTuples());
    }
}

public readonly struct AllSubgroups<T> : IEnumerable<SubgroupConjugates<T>>, IEquatable<AllSubgroups<T>> where T : struct, IElt<T>
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

    private AllSubgroups(SubgroupConjugates<T>[] all)
    {
        AllSubgroupConjugates = all.Order().ToArray();
        var conjs = AllSubgroupConjugates.Length;
        var subs = AllSubgroupConjugates.Sum(sc => sc.Size);
        var norms = AllSubgroupConjugates.Count(sc => sc.IsNormal);
        Infos = new(subs, conjs, norms);
    }

    public bool Equals(AllSubgroups<T> other) => Parent.SetEquals(other.Parent);
    public override int GetHashCode() => Parent.Hash;

    public bool IsSimple()
    {
        return Infos.AllNorms == 2;
    }

    public AllSubgroups<T> Restriction(ConcreteGroup<T> g)
    {
        return new(AllSubgroupConjugates
            .Where(sc => sc.Conjugates.Any(e => e.SubSetOf(g)))
            .SelectMany(sc => sc.Restriction(g))
            .Distinct()
            .ToArray());
    }

    public IEnumerator<SubgroupConjugates<T>> GetEnumerator() => AllSubgroupConjugates.AsEnumerable().GetEnumerator();

    IEnumerator IEnumerable.GetEnumerator()
    {
        return GetEnumerator();
    }

    public SubGroupsInfos Infos { get; }

    public Dictionary<int, SubgroupConjugates<T>[]> AllSylows
    {
        get
        {
            var ord = Parent.Count();
            var allConjs = AllSubgroupConjugates;
            return IntExt.PrimesDec(ord).Select(e => e.Key.Pow(e.Value))
                .ToDictionary(q => q, q => allConjs.Where(cj => cj.Order == q).ToArray());
        }
    }

    public string Name => $"SubGroups of {Parent}";

    public AllSubgroups<WElt> ToGroupWrapper()
    {
        return new AllSubgroups<WElt>(AllSubgroupConjugates.Select(sc => sc.ToGroupWrapper()).ToArray());
    }

    public SubgroupConjugates<T>[] MaximalSubgroups()
    {
        var allMax = new List<SubgroupConjugates<T>>();
        foreach (var h in AllSubgroupConjugates.Where(cj => cj.IsProper))
        {
            allMax = allMax.Except(allMax.Where(h0 => h0.IsSubClassOf(h))).ToList();
            if (allMax.All(h0 => !h.IsSubClassOf(h0)))
                allMax.Add(h);
        }

        return allMax.Order().ToArray();
    }

    public ConcreteGroup<T> FrattiniSubGroup
    {
        get
        {
            var max = MaximalSubgroups().SelectMany(cj => cj.Conjugates)
                .Aggregate(Parent.AsEnumerable(), (acc, sg) => acc.Intersect(sg));
            return AllRepresentatives.First(sg => sg.SetEquals(max));
        }
    }

    public ConcreteGroup<T> FittingSubGroup
    {
        get
        {
            if (Parent.Count() == 1)
                return AllSubgroupConjugates.First(cj => cj.IsTrivial).Representative;

            return AllSubgroupConjugates.Where(cj => cj.IsProperNormal)
                .Descending()
                .Select(cj => cj.Representative)
                .First(sg => Group.IsNilpotent(sg));
        }
    }

    public IEnumerable<ConcreteGroup<T>> All => AllSubgroupConjugates.SelectMany(sc => sc.Conjugates);
    public IEnumerable<ConcreteGroup<T>> AllRepresentatives => AllSubgroupConjugates.Select(sc => sc.Representative);
    public override string ToString() => Name;
}

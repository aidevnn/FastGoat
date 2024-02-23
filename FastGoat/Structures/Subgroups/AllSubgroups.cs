using System.Collections;
using FastGoat.Commons;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.Naming;

namespace FastGoat.Structures.Subgroups;

public readonly struct AllSubgroups<T> : IEnumerable<SubgroupConjugates<T>>, IEquatable<AllSubgroups<T>> where T : struct, IElt<T>
{
    public SubgroupConjugates<T>[] AllSubgroupConjugates { get; }
    public ConcreteGroup<T> Parent => AllSubgroupConjugates.First().Parent;

    public AllSubgroups(ConcreteGroup<T> g)
    {
        var allSubs = Group.AllSubGroups(g);
        Infos = new(allSubs.Values.Sum(s => s.Count), allSubs.Count, allSubs.Count(s => s.Value.Count == 1));
        AllSubgroupConjugates = allSubs.Values.Select(l => new SubgroupConjugates<T>(g, l)).Order().ToArray();
        AllSubgroupConjugates.First().Conjugates[0].Name = "C1";
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

    public bool IsSylow()
    {
        return Factors.Count == 1;
    }

    public AllSubgroups<T> Restriction(ConcreteGroup<T> g)
    {
        if (g.Count() == Parent.Count())
            return this;

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
    public Dictionary<int, int> Factors => AllSubgroupConjugates.Last().Factors;
    
    public Dictionary<int, SubgroupConjugates<T>[]> AllSylows()
    {
        var allConjs = AllSubgroupConjugates;
        return Factors.Select(e => e.Key.Pow(e.Value))
            .ToDictionary(q => q, q => allConjs.Where(cj => cj.Order == q).ToArray());
    }

    public List<(int p, int q, SubgroupConjugates<T>)> AllPGroups()
    {
        return AllSubgroupConjugates.Where(cj => cj.IsPGroup())
            .Select(cj => (cj.Factors.First().Key, cj.Order, cj))
            .ToList();
    }

    public string Name => $"SubGroups of {Parent}";

    public AllSubgroups<WElt> ToGroupWrapper()
    {
        if (this is AllSubgroups<WElt> subgroups)
            return subgroups;

        return new AllSubgroups<WElt>(AllSubgroupConjugates.Select(sc => sc.ToGroupWrapper()).ToArray());
    }

    public SubgroupConjugates<T>[] MaximalSubgroups() => MaximalSubgroups(AllSubgroupConjugates.Last());

    public SubgroupConjugates<T>[] MaximalSubgroups(SubgroupConjugates<T> top)
    {
        var allMax = new List<SubgroupConjugates<T>>();
        foreach (var h in AllSubgroupConjugates.Where(cj => cj.IsSubClassOf(top)))
        {
            allMax = allMax.Except(allMax.Where(h0 => h0.IsSubClassOf(h))).ToList();
            if (allMax.All(h0 => !h.IsSubClassOf(h0)))
                allMax.Add(h);
        }

        return allMax.Order().ToArray();
    }

    public SubgroupConjugates<T>[] MaximalNormalSubgroups(SubgroupConjugates<T> top)
    {
        var allMax = new List<SubgroupConjugates<T>>();
        foreach (var h in AllSubgroupConjugates.Where(cj => cj.IsNormal && cj.IsSubClassOf(top)))
        {
            allMax = allMax.Except(allMax.Where(h0 => h0.IsSubClassOf(h))).ToList();
            if (allMax.All(h0 => !h.IsSubClassOf(h0)))
                allMax.Add(h);
        }

        return allMax.Order().ToArray();
    }

    public void Naming()
    {
        var subscript = "₀₁₂₃₄₅₆₇₈₉";
        if (All.All(sg => !sg.Name.Contains("SubGr")))
            return;

        if (this is AllSubgroups<WElt> subgroups)
        {
            Parent.Name = NamesTree.BuildName(subgroups)[0].Name;
            foreach (var cj in subgroups.AllSubgroupConjugates)
            {
                if (cj.Conjugates.All(sg => !sg.Name.Contains("SubGr")))
                    continue;

                var r = subgroups.Restriction(cj.Representative);
                var name = NamesTree.BuildName(r)[0].Name;
                cj.Conjugates.ForEach(sg => sg.Name = name);
            }

            foreach (var gr in AllSubgroupConjugates.Select((cj, i) => (cj.Representative.Name, i)).GroupBy(cj => cj.Name).ToArray())
            {
                if (gr.Count() == 1)
                    continue;

                var arrCjs = gr.ToArray();
                for (int k = 0; k < arrCjs.Length; ++k)
                {
                    var strK = $"{k + 1}";
                    for (int i = 0; i < 10; ++i)
                        strK = strK.Replace($"{i}", $"{subscript[i]}");

                    AllSubgroupConjugates[arrCjs[k].i].Subscript = strK;
                }
            }
        }
        else
            throw new("Naming avalaible only for GroupWrapper");
    }

    public List<SubgroupConjugates<T>>[] MaximalSubgroupSeries()
    {
        var g = AllSubgroupConjugates.Last();
        var all = MaximalSubgroups().Select(m => new List<SubgroupConjugates<T>>() { g, m }).ToList();
        var sz = 0;
        while (all.Sum(s => s.Count) != sz)
        {
            sz = all.Sum(s => s.Count);
            var tmp = all.ToList();
            all.Clear();
            foreach (var serie in tmp)
            {
                var last = serie.Last();
                if (last.Order == 1)
                {
                    all.Add(serie);
                    continue;
                }

                foreach (var m in MaximalSubgroups(last))
                    all.Add(serie.Append(m).ToList());
            }
        }

        return all.OrderBy(s => s.Count).ToArray();
    }

    private List<SubgroupConjugates<T>> GetChiefSerie(List<SubgroupConjugates<T>> serie) => serie.Where(cj => cj.IsNormal).ToList();

    private List<SubgroupConjugates<T>> GetCompositionSerie(List<SubgroupConjugates<T>> serie)
    {
        var compSerie = new List<SubgroupConjugates<T>>();
        var queue = new Queue<SubgroupConjugates<T>>(serie);
        while (queue.TryDequeue(out var hi))
        {
            compSerie.Add(hi);
            while (queue.TryPeek(out var hj))
            {
                var K = hj.Representative;
                var setK = K.ToSet();
                var G = hi.Conjugates.First(sg => sg.SuperSetOf(K));
                var op = Group.ByConjugateSet(G);
                if (G.All(g => setK.SetEquals(op(g, setK))))
                    break;
                else
                    queue.Dequeue();
            }
        }

        return compSerie;
    }

    public Dictionary<SerieType, Serie<T>[]> AllSeries()
    {
        Naming();
        var digits = AllSubgroupConjugates.Max(sg => $"{sg.FullName}".Length) + 1;
        var series = MaximalSubgroupSeries().Select(s => new Serie<T>(s, SerieType.Serie, digits)).ToArray();
        var chiefSeries = new HashSet<Serie<T>>();
        var compSeries = new HashSet<Serie<T>>();
        foreach (var serie in series)
        {
            var chief = new Serie<T>(GetChiefSerie(serie.Content), SerieType.Chief, digits);
            chief.AddTo(chiefSeries);
            var comp = new Serie<T>(GetCompositionSerie(serie.Content), SerieType.Composition, digits);
            comp.AddTo(compSeries);
        }

        compSeries.ExceptWith(chiefSeries);

        var allSeries = new Dictionary<SerieType, Serie<T>[]>()
        {
            [SerieType.Serie] = series,
            [SerieType.Chief] = chiefSeries.OrderBy(s => s.Count).ToArray(),
            [SerieType.Composition] = compSeries.OrderBy(s => s.Count).ToArray()
        };
        return allSeries;
    }

    public void DisplayAllSeries()
    {
        var allSeries = AllSeries();
        
        Console.WriteLine(Parent.ShortName);
        Console.WriteLine(Infos);
        
        Console.WriteLine();
        var conjSeries = allSeries[SerieType.Serie];
        var chiefSeries = allSeries[SerieType.Chief];
        var compSeries = allSeries[SerieType.Composition];
        
        conjSeries.OrderBy(e => e.Count).Println("Lattice Subgroups");
        Console.WriteLine($"Total:{conjSeries.Length}");
        Console.WriteLine();
        
        if (chiefSeries.Length != 0)
            chiefSeries.OrderBy(s => s.Count).Println("Chief Series");
        
        if (compSeries.Length != 0)
            compSeries.OrderBy(s => s.Count).Println("Composition Series");
        
        Console.WriteLine();
    }

    public ConcreteGroup<T> FrattiniSubGroup()
    {
        var elts = Parent.ToHashSet();
        foreach (var sg in MaximalSubgroups().SelectMany(cj => cj.Conjugates))
            elts.IntersectWith(sg);
        return AllRepresentatives.First(sg => sg.SetEquals(elts));
    }

    public ConcreteGroup<T> FittingSubGroup()
    {
        if (Parent.Count() == 1)
            return AllSubgroupConjugates.First(cj => cj.IsTrivial).Representative;

        return AllSubgroupConjugates.Where(cj => cj.IsProperNormal)
            .Descending()
            .Select(cj => cj.Representative)
            .First(sg => Group.IsNilpotent(sg));
    }

    public IEnumerable<ConcreteGroup<T>> All => AllSubgroupConjugates.SelectMany(sc => sc.Conjugates);
    public IEnumerable<ConcreteGroup<T>> AllRepresentatives => AllSubgroupConjugates.Select(sc => sc.Representative);
    public override string ToString() => Name;
}
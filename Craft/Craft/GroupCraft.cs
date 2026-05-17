using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.Subgroups;
using FastGoat.UserGroup.Perms;

namespace Craft.Craft;

public static class GroupCraft
{
    public static List<GroupSubset<T>> SubGroupsConjugates<T>(ConcreteGroup<T> g, GroupSubset<T> h)
        where T : struct, IElt<T>
    {
        if (!h.SubSetOf(g))
            throw new GroupException(GroupExceptionType.NotSubGroup);

        if (g.GroupType == GroupType.AbelianGroup)
            return [h];

        var all = new List<GroupSubset<T>>();
        var gens = g.GetGenerators().ToHashSet();
        var conjsH = Group.Orbits(gens, Group.ByConjugateSet(g), h);
        foreach (var cj in conjsH)
        {
            var gensConjs = cj.Generators.ToArray();
            var elts = Group.GenerateElements(g, gensConjs);
            all.Add(new([..gensConjs], elts));
        }

        return all;
    }

    public static T[] RecreateGenerators<T>(ConcreteGroup<T> g, T[] generators, Comparer<T> comp)
        where T : struct, IElt<T>
    {
        HashSet<T> tmpElements = new() { g.Neutral() };
        List<T> newGens = new();
        var set = generators.ToHashSet();
        set.ExceptWith(tmpElements);
        var gens = g.GetGenerators().ToHashSet();
        do
        {
            newGens.Add(set.OrderBy(e => e, comp).First());
            tmpElements = Group.GenerateElements(g, tmpElements, newGens);
            set.ExceptWith(tmpElements.SelectMany(x => Group.Orbits(gens, Group.ByConjugate(g), x)));
        } while (g.Count() != tmpElements.Count && set.Any());

        if (g.Count() != tmpElements.Count)
            return [];

        return newGens.ToArray();
    }

    public static T[] RecreateGenerators<T>(ConcreteGroup<T> g, T[] generators) where T : struct, IElt<T>
    {
        var idx = generators.Index().ToDictionary(e => e.Item, e => e.Index);
        var comp = Comparer<T>.Create((a, b) => idx[a].CompareTo(idx[b]));
        return RecreateGenerators(g, generators, comp);
    }

    public static T[] RecreateGenerators<T>(ConcreteGroup<T> g) where T : struct, IElt<T>
    {
        return RecreateGenerators(g, g.ToArray());
    }

    public static T[] RecreateGenerators<T>(ConcreteGroup<T> g, Comparer<T> comp) where T : struct, IElt<T>
    {
        return RecreateGenerators(g, g.ToArray(), comp);
    }

    public static HashSet<T> GenerateElementsLimited<T>(IGroup<T> bg, T[] generators, int limits)
        where T : struct, IElt<T>
    {
        var elements = new HashSet<T>([bg.Neutral()]);
        var q = new Queue<T>(elements);
        HashSet<T> generatedElements = new HashSet<T>(elements);
        while (q.Count != 0 && limits >= generatedElements.Count)
        {
            var e1 = q.Dequeue();
            foreach (var e2 in generators)
            {
                var e3 = bg.Op(e2, e1);
                if (generatedElements.Add(e3))
                    q.Enqueue(e3);

                if (limits < generatedElements.Count)
                    return [];
            }
        }

        if (limits < generatedElements.Count)
            return [];

        return generatedElements;
    }

    public static HashSet<T> GenerateElementsLimited<T>(IGroup<T> bg, HashSet<T> elements, T[] generators, int limits)
        where T : struct, IElt<T>
    {
        if (!elements.Contains(bg.Neutral()))
            throw new GroupException(GroupExceptionType.BaseGroup);

        var q = new Queue<T>(elements);
        HashSet<T> generatedElements = new HashSet<T>(elements);
        while (q.Count != 0)
        {
            var e1 = q.Dequeue();
            foreach (var e2 in generators)
            {
                var e3 = bg.Op(e2, e1);
                if (generatedElements.Add(e3))
                    q.Enqueue(e3);

                if (limits < generatedElements.Count)
                    return [];
            }
        }

        if (limits < generatedElements.Count)
            return [];

        return generatedElements;
    }

    public static Dictionary<ConcreteGroup<T>, List<GroupSubset<T>>> AllFactors<T>(ConcreteGroup<T> g,
        ConcreteGroup<T> h)
        where T : struct, IElt<T>
    {
        var og = g.Count();
        var oh = h.Count();
        var ok = og / oh;
        var h0 = h.Except([h.Neutral()]).ToHashSet();
        var allSubGrs = new HashSet<GroupSubset<T>>(og * og);
        var table = new Dictionary<ConcreteGroup<T>, List<GroupSubset<T>>>(og * og, new GroupSetEquality<T>());
        var tablePCycles = new Dictionary<GroupSubset<T>, List<GroupSubset<T>>>(og * og);

        foreach (var e0 in g.Except(h).Where(e => ok % g.ElementsOrders[e] == 0))
        {
            var cyc = Group.Generate(g, e0);
            var oc = cyc.Count();
            var setCyc = cyc.ToSet();
            if (!allSubGrs.Contains(setCyc))
            {
                var conjs = table[cyc] = GroupCraft.SubGroupsConjugates(g, setCyc);
                allSubGrs.UnionWith(conjs);
                if (IntExt.PrimesDec(oc).Count == 1)
                {
                    foreach (var c in conjs)
                        tablePCycles[c] = conjs.ToList();
                }
            }
        }

        var sgsRem = table.Keys.Select(c => c.ToSet()).ToHashSet();
        var cycRem = tablePCycles.Values.SelectMany(e => e).ToHashSet();
        while (sgsRem.Count != 0)
        {
            var sgs = sgsRem.OrderBy(sg => sg.Count()).ToArray();
            var cyc = cycRem.OrderBy(sg => sg.Count()).ToArray();
            sgsRem.Clear();
            cycRem.Clear();
            foreach (var (sg0, sg1) in cyc.Grid2D(sgs))
            {
                if (sg0.SuperSetOf(sg1) || sg0.SubSetOf(sg1) || sg0.Overlaps(h0) || sg1.Overlaps(h0))
                    continue;

                var gens = sg0.Generators.Union(sg1.Generators).ToArray();
                var elts = GroupCraft.GenerateElementsLimited(g, sg1.Concat(sg0).ToHashSet(), gens, ok);
                if (elts.Count == 0 || ok % elts.Count != 0 || elts.Overlaps(h0))
                    continue;

                var set = new GroupSubset<T>(gens.ToHashSet(), elts);
                if (!allSubGrs.Contains(set))
                {
                    var sg2 = Group.Generate("K", g, gens);
                    var conjsSg2 = table[sg2] = GroupCraft.SubGroupsConjugates(g, set);
                    allSubGrs.UnionWith(conjsSg2);
                    sgsRem.Add(sg2.ToSet());
                    cycRem.UnionWith(tablePCycles[sg0]);
                }
            }
        }

        return table.Where(e => e.Key.Count() == ok && !e.Key.ToSet().Overlaps(h0)).ToDictionary();
    }
}
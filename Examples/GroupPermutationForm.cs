using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.Subgroups;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Perms;

namespace Examples;

public static class GroupPermutationForm
{
    static GroupPermutationForm()
    {
        Perm.Style = DisplayPerm.CyclesComma;
    }

    static HashSet<T> GenerateElementsLimited<T>(IGroup<T> bg, T[] generators, int limits)
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

    static Dictionary<ConcreteGroup<T>, List<GroupSubset<T>>> RestrictionConjugates<T>(ConcreteGroup<T> g,
        List<GroupSubset<T>> conjugates) where T : struct, IElt<T>
    {
        var lt = conjugates.Where(h => h.SubSetOf(g)).ToHashSet();
        var all = new Dictionary<ConcreteGroup<T>, List<GroupSubset<T>>>(conjugates.Count);
        var gens = g.GetGenerators().ToHashSet();
        var act = Group.ByConjugateSet(g);
        while (lt.Count != 0)
        {
            var sg = lt.First();
            var conjs = Group.Orbits(gens, act, sg);
            var subConjs = lt.Where(c0 => conjs.Any(c1 => c0.SetEquals(c1))).ToList();
            var g1 = Group.Generate("H", g, sg.Generators.ToArray());
            all[g1] = subConjs;
            lt.ExceptWith(subConjs);
        }

        return all;
    }

    static Dictionary<ConcreteGroup<T>, List<GroupSubset<T>>> Restriction<T>(ConcreteGroup<T> G,
        Dictionary<ConcreteGroup<T>, List<GroupSubset<T>>> allSubgroups) where T : struct, IElt<T>
    {
        var setG = G.ToSet();
        return allSubgroups.Where(e => e.Value.Any(f => f.SubSetOf(setG)))
            .SelectMany(e => RestrictionConjugates(G, e.Value)).ToDictionary();
    }

    static (Sn sn, Perm[] gens) RemovePassiveCycles(Perm[] generators)
    {
        var gens = generators.ToArray();
        while (true)
        {
            var prepGens = gens.ToDictionary(g => g, g => g.DisjoinCycles);
            var sn0 = gens[0].Sn;
            var dijCycles = gens.Select(e => (Key: e, Value: e.DisjoinCycles.ToXSet())).ToArray();
            var candidates = dijCycles.Select(e => (e.Key,
                e.Value.Where(f => e.Value.Any(g => !g.Equals(f) && g.Order % f.Order == 0)).ToXSet()));
            var passiveCycles = new Dictionary<Perm, XSet<int>>();
            foreach (var (k, v) in candidates)
            {
                if (v.Count == 0)
                    continue;

                foreach (var c in v)
                {
                    var set = c.Orbits.Where(o => o.Length > 1).SelectMany(o => o).ToXSet();
                    if (prepGens.All(g => g.Value.Contains(c) || set.All(i => g.Key.Table[i] == i)))
                    {
                        passiveCycles[k] = set;
                        break;
                    }
                }

                if (passiveCycles.Count > 0)
                    break;
            }

            if (passiveCycles.Count > 0)
            {
                var (g, s) = passiveCycles.First();
                var (N, n) = (g.Sn.N, s.Count);
                var swap = N.SeqLazy().Where(i => !s.Contains(i)).Concat(s).Index()
                    .ToDictionary(e => e.Item, e => e.Index);

                var sn1 = new Sn(N - n);
                gens = gens.Select(e =>
                        e.Orbits.Where(o => o.Length > 1 && !s.Equals(new XSet<int>(o)))
                            .Select(o => o.Select(i => swap[i]).ToArray()).ToArray())
                    .Select(e => e.Aggregate(sn1.Neutral(), (acc, o) => acc * sn1.CycleP1(o)))
                    .ToArray();
            }

            var sn2 = gens[0].Sn;
            if (sn2.N == sn0.N)
                return (sn2, gens);
        }
    }

    static bool IsIsomorphicTo<T1, T2>(Dictionary<ConcreteGroup<T1>, List<GroupSubset<T1>>> subGroups,
        ConcreteGroup<T2> g2)
        where T1 : struct, IElt<T1>
        where T2 : struct, IElt<T2>
    {
        var g1 = subGroups.Keys.MaxBy(e => e.Order)!;
        if (g2.Order != g1.Order || g1.GroupType != g2.GroupType)
            return false;

        if (!g1.ElementsOrdersList().Ascending().SequenceEqual(g2.ElementsOrdersList().Ascending()))
            return false;

        if (g1.GroupType == GroupType.AbelianGroup && g2.GroupType == GroupType.AbelianGroup)
        {
            var gType = Group.AbelianGroupType(g1);
            var guType = Group.AbelianGroupType(g2);
            return gType.SequenceEqual(guType);
        }

        var Infos = (subGroups.Values.Sum(s => s.Count), subGroups.Count, subGroups.Count(s => s.Value.Count == 1));
        if (!Infos.Equals(g2.AllSubgroups().Infos.ToTuples()))
            return false;

        return Group.AllMorphisms(g1, g2, Group.MorphismType.Isomorphism).Any();
    }

    public static ConcreteGroup<Perm> PermutationForm<T>(ConcreteGroup<T> G) where T : struct, IElt<T>
    {
        return PermutationForm(Group.AllSubGroups(G)
            .ToDictionary(e => e.Key, e => e.Value.Select(f => f.ToSet()).ToList()));
    }

    public static ConcreteGroup<Perm> PermutationForm<T>(Dictionary<ConcreteGroup<T>, List<GroupSubset<T>>> subGroups)
        where T : struct, IElt<T>
    {
        var G = subGroups.MaxBy(g => g.Key.Order).Key;
        if (G.GroupType == GroupType.AbelianGroup)
            return FG.AbelianPerm(Group.AbelianGroupType(G));

        var ordG = G.Order;
        var byOrder = subGroups.Keys.GroupBy(e => e.Order).ToDictionary(
            e => e.Key,
            e => e.OrderBy(f => subGroups[f].Count).ToArray());
        var (H, Facts) = subGroups.Keys
            .Where(e => e.Order != 1 && e.Order != ordG && subGroups[e].Count == 1).OrderBy(g => g.Order).ToArray()
            .ToDictionary(d => d, d =>
            {
                if (byOrder.ContainsKey(ordG / d.Order))
                    return byOrder[ordG / d.Order].Where(e => e.Intersect(d).Count() == 1).Take(1).ToList();

                return [];
            })
            .Where(e => e.Value.Count > 0)
            .OrderByDescending(f => f.Value.Take(1).Count(g => subGroups[g].Count == 1))
            .ThenByDescending(f => f.Value.Take(1).Sum(g => g.Order))
            .ThenBy(f => f.Key.Order)
            .FirstOrDefault(_ => true, new(G, new()));

        if (Facts.Count == 0)
        {
            if (G.Order == 60)
            {
                var A5 = FG.Alternate(5);
                if (G.IsIsomorphicTo(A5))
                    return A5;
            }

            return G.ToPermGroup().Item1;
        }

        var K = Facts.First();
        if (Group.IsNormalSubgroup(G, K))
        {
            var h = PermutationForm(Restriction(H, subGroups));
            var k = PermutationForm(Restriction(K, subGroups));
            var nh = h.Neutral().Sn.N;
            var nk = k.Neutral().Sn.N;
            var hGens = h.GetGenerators().Select(e => FG.PaddingRight(e, nk));
            var kGens = k.GetGenerators().Select(e => FG.PaddingLeft(e, nh));
            var (sn, hkGens) = RemovePassiveCycles(hGens.Concat(kGens).ToArray());
            var HK = Group.Generate(G.Name, sn, hkGens);
            return HK;
        }
        else
        {
            var k = PermutationForm(Restriction(K, subGroups));
            var (h, autHpg, _) = FG.RegPermAutGroup(PermutationForm(Restriction(H, subGroups)));
            var snH = h.Neutral().Sn;
            var nK = k.Neutral().Sn.N;

            var homKAutH = Group.AllHomomorphisms(k, autHpg);
            foreach (var hom in homKAutH.OrderBy(e => e.Kernel().Count()))
            {
                if (hom.Kernel().Count() == 1)
                {
                    var hkGens = h.GetGenerators().Concat(hom.Image()).ToArray();
                    if (GenerateElementsLimited(snH, hkGens, ordG).Count == ordG)
                    {
                        var HK = Group.Generate(G.Name, snH, hkGens);
                        if (IsIsomorphicTo(subGroups, HK))
                            return HK;
                    }
                }

                {
                    var kGens = k.GetGenerators().Select(e => FG.ConcatPerm(hom[e], e)).ToArray();
                    var hGens = h.GetGenerators().Select(f => FG.PaddingRight(f, nK)).ToArray();
                    var (sn, hkGens) = RemovePassiveCycles(hGens.Concat(kGens).ToArray());
                    if (GenerateElementsLimited(sn, hkGens, ordG).Count != ordG)
                        continue;

                    var HK = Group.Generate(G.Name, sn, hkGens);
                    if (IsIsomorphicTo(subGroups, HK))
                        return HK;
                }
            }

            Console.WriteLine($"NOT FOUND  H:{H.ShortName} K:{K.ShortName} {autHpg.ShortName}");
            Console.WriteLine($"H in {h.Neutral().Sn} K in {k.Neutral().Sn}");
            throw new($"###2 {G.ShortName}");
        }
    }

    public static void ExampleAbelian()
    {
        for (int ord = 1; ord <= 128; ord++)
        {
            foreach (var abType in AbelianExt.AllAbTypes(ord))
            {
                var ab = FG.AbelianPerm(abType.can);
                DisplayGroup.HeadGenerators(ab);
                var abType2 = Group.AbelianGroupType(ab);
                if (!abType.can.SequenceEqual(abType2))
                    throw new();
            }
        }
    }

    public static void ExampleDihedral()
    {
        for (int n = 3; n < 33; n++)
        {
            var D2pg = FG.Dihedral(n);
            var D2n = FG.DihedralGL2p(n);
            DisplayGroup.HeadGenerators(D2pg);
            DisplayGroup.AreIsomorphics(D2pg, D2n);
            Console.WriteLine();
            if (!D2n.IsIsomorphicTo(D2pg))
                throw new();
        }
    }

    public static void ExampleDicyclic()
    {
        for (int m = 3; m < 33; m++)
        {
            var Dic_m = FG.DicyclicGL2p(m);
            var Dic_m_pg = FG.DicyclicPg(m);
            DisplayGroup.HeadGenerators(Dic_m_pg);
            DisplayGroup.AreIsomorphics(Dic_m_pg, Dic_m);
            Console.WriteLine();
            if (!Dic_m_pg.IsIsomorphicTo(Dic_m))
                throw new();
        }
    }

    public static void ExampleSemiDihedralAndModularMax()
    {
        for (int n = 4; n < 9; n++)
        {
            var qdpg = FG.SemiDihedralPg(n);
            DisplayGroup.HeadGenerators(qdpg);
            if (!qdpg.IsIsomorphicTo(FG.SemiDihedralGL2p(n)))
                throw new();

            var mmpg = FG.ModularMaxPg(n);
            DisplayGroup.HeadGenerators(mmpg);
            if (!mmpg.IsIsomorphicTo(FG.ModularMaxGL2p(n)))
                throw new();

            Console.WriteLine();
        }
    }

    public static void ExampleAllMetaCyclicSemiDirectProducts()
    {
        GlobalStopWatch.Restart();
        for (int ord = 2; ord < 129; ord++)
        {
            foreach (var (m, n, r) in FG.MetaCyclicCoefs(ord))
            {
                var mtpg = FG.MetaCyclicPg(m, n, r);
                var mt = FG.MetaCyclicSdp(m, n, r);
                DisplayGroup.HeadGenerators(mtpg);
                if (!mt.IsIsomorphicTo(mtpg))
                    throw new();
            }
        }

        GlobalStopWatch.Show();
    }

    static void MetaCyclicPermFormUpTo(int maxOrd)
    {
        var (errors, total) = (0, 0);
        GlobalStopWatch.Restart();

        for (int ord = 1; ord <= maxOrd; ord++)
        {
            foreach (var (m, n, r) in FG.MetaCyclicCoefs(ord))
            {
                ++total;
                var G = FG.MetaCyclicSdpWg(m, n, r);
                var Gpg = PermutationForm(G);
                DisplayGroup.HeadOrders(Gpg);
                DisplayGroup.Generators(Gpg, showBaseGroup: true);
                var sn = Gpg.Neutral().Sn;
                var sn0 = FG.MetaCyclicGens(m, n, r).sn;
                var test = sn.N == sn0.N;
                Console.WriteLine($"{Gpg.ShortName} in {sn} expected {sn0} Test:{test}");
                if (!test)
                    ++errors;

                Console.WriteLine();
            }
        }

        GlobalStopWatch.Show($"Errors:{errors}/{total} Max Order:{maxOrd}");
    }

    static void AllGroupPermFormUpTo(int maxOrd)
    {
        var total = 0;
        GlobalStopWatch.Restart();
        for (int ord = 1; ord <= maxOrd; ord++)
        {
            foreach (var G in FG.AllGroupsOfOrder(ord))
            {
                ++total;
                var Gpg = PermutationForm(G);
                DisplayGroup.HeadGenerators(Gpg);
                var sn = Gpg.Neutral().Sn;
                Console.WriteLine($"{G.ShortName} in {sn} Ratio:{G.Count() / (sn.N + 0.0):f3}");
                Console.WriteLine();
            }
        }

        GlobalStopWatch.Show($"Total:{total}  Max Order:{maxOrd}");
    }

    public static void ExampleAllMetaCyclicPermForm()
    {
        MetaCyclicPermFormUpTo(128);
    }

    public static void ExampleAllGroupsPermFormUpTo96()
    {
        AllGroupPermFormUpTo(96); // # Total:1024  Max Order:96 Time:5m13s
    }
}

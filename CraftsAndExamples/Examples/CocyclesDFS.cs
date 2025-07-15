using System.Numerics;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.UserGroup.GModuleN;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Matrix;

namespace CraftsAndExamples.Examples;

public static class CocyclesDFS
{
    public class OneCocyclesDFS<Tn, Tg> where Tn : struct, IElt<Tn> where Tg : struct, IElt<Tg>
    {
        public ConcreteGroup<Tg> G { get; }
        public ConcreteGroup<Tn> N { get; }
        public MapElt<Tg, Automorphism<Tn>> L { get; }

        public OneCocyclesDFS(ConcreteGroup<Tg> G0, ConcreteGroup<Automorphism<Tn>> AutN, Homomorphism<Tg, Automorphism<Tn>> L0)
        {
            (N, G) = (AutN.Neutral().AutGroup.G, G0);
            L = new MapElt<Tg, Automorphism<Tn>>(G, AutN, L0.HomMap.ToDictionary(e => e.Key, e => e.Value));
        }


        public OneCocyclesDFS(ConcreteGroup<Tn> N0, ConcreteGroup<Tg> G0, MapElt<Tg, Automorphism<Tn>> L0)
        {
            (N, G, L) = (N0, G0, L0);
        }

        public (Tg, Tn) OneCocycleExpression((Tg, Tn) e0, (Tg, Tn) e1)
        {
            // ω(rs) = L(r)[ω(s)] ω(r)

            var ((g0, n0), (g1, n1)) = (e0, e1);
            var g2 = G.Op(g0, g1);
            var n2 = N.Op(L[g0][n1], n0);
            return (g2, n2);
        }

        public List<MapElt<Tg, Tn>> All1Cocycles()
        {
            return SearchOneCocycles(new() { [G.Neutral()] = N.Neutral() }).Distinct().ToList();
        }

        private IEnumerable<MapElt<Tg, Tn>> SearchOneCocycles(Dictionary<Tg, Tn> current)
        {
            if (current.Count == G.Count())
            {
                var sol = new MapElt<Tg, Tn>(G, N, current);
                yield return sol;
            }

            foreach (var next in NextChildsOneCocycles(current))
            {
                foreach (var sol in SearchOneCocycles(next))
                    yield return sol;
            }
        }

        private Dictionary<Tg, Tn> OneCocycleUpdate(Dictionary<Tg, Tn> prev, Tg g0, Tn n0)
        {
            if (prev.ContainsKey(g0))
                return prev;

            var next = prev.ToDictionary(e => e.Key, e => e.Value);
            next[g0] = n0;
            var g1g2 = next.Grid2D(next).ToArray();
            foreach (var ((g1, n1), (g2, n2)) in g1g2)
            {
                var (g3, n3) = OneCocycleExpression((g1, n1), (g2, n2));
                if (next.TryGetValue(g3, out Tn n4))
                {
                    if (!n4.Equals(n3))
                        return prev;
                }
                else
                {
                    next[g3] = n3;
                }
            }

            return next;
        }

        private IEnumerable<Dictionary<Tg, Tn>> NextChildsOneCocycles(Dictionary<Tg, Tn> current)
        {
            var rems = G.Except(current.Keys).OrderByDescending(e => G.ElementsOrders[e]).ThenAscending().ToArray();
            var nb = current.Count;
            foreach (var (g0, n0) in rems.Grid2D(N))
            {
                var next = OneCocycleUpdate(current, g0, n0);
                if (next.Count == nb)
                    continue;
                else if (next.Count > nb)
                {
                    yield return next;
                }
                else
                    throw new();
            }
        }
    }

    static void All_1_Cocycles<Tn, Tg>(ConcreteGroup<Tn> N, ConcreteGroup<Tg> G)
        where Tg : struct, IElt<Tg>
        where Tn : struct, IElt<Tn>
    {
        var autN = Group.AutomorphismGroup(N);
        var allOps = Group.AllHomomorphisms(G, autN);
        Console.WriteLine($"N:{N.ShortName} and G:{G.ShortName}");
        var all = new Dictionary<MapElt<Tg, Automorphism<Tn>>, HashSet<MapElt<Tg, Tn>>>();
        foreach (var L in allOps)
        {
            var L0 = new MapElt<Tg, Automorphism<Tn>>(G, autN, new(L.HomMap));
            var homol = new OneCocyclesDFS<Tn, Tg>(N, G, L0);
            var allOneCocycles = homol.All1Cocycles();
            Console.WriteLine($"L:{L0}");
            Console.WriteLine($"Count:{allOneCocycles.Count}");
            all[L0] = allOneCocycles.ToHashSet();
            foreach (var co in allOneCocycles)
                Console.WriteLine($"  1co:{co}");

            Console.WriteLine();
        }

        Console.WriteLine($"Total:{all.Values.Sum(v => v.Count)}");
        Console.WriteLine();
    }

    public static void ExampleAll1Cocycle()
    {
        var (c2, c4, c2c2, c2c2c2) = (new Cn(2), new Cn(4), FG.Abelian(2, 2), FG.Abelian(2, 2, 2));
        All_1_Cocycles(c4, c4);
        // All_1_Cocycles(c2c2, c2c2c2);
    }

    public static void DisplayMapElt<Tn, Tg>(string title, params MapElt<Ep<Tg>, Tn>[] crMaps)
        where Tg : struct, IElt<Tg>
        where Tn : struct, IElt<Tn>
    {
        var c0 = crMaps[0];
        var G = (ConcreteGroup<Tg>)(((Gp<Tg>)c0.Domain.BaseGroup).Gi[0]);
        var table = c0.map.OrderKeys(G).SelectMany(e => crMaps.Select(c => c[e.Key]).Cast<object>().Prepend(e.Key)).ToArray();
        var mat = Ring.Matrix(c0.Map.Count, table);
        Console.WriteLine(title);
        Ring.DisplayMatrix(mat, sep: " ");
        Console.WriteLine();
    }

    public static Dictionary<MapElt<Tg, Automorphism<Tn>>, HashSet<MapElt<Ep<Tg>, Tn>>> All_2_Cocycles_N_G<Tn, Tg>(ConcreteGroup<Tn> N,
        ConcreteGroup<Tg> G, bool trivialActionOnly = true)
        where Tg : struct, IElt<Tg>
        where Tn : struct, IElt<Tn>
    {
        var autN = Group.AutomorphismGroup(N);
        var allOps = Group.AllHomomorphisms(G, autN);
        Console.WriteLine("Start Search All ext N -> E -> G");
        Console.WriteLine($"    with N:{N.ShortName} and G:{G.ShortName}");
        Ring.MatrixDisplayForm = Ring.MatrixDisplay.SquareBracket;
        var all = new Dictionary<MapElt<Tg, Automorphism<Tn>>, HashSet<MapElt<Ep<Tg>, Tn>>>();
        var lbl = 0;
        var trivL = new Homomorphism<Tg, Automorphism<Tn>>(G, G.ToDictionary(e => e, _ => autN.Neutral()));
        if (trivialActionOnly)
            allOps = new() { trivL };

        foreach (var L in allOps)
        {
            var L0 = new MapElt<Tg, Automorphism<Tn>>(G, autN, new(L.HomMap));
            var homol = new TwoCocyclesDFS<Tn, Tg>(N, G, L0, $"Lbl{++lbl}/{allOps.Count}");
            var all2Cocycles = homol.AllTwoCocycles();
            all[L0] = all2Cocycles.ToHashSet();
        }

        Console.WriteLine($"N:{N.ShortName} and G:{G.ShortName}");
        Console.WriteLine($"Total:{all.Values.Sum(v => v.Count)}");
        Console.WriteLine();

        return all;
    }

    public static void TwoCocyclesExamples()
    {
        var (c2, c4, c2c2, c8, c2c2c2, c4c4, c2c4, d8, q8, c16) = (new Cn(2), new Cn(4), FG.Abelian(2, 2), new Cn(8),
            FG.Abelian(2, 2, 2), FG.Abelian(4, 4), FG.Abelian(4, 2), FG.Dihedral(4), FG.Quaternion(8), new Cn(16));

        // Twisted Actions for trivial action of G by N
        // it takes longuer time than before but it is more accure
        All_2_Cocycles_N_G(c4, c4);
        All_2_Cocycles_N_G(c2c2, c4);
        All_2_Cocycles_N_G(c2c4, c4);

        All_2_Cocycles_N_G(c4, c2c2);
        All_2_Cocycles_N_G(c2c2, c2c2);
        All_2_Cocycles_N_G(c2c4, c2c2);

        All_2_Cocycles_N_G(c8, c2);
        All_2_Cocycles_N_G(d8, c2);
        All_2_Cocycles_N_G(c16, c2);
    }

    static HashSet<ExtensionGroup<Tn, Tg>> AllExt<Tn, Tg>(ConcreteGroup<Tn> N, ConcreteGroup<Tg> G)
        where Tg : struct, IElt<Tg>
        where Tn : struct, IElt<Tn>
    {
        var all = All_2_Cocycles_N_G(N, G, trivialActionOnly: false);
        var allExt = new HashSet<ExtensionGroup<Tn, Tg>>();
        foreach (var (L, tws) in all)
        {
            foreach (var w in tws)
            {
                var ext = Group.ExtensionGroup(N, L, w, G);
                if (!ext.ExtBase.IsGroup)
                    continue;

                allExt.Add(ext);
            }
        }

        return allExt;
    }

    public static IEnumerable<(ConcreteGroup<Ep2<Tn, Ep<ZnInt>>>, (int, int, int))> BuildExtensions<Tn>(HashSet<ConcreteGroup<Tn>> gr)
        where Tn : struct, IElt<Tn>
    {
        foreach (var tuple in BuildExtensions(new HashSet<ConcreteGroup<Ep<ZnInt>>>() { FG.Abelian(2) }, gr))
            yield return tuple;
    }

    public static IEnumerable<(ConcreteGroup<Ep2<Tn, Tg>>, (int, int, int))> BuildExtensions<Tn, Tg>(HashSet<ConcreteGroup<Tg>> Gs,
        HashSet<ConcreteGroup<Tn>> Ns)
        where Tn : struct, IElt<Tn>
        where Tg : struct, IElt<Tg>
    {
        var allNg = Ns.Grid2D(Gs).ToArray();
        foreach (var tuple in BuildExtensions(allNg))
            yield return tuple;
    }

    public static (int, int, int) SubGroupsDetails<T>(Dictionary<ConcreteGroup<T>, List<ConcreteGroup<T>>> sg0)
        where T : struct, IElt<T>
    {
        return (sg0.Values.Sum(s => s.Count), sg0.Count, sg0.Count(s => s.Value.Count == 1));
    }

    public static (int, int, int) SubGroupsDetails<T>(ConcreteGroup<T> g) where T : struct, IElt<T>
    {
        return SubGroupsDetails(Group.AllSubGroups(g));
    }

    public static IEnumerable<(ConcreteGroup<Ep2<Tn, Tg>>, (int, int, int))> BuildExtensions<Tn, Tg>(
        (ConcreteGroup<Tn> N, ConcreteGroup<Tg> G)[] AllNG)
        where Tn : struct, IElt<Tn>
        where Tg : struct, IElt<Tg>
    {
        var allExts = AllNG.SelectMany(e => AllExt(e.N, e.G)).ToList();
        foreach (var ext in BuildExtensions(allExts))
            yield return ext;
    }

    public static IEnumerable<(ConcreteGroup<Ep2<Tn, Tg>>, (int, int, int))> BuildExtensions<Tn, Tg>(
        List<ExtensionGroup<Tn, Tg>> allExts)
        where Tn : struct, IElt<Tn>
        where Tg : struct, IElt<Tg>
    {
        var infosSubGroups =
            new Dictionary<ConcreteGroup<Ep2<Tn, Tg>>,
                Dictionary<ConcreteGroup<Ep2<Tn, Tg>>, List<ConcreteGroup<Ep2<Tn, Tg>>>>>();

        var abTypes = new HashSet<int[]>(new SequenceEquality<int>());
        var listIsos = new HashSet<ConcreteGroup<Ep2<Tn, Tg>>>(new IsomorphEquality<Ep2<Tn, Tg>>());

        foreach (var (ext, i) in allExts.Select((ext, i) => (ext, i + 1)))
        {
            Console.Write($"Progress:{i,5} / {allExts.Count} Found:{listIsos.Count}");
            Console.CursorLeft = 0;

            if (ext.GroupType == GroupType.AbelianGroup)
            {
                try
                {
                    var abType = AbelianInvariantsFactors.Reduce(ext).Order().ToArray();
                    if (!abTypes.Add(abType))
                        continue;
                }
                catch (Exception e)
                {
                    Console.WriteLine("@???????????????????????????????????????????????????????????");
                    // TODO why some twisted actions wont give a group
                }
            }

            if (!Group.IsGroup(ext))
            {
                Console.WriteLine("@@??????????????????????????????????????????????????????????");
                continue;
            }

            if (!listIsos.Add(ext))
                continue;

            var allSubs = Group.AllSubGroups(ext);
            var info0 = SubGroupsDetails(allSubs);

            yield return (ext, info0);

            infosSubGroups[ext] = allSubs;
        }

        Console.WriteLine();
        Console.WriteLine($"AllExts Found : {infosSubGroups.Count}");
    }

    public static void DisplayInfosGroups<Tg>((ConcreteGroup<Tg>, (int, int, int))[] elts, int countStart = 0, bool naming = true,
        string prefix = "Sm")
        where Tg : struct, IElt<Tg>
    {
        var names = elts.Select((e, i) => (e.Item1, i + 1))
            .ToDictionary(e => e.Item1, e => naming ? $"{prefix}{e.Item1.Count()}[{e.Item2}]" : e.Item1.Name);
        var maxLt = names.Max(e => e.Value.Length);
        var lt = Enumerable.Repeat('#', maxLt + 4).Glue();
        var line = $"#################{lt}#################";
        var fmt = $"#################  {{0,{-maxLt}}}  #################";
        foreach (var (g, k, infos) in elts.Select((e, k) => (e.Item1, countStart + k + 1, e.Item2)))
        {
            var name = names[g];
            Console.WriteLine(line);
            Console.WriteLine(fmt, name);
            Console.WriteLine(line);
            g.Name = name;
            DisplayGroup.HeadOrders(g);
            Console.WriteLine($"AllSubGr:{infos.Item1} AllConjsCl:{infos.Item2} AllNorms:{infos.Item3}");
            Console.WriteLine();
        }
    }

    public static void DisplayInfosGroups<Tg>(IEnumerable<ConcreteGroup<Tg>> gs, bool naming = false, string prefix = "Sm")
        where Tg : struct, IElt<Tg>
    {
        DisplayInfosGroups(
            gs.Select(e => (e, SubGroupsDetails(e)))
                .OrderBy(e => e.Item1.GroupType)
                .ThenByDescending(e => e.Item1.ElementsOrders.Values.Max())
                .ThenBy(e => e.Item2).ToArray()
            , naming: naming, prefix: prefix);
    }

    public static void ExampleAll16Orders()
    {
        var allOrder4 = BuildExtensions(new HashSet<ConcreteGroup<ZnInt>>() { new Cn(2) })
            .OrderBy(e => e.Item1.GroupType)
            .ThenByDescending(e => e.Item1.ElementsOrders.Values.Max())
            .ThenBy(e => e.Item2).ToArray();

        var allOrder8 = BuildExtensions(allOrder4.Select(e => e.Item1)
                .ToHashSet(new IsomorphEquality<Ep2<ZnInt, Ep<ZnInt>>>()))
            .OrderBy(e => e.Item1.GroupType)
            .ThenByDescending(e => e.Item1.ElementsOrders.Values.Max())
            .ThenBy(e => e.Item2).ToArray();

        var allOrder16 = BuildExtensions(allOrder8.Select(e => e.Item1)
                .ToHashSet(new IsomorphEquality<Ep2<Ep2<ZnInt, Ep<ZnInt>>, Ep<ZnInt>>>()))
            .OrderBy(e => e.Item1.GroupType)
            .ThenByDescending(e => e.Item1.ElementsOrders.Values.Max())
            .ThenBy(e => e.Item2).ToArray();

        Console.Clear();

        DisplayInfosGroups(allOrder4);
        DisplayInfosGroups(allOrder8);
        DisplayInfosGroups(allOrder16);

        Console.Beep();
    }

    public static void ExampleAll32Orders()
    {
        var allAb16 = IntExt.Partitions32[4].Select(l => FG.Abelian(l.Select(k => 2.Pow(k)).ToArray()))
            .ToHashSet(); // 5 groups

        var allSdp16 = Group.AllSemiDirectProd("C8x:C2", FG.Abelian(8), new Cn(2))
            .Concat(Group.AllSemiDirectProd("(C4xC2)x:C2", FG.Abelian(4, 2), new Cn(2)))
            .Concat(Group.AllSemiDirectProd("C4x:C4", FG.Abelian(4), new Cn(4)))
            .ToHashSet(new IsomorphEquality<Ep2<Ep<ZnInt>, ZnInt>>()); // 7 groups

        var c2 = new Cn(2);
        var c2q8 = Product.Generate(new Cn(2), FG.Quaternion(8));
        var q16 = Product.Generate(Group.Generate("C1", c2, c2[0]), FG.Quaternion(16));
        q16.Name = "Q16";

        var l0 = new List<(ConcreteGroup<Ep2<Ep2<ZnInt, Mat>, Ep<ZnInt>>>, (int, int, int))>();
        var l1 = new List<(ConcreteGroup<Ep2<Ep2<Ep<ZnInt>, ZnInt>, Ep<ZnInt>>>, (int, int, int))>();
        var l2 = new List<(ConcreteGroup<Ep2<Ep<ZnInt>, Ep<ZnInt>>>, (int, int, int))>();
        foreach (var (g, infos) in BuildExtensions(new HashSet<ConcreteGroup<Ep2<ZnInt, Mat>>>() { q16, c2q8 }))
            l0.Add((g, infos));

        foreach (var (g, infos) in BuildExtensions(allSdp16))
        {
            if (l0.Count + l1.Count == 51)
                break;

            if (l0.Any(g0 => g.IsIsomorphicTo(g0.Item1)))
                continue;

            l1.Add((g, infos));
        }

        foreach (var (g, infos) in BuildExtensions(allAb16))
        {
            if (l0.Any(g0 => g.IsIsomorphicTo(g0.Item1)))
                continue;

            if (l1.Any(g0 => g.IsIsomorphicTo(g0.Item1)))
                continue;

            l2.Add((g, infos));
            if (l0.Count + l1.Count + l2.Count == 51)
                break;
        }

        DisplayInfosGroups(l0.ToArray(), countStart: 0);
        DisplayInfosGroups(l1.ToArray(), countStart: l0.Count);
        DisplayInfosGroups(l2.ToArray(), countStart: l0.Count + l1.Count);

        Console.Beep(); // ~50min
    }

    public static void AllGroupsOrder_12_24()
    {
        var c2 = FG.Abelian(2);
        var c3 = FG.Abelian(3);
        var ab4 = new[] { FG.Abelian(2, 2), FG.Abelian(4) };

        var allOrder8 = BuildExtensions(ab4.ToHashSet())
            .OrderBy(e => e.Item1.GroupType)
            .ThenByDescending(e => e.Item1.ElementsOrders.Values.Max())
            .ThenBy(e => e.Item2).ToArray();

        var tuple12 = ab4.SelectMany(e => new[] { (c3, e), (e, c3) }).ToArray();
        var allOrder12 = BuildExtensions(tuple12)
            .OrderBy(e => e.Item1.GroupType)
            .ThenByDescending(e => e.Item1.ElementsOrders.Values.Max())
            .ThenBy(e => e.Item2).ToArray();

        var tuple24 = allOrder12.Select(e => (e.Item1, c2))
            .Concat(allOrder8.Select(e => (e.Item1, c3)))
            .ToArray();
        var allOrder24 = BuildExtensions(tuple24)
            .OrderBy(e => e.Item1.GroupType)
            .ThenByDescending(e => e.Item1.ElementsOrders.Values.Max())
            .ThenBy(e => e.Item2).ToArray();

        Console.Clear();

        DisplayInfosGroups(allOrder8);
        DisplayInfosGroups(allOrder12);
        DisplayInfosGroups(allOrder24);

        Console.Beep();
    }

    public static void AllGroupsOrder_10_20_40()
    {
        var c5 = FG.Abelian(5);
        var c10 = FG.Abelian(2, 5);
        var ab4 = new[] { FG.Abelian(2, 2), FG.Abelian(4) };
        var tuples20 = ab4.Select(g => (c5, g)).ToArray();
        var tuples40 = ab4.Select(g => (c10, g)).ToArray();

        var allOrder10 = BuildExtensions(new[] { (c5, FG.Abelian(2)) })
            .OrderBy(e => e.Item1.GroupType)
            .ThenByDescending(e => e.Item1.ElementsOrders.Values.Max())
            .ThenBy(e => e.Item2).ToArray();

        var allOrder20 = BuildExtensions(tuples20)
            .OrderBy(e => e.Item1.GroupType)
            .ThenByDescending(e => e.Item1.ElementsOrders.Values.Max())
            .ThenBy(e => e.Item2).ToArray();

        var allOrder40 = BuildExtensions(tuples40)
            .OrderBy(e => e.Item1.GroupType)
            .ThenByDescending(e => e.Item1.ElementsOrders.Values.Max())
            .ThenBy(e => e.Item2).ToArray();

        Console.Clear();

        DisplayInfosGroups(allOrder10);
        DisplayInfosGroups(allOrder20);
        DisplayInfosGroups(allOrder40);

        Console.Beep();
    }

    public static void ExampleAll81Orders()
    {
        var c3 = (ConcreteGroup<ZnInt>)(new Cn(3));
        var allOrder9 = BuildExtensions(new[] { (c3, c3) })
            .OrderBy(e => e.Item1.GroupType)
            .ThenByDescending(e => e.Item1.ElementsOrders.Values.Max())
            .ThenBy(e => e.Item2).ToArray();

        var allOrder27 = BuildExtensions(allOrder9.Select(e => (e.Item1, c3)).ToArray())
            .OrderBy(e => e.Item1.GroupType)
            .ThenByDescending(e => e.Item1.ElementsOrders.Values.Max())
            .ThenBy(e => e.Item2).ToArray();

        (allOrder27[2], allOrder27[3]) = (allOrder27[3], allOrder27[2]);
        var allOrder81 = BuildExtensions(allOrder27.Select(e => (e.Item1, c3)).ToArray())
            .Take(15)
            .OrderBy(e => e.Item1.GroupType)
            .ThenByDescending(e => e.Item1.ElementsOrders.Values.Max())
            .ThenBy(e => e.Item2).ToArray();

        Console.Clear();

        DisplayInfosGroups(allOrder9);
        DisplayInfosGroups(allOrder27);
        DisplayInfosGroups(allOrder81);

        Console.Beep(); // ~30min
    }

    public static void ExampleAll36Order()
    {
        var c2 = FG.Abelian(2);
        var c3 = FG.Abelian(3);
        var ab4 = new[] { FG.Abelian(2, 2), FG.Abelian(4) };
        var ab9 = new[] { FG.Abelian(3, 3), FG.Abelian(9) };

        var tuple12 = ab4.SelectMany(e => new[] { (c3, e), (e, c3) }).ToArray();
        var allOrder12 = BuildExtensions(tuple12)
            .OrderBy(e => e.Item1.GroupType)
            .ThenByDescending(e => e.Item1.ElementsOrders.Values.Max())
            .ThenBy(e => e.Item2).ToArray();

        var allOrder18 = BuildExtensions(ab9.ToHashSet())
            .OrderBy(e => e.Item1.GroupType)
            .ThenByDescending(e => e.Item1.ElementsOrders.Values.Max())
            .ThenBy(e => e.Item2).ToArray();

        var ab9a = ab9.Select(e9 => Product.Generate(FG.Abelian(1), e9)).ToList();
        ab9a.ForEach(e => e.Name = ((Gp2<Ep<ZnInt>, Ep<ZnInt>>)e.BaseGroup).G2.Name);
        var tuple36 = ab9a.Select(e9 => (e9, FG.Abelian(4)))
            .Concat(allOrder12.Select(e => (e.Item1, c3)))
            .Concat(allOrder18.Select(e => (e.Item1, c2)))
            .ToArray();

        var allOrder36 = BuildExtensions(tuple36)
            // .Take(14) // stop 
            .OrderBy(e => e.Item1.GroupType)
            .ThenByDescending(e => e.Item1.ElementsOrders.Values.Max())
            .ThenBy(e => e.Item2).ToArray();

        Console.Clear();

        DisplayInfosGroups(allOrder12);
        DisplayInfosGroups(allOrder18);
        DisplayInfosGroups(allOrder36);

        Console.Beep();
    }

    public static void ExampleAll42Order()
    {
        var (c7, c3, c2) = (FG.Abelian(7), FG.Abelian(3), FG.Abelian(2));
        var allOrder14 = BuildExtensions(new[] { (c7, c2) })
            .OrderBy(e => e.Item1.GroupType)
            .ThenByDescending(e => e.Item1.ElementsOrders.Values.Max())
            .ThenBy(e => e.Item2).ToArray();

        var allOrder21 = BuildExtensions(new[] { (c7, c3) })
            .OrderBy(e => e.Item1.GroupType)
            .ThenByDescending(e => e.Item1.ElementsOrders.Values.Max())
            .ThenBy(e => e.Item2).ToArray();

        var tuple42 = allOrder14.Select(e => (e.Item1, c3)).Concat(allOrder21.Select(e => (e.Item1, c2))).ToArray();
        var allOrder42 = BuildExtensions(tuple42)
            .OrderBy(e => e.Item1.GroupType)
            .ThenByDescending(e => e.Item1.ElementsOrders.Values.Max())
            .ThenBy(e => e.Item2).ToArray();

        DisplayInfosGroups(allOrder14);
        DisplayInfosGroups(allOrder21);
        DisplayInfosGroups(allOrder42);
    }

    public static void ExampleAll56Order()
    {
        var allOrder14 = new[] { Group.SemiDirectProd(FG.Abelian(7), FG.Abelian(2)), Product.Generate(FG.Abelian(2), FG.Abelian(7)) }
            .Select(e => (e, SubGroupsDetails(Group.AllSubGroups(e))))
            .OrderBy(e => e.Item1.GroupType)
            .ThenByDescending(e => e.Item1.ElementsOrders.Values.Max())
            .ThenBy(e => e.Item2).ToArray();

        var allOrder28 = BuildExtensions(allOrder14.Select(e => e.Item1).ToHashSet())
            .OrderBy(e => e.Item1.GroupType)
            .ThenByDescending(e => e.Item1.ElementsOrders.Values.Max())
            .ThenBy(e => e.Item2).ToArray();

        var allOrder56 = BuildExtensions(allOrder28.Select(e => e.Item1).ToHashSet())
            .OrderBy(e => e.Item1.GroupType)
            .ThenByDescending(e => e.Item1.ElementsOrders.Values.Max())
            .ThenBy(e => e.Item2).ToArray();

        DisplayInfosGroups(allOrder14);
        DisplayInfosGroups(allOrder28);
        DisplayInfosGroups(allOrder56);

        Console.Beep();
    }

    public static void ExampleAll60Order()
    {
        var ab4 = new[] { Product.Generate(new Cn(2), new Cn(2)), Product.Generate("C4", new Cn(1), new Cn(4)) };
        var c15 = Product.Generate("C15", FG.Abelian(1), Product.Generate(new Cn(3), new Cn(5)));
        var allOrder20 = BuildExtensions(
                ab4.Select(e => (FG.Abelian(5), e))
                    .Append((FG.Abelian(10), Product.Generate("C2", new Cn(1), new Cn(2)))).ToArray())
            .OrderBy(e => e.Item1.GroupType)
            .ThenByDescending(e => e.Item1.ElementsOrders.Values.Max())
            .ThenBy(e => e.Item2).ToArray();

        var tuple6 = new[] { FG.DihedralSdp(3), Product.Generate(new Cn(2), new Cn(3)) };
        var allOrder30 = BuildExtensions(tuple6.Select(e => (FG.Abelian(5), e)).ToArray())
            .OrderBy(e => e.Item1.GroupType)
            .ThenByDescending(e => e.Item1.ElementsOrders.Values.Max())
            .ThenBy(e => e.Item2).ToArray();

        var tuple60 = allOrder20.Select(e => (e.Item1, FG.Abelian(3)))
            .Concat(allOrder30.Select(e => (e.Item1, FG.Abelian(2))))
            .Append((c15, FG.Abelian(4))).ToArray();
        var allOrder60 = BuildExtensions(tuple60)
            .OrderBy(e => e.Item1.GroupType)
            .ThenByDescending(e => e.Item1.ElementsOrders.Values.Max())
            .ThenBy(e => e.Item2).ToArray();

        DisplayInfosGroups(allOrder20);
        DisplayInfosGroups(allOrder30);
        DisplayInfosGroups(allOrder60); // Except A5 which is Simple

        Console.Beep();
    }
}
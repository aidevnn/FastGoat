using FastGoat.Commons;
using FastGoat.Structures.GenericGroup;

namespace FastGoat.Structures;

public delegate T2 GroupAction<in T1, T2>(T1 g, T2 x) where T1 : struct, IElt<T1> where T2 : struct, IElt<T2>;

public static partial class Group
{
    public static GroupAction<T, T> ByTranslate<T>(IGroup<T> gr) where T : struct, IElt<T> => gr.Op;

    public static GroupAction<T, T> ByConjugate<T>(IGroup<T> gr) where T : struct, IElt<T>
    {
        return (T g, T x) => gr.Op(g, gr.Op(x, gr.Invert(g)));
    }

    public static GroupAction<T1, T2> ByAutomorphism<T1, T2>(Homomorphism<T1, Automorphism<T2>> aut)
        where T1 : struct, IElt<T1> where T2 : struct, IElt<T2>
    {
        return (T1 g, T2 x) => aut[g][x];
    }

    public static GroupAction<T, Coset<T>> ByLeftCoset<T>(ConcreteGroup<T> grG, ConcreteGroup<T> grH)
        where T : struct, IElt<T>
    {
        var cosets = Cosets(grG, grH, CosetType.Left);
        return (g, xH) => cosets[grG.Op(g, xH.X)];
    }

    public static HashSet<T1> Stabs<T1, T2>(ConcreteGroup<T1> gr, GroupAction<T1, T2> act, T2 x)
        where T1 : struct, IElt<T1>
        where T2 : struct, IElt<T2>
    {
        return gr.Where(g => act(g, x).Equals(x)).ToHashSet();
    }

    public static HashSet<T2> Orbits<T1, T2>(ConcreteGroup<T1> gr, GroupAction<T1, T2> act, T2 x)
        where T1 : struct, IElt<T1>
        where T2 : struct, IElt<T2>
    {
        HashSet<T2> set = new() { x };
        Queue<T2> q = new Queue<T2>();
        var gens = gr.GetGenerators().ToArray();
        q.Enqueue(x);
        while (q.Count != 0)
        {
            var x0 = q.Dequeue();
            foreach (var g in gens)
            {
                var x1 = act(g, x0);
                if (set.Add(x1))
                    q.Enqueue(x1);
            }
        }

        return set;
    }

    public static Dictionary<T2, (HashSet<T1> Stabx, HashSet<T2> Orbx)> AllOrbits<T1, T2>(ConcreteGroup<T1> gr,
        T2[] set, GroupAction<T1, T2> act)
        where T1 : struct, IElt<T1>
        where T2 : struct, IElt<T2>
    {
        var allSets = new HashSet<HashSet<T2>>(gr.Count(), new SetEquality<T2>());
        var allStabsOrbits = new Dictionary<T2, (HashSet<T1> Stabx, HashSet<T2> Orbx)>();
        foreach (var x in set)
        {
            var orbx = Orbits(gr, act, x);
            if (allSets.Add(orbx))
            {
                var stabx = Stabs(gr, act, x);
                allStabsOrbits[x] = (stabx, orbx);
            }
        }

        return allStabsOrbits;
    }

    public static Dictionary<T, (HashSet<T> Stabx, HashSet<T> Orbx)> AllOrbits<T>(ConcreteGroup<T> gr, GroupAction<T, T> act)
        where T : struct, IElt<T>
    {
        return AllOrbits(gr, gr.ToArray(), act);
    }

    public static (string name, T repr, HashSet<T> stabx, HashSet<T> orbx)[] AllConjugacyClassesNames<T>(ConcreteGroup<T> gr)
        where T : struct, IElt<T>
    {
        var act = ByConjugate(gr);
        var alpha = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";
        var allClasses = AllOrbits(gr, gr.Order().ToArray(), act);
        var classNames = allClasses.GroupBy(e => gr.ElementsOrders[e.Key])
            .ToDictionary(
                e => e.Key,
                e => e.OrderBy(f => f.Value.Orbx.Count)
                    .Select((f, i) => e.Count() == 1 ? ($"{e.Key}", f) : ($"{e.Key}{alpha[i]}", f)).ToArray()
            ).SelectMany(e => e.Value).Select(e => (e.Item1, e.f.Key, e.f.Value.Stabx, e.f.Value.Orbx))
            .OrderBy(e => gr.ElementsOrders[e.Key]).ThenBy(e => e.Orbx.Count).ThenBy(e => e.Item1)
            .ToArray();

        return classNames;
    }

    public static Dictionary<T, (string name, T repr, HashSet<T> stabx, HashSet<T> orbx)>
        AllConjugacyClasses<T>(ConcreteGroup<T> gr)
        where T : struct, IElt<T>
    {
        var clNames = AllConjugacyClassesNames(gr);
        return clNames.SelectMany(e => e.orbx.Select(ei => (ei, e))).ToDictionary(e => e.ei, e => e.e);
    }

    public static void DisplayOrbx<T1, T2>(Dictionary<T2, (HashSet<T1> Stabx, HashSet<T2> Orbx)> allClasses, bool details = false)
        where T1 : struct, IElt<T1>
        where T2 : struct, IElt<T2>
    {
        var i = 0;
        foreach (var kp in allClasses.OrderBy(p => p.Value.Orbx.Count))
        {
            ++i;
            var x = kp.Key;
            var (stabx, orbx) = kp.Value;
            if (details)
                Console.WriteLine($"x{i} = {x,-40} Stab(x{i}):{stabx.Count,-4} Orb(x{i}):{orbx.Count}   {orbx.Glue(", ")}");
            else
                Console.WriteLine($"x{i} = {x,-40} Stab(x{i}):{stabx.Count,-4} Orb(x{i}):{orbx.Count}");
        }

        Console.WriteLine();
    }

    public static void DisplayOrbxSelf<T>(ConcreteGroup<T> gr, bool details = false)
        where T : struct, IElt<T>
    {
        var clNames = AllConjugacyClassesNames(gr);
        var digits = clNames.Max(e => $"{e.repr}".Length);
        var fmt = $"{{0,-{digits}}}";

        foreach (var e in clNames)
        {
            if (details)
                Console.WriteLine(
                    $"{e.name,-3} = {string.Format(fmt, e.repr)} {$"Stab({e.name})",-10}:{e.stabx.Count,-4} {$"Orb({e.name})",-10}:{e.orbx.Count,-4}  {e.orbx.Glue(", ")}");
            else
                Console.WriteLine(
                    $"{e.name,-3} = {string.Format(fmt, e.repr)} {$"Stab({e.name})",-10}:{e.stabx.Count,-4} {$"Orb({e.name})",-10}:{e.orbx.Count,-4}");
        }

        Console.WriteLine($"Nb Classes:{clNames.Length}");
        Console.WriteLine();
    }

    public static void DisplayOrbx<T1, T2>(ConcreteGroup<T1> gr, T2[] set, GroupAction<T1, T2> act, bool details = false)
        where T1 : struct, IElt<T1>
        where T2 : struct, IElt<T2>
    {
        DisplayGroup.Head(gr);
        Console.WriteLine($"Classes for action {act.Method.Name}");
        DisplayOrbx(AllOrbits(gr, set, act), details);
    }

    public static void DisplayOrbx<T>(ConcreteGroup<T> gr, GroupAction<T, T> act)
        where T : struct, IElt<T>
    {
        DisplayOrbx(gr, gr.ToArray(), act);
    }

    public static void DisplayConjugacyClasses<T>(ConcreteGroup<T> gr, bool details = false) where T : struct, IElt<T>
    {
        DisplayGroup.Head(gr);
        Console.WriteLine($"Classes for action {ByConjugate(gr).Method.Name}");

        var clNames = AllConjugacyClassesNames(gr);
        var digits = clNames.Max(e => $"{e.repr}".Length);
        var fmt = $"{{0,-{digits}}}";

        foreach (var e in clNames)
        {
            if (details)
                Console.WriteLine(
                    $"{e.name,-3} = {string.Format(fmt, e.repr)} {$"Stab({e.name})",-10}:{e.stabx.Count,-4} {$"Orb({e.name})",-10}:{e.orbx.Count,-4}  {e.orbx.Glue(", ")}");
            else
                Console.WriteLine(
                    $"{e.name,-3} = {string.Format(fmt, e.repr)} {$"Stab({e.name})",-10}:{e.stabx.Count,-4} {$"Orb({e.name})",-10}:{e.orbx.Count,-4}");
        }

        Console.WriteLine($"Nb Classes:{clNames.Length}");
        Console.WriteLine();
    }

    public static bool AreConjugate<T>(ConcreteGroup<T> g, T a, T b) where T : struct, IElt<T>
    {
        return g.Contains(a) && g.Contains(b) && Orbits(g, ByConjugate(g), a).Contains(b);
    }

    public static ConjugacyClasses<T> ConjugacyClasses<T>(ConcreteGroup<T> gr) where T : struct, IElt<T>
    {
        return new ConjugacyClasses<T>(gr);
    }

    public static List<ConcreteGroup<T>> SubGroupsConjugates<T>(ConcreteGroup<T> g, ConcreteGroup<T> h) where T : struct, IElt<T>
    {
        if (!h.SubSetOf(g))
            throw new GroupException(GroupExceptionType.NotSubGroup);

        var all = new HashSet<ConcreteGroup<T>>(new GroupSetEquality<T>());
        foreach (var s in g)
        {
            var si = g.Invert(s);
            var set = h.GetGenerators().Select(x => g.Op(s, g.Op(x, si))).ToHashSet();
            all.Add(Generate(g, set.ToArray()));
        }

        foreach (var (sg, i) in all.Select((sg, i) => (sg, i)).OrderByDescending(e => e.sg.Count()))
            sg.SetName($"{h.Name}[{i}]");

        return all.ToList();
    }

    public static ConcreteGroup<T1> IsomorphicSubgroup<T1, T2>(ConcreteGroup<T1> g, ConcreteGroup<T2> sg, string name)
        where T1 : struct, IElt<T1> where T2 : struct, IElt<T2>
    {
        var iso = Group.AllMorphisms(sg, g, MorphismType.Isomorphism).FirstOrDefault();
        if (iso.HomMap is null)
            throw new GroupException(GroupExceptionType.GroupDef);

        return Generate($"{name}", g, iso.Image().ToArray());
    }

    public static List<ConcreteGroup<T1>> IsomorphicsSubgroupsAll<T1, T2>(ConcreteGroup<T1> g, ConcreteGroup<T2> sg, string name)
        where T1 : struct, IElt<T1> where T2 : struct, IElt<T2>
    {
        var isos = Group.AllIsomorphisms(sg, g).Select(h => h.Image().ToHashSet()).ToHashSet(new SetEquality<T1>());
        return isos.Select((h, i) => Group.Generate($"{name}[{i + 1}]", g, h.ToArray())).ToList();
    }

    public static List<ConcreteGroup<T1>> IsomorphicsSubgroupsAll<T1, T2>(ConcreteGroup<T1> g, ConcreteGroup<T2> sg)
        where T1 : struct, IElt<T1> where T2 : struct, IElt<T2>
    {
        return IsomorphicsSubgroupsAll(g, sg, sg.Name);
    }

    public static Dictionary<ConcreteGroup<T>, List<ConcreteGroup<T>>> AllSubGroups<T>(ConcreteGroup<T> g)
        where T : struct, IElt<T>
    {
        var og = g.Count();
        var allSubGrs = new HashSet<ConcreteGroup<T>>(20 * og, new GroupSetEquality<T>());
        var table = new Dictionary<ConcreteGroup<T>, List<ConcreteGroup<T>>>(og);
        var count = 0;
        foreach (var e0 in g)
        {
            ++count;
            var cyc = Generate(g, e0);
            if (!allSubGrs.Contains(cyc))
            {
                var conjs = table[cyc] = SubGroupsConjugates(g, cyc);
                allSubGrs.UnionWith(conjs);
            }
        }

        var clsRem = table.Keys.ToHashSet(new GroupSetEquality<T>());
        var cyclesRem = allSubGrs.ToHashSet(new GroupSetEquality<T>());
        while (clsRem.Count != 0)
        {
            var cls = clsRem.OrderBy(sg0 => sg0.Count()).ToArray();
            var cycles = cyclesRem.OrderBy(sg0 => sg0.Count()).ToArray();
            clsRem.Clear();
            cyclesRem.Clear();
            foreach (var (cycle, clRepr) in cycles.Grid2D(cls))
            {
                if (clRepr.SuperSetOf(cycle) || clRepr.SubSetOf(cycle))
                    continue;

                ++count;
                var gens = clRepr.GetGenerators().Union(cycle.GetGenerators()).ToArray();
                var elts = GenerateElements(g, gens);
                if (allSubGrs.All(sg => !sg.SetEquals(elts)))
                {
                    var sg1 = Generate(g, gens);
                    // var sg1 = DirectProduct(clRepr, cycle);
                    var conjs = table[sg1] = SubGroupsConjugates(g, sg1);
                    allSubGrs.UnionWith(conjs);
                    clsRem.Add(sg1);
                    cyclesRem.Add(cycle);
                }
            }
        }

        foreach (var (g0, i) in table.Keys.OrderBy(g0 => g0.Count()).Select((g0, i) => (g0, i)))
            g0.SetName($"SubGr{i + 1}");

        // Console.WriteLine($"{g} NbLoops {count}");
        return table;
    }

    public static List<ConcreteGroup<T>> MaximalSubGroups<T>(List<ConcreteGroup<T>> allSubGr, ConcreteGroup<T> g)
        where T : struct, IElt<T>
    {
        var allMax = new List<ConcreteGroup<T>>();
        foreach (var h in allSubGr)
        {
            if (g.Count() <= h.Count() || !h.SubSetOf(g))
                continue;

            allMax = allMax.Except(allMax.Where(h0 => h0.SubSetOf(h)), new GroupSetEquality<T>()).ToList();
            if (allMax.All(h0 => !h.SubSetOf(h0)))
                allMax.Add(h);
        }

        return allMax;
    }

    public static IEnumerable<List<ConcreteGroup<T>>> SubGroupsLattice<T>(List<ConcreteGroup<T>> allSubGr) where T : struct, IElt<T>
    {
        var g0 = allSubGr.MaxBy(g => g.Count());
        if (g0 is not null)
        {
            var all = new List<List<ConcreteGroup<T>>>() { new() { g0 } };
            while (all.Count != 0)
            {
                var tmp = all.ToList();
                all.Clear();
                foreach (var lt in tmp)
                {
                    var g = lt.Last();
                    if (g.Count() == 1)
                    {
                        yield return lt;
                        continue;
                    }

                    foreach (var h in MaximalSubGroups(allSubGr, g))
                    {
                        var lt2 = lt.Append(h).ToList();
                        all.Add(lt2);
                    }
                }
            }
        }
    }

    public static ConcreteGroup<T> FrattiniSubGroup<T>(List<ConcreteGroup<T>> subGroups, ConcreteGroup<T> g) where T : struct, IElt<T>
    {
        var fratGens = MaximalSubGroups(subGroups, g).Aggregate(g.ToArray(), (acc, g0) => acc.Intersect(g0).ToArray());
        return Generate($"Φ({g})", g, fratGens);
    }

    public static ConcreteGroup<T> FrattiniSubGroup<T>(ConcreteGroup<T> g) where T : struct, IElt<T>
    {
        var subGroups = AllSubGroups(g).Values.SelectMany(e => e).ToList();
        return FrattiniSubGroup(subGroups, g);
    }

    public static ConcreteGroup<T> Commutator<T>(string name, ConcreteGroup<T> grG, ConcreteGroup<T> grH,
        ConcreteGroup<T> grK)
        where T : struct, IElt<T>
    {
        if (!grH.SubSetOf(grG) || !grK.SubSetOf(grG))
            throw new GroupException(GroupExceptionType.NotSubGroup);

        var set = new HashSet<T>();
        foreach (var x in grH)
        {
            var xi = grG.Invert(x);
            foreach (var y in grK)
            {
                var yi = grG.Invert(y);
                var d = grG.Op(grG.Op(x, y), grG.Op(xi, yi));
                set.Add(d);
            }
        }

        var com = Generate(name, grG, set.ToArray());
        return com;
    }

    public static ConcreteGroup<T> Commutator<T>(ConcreteGroup<T> grG, ConcreteGroup<T> grH, ConcreteGroup<T> grK)
        where T : struct, IElt<T>
    {
        return Commutator($"[{grH}, {grK}]", grG, grH, grK);
    }

    static void CommutatorsChain<T>(ConcreteGroup<T> gr, List<ConcreteGroup<T>> chain) where T : struct, IElt<T>
    {
        var i = chain.Count;
        var di = chain.Last();
        var di1 = Commutator($"D{i}", gr, di, gr);
        if (di.SetEquals(di1))
            return;

        chain.Add(di1);
        CommutatorsChain<T>(gr, chain);
    }

    public static List<ConcreteGroup<T>> CommutatorsChain<T>(ConcreteGroup<T> gr) where T : struct, IElt<T>
    {
        var d0 = Create("D0", gr);
        List<ConcreteGroup<T>> chain = new() { d0 };
        CommutatorsChain(gr, chain);
        return chain;
    }

    static void DerivedChain<T>(List<ConcreteGroup<T>> chain) where T : struct, IElt<T>
    {
        var i = chain.Count;
        var gr = chain.Last();
        var comGr = Commutator($"G({i})", gr, gr, gr);
        if (gr.SetEquals(comGr))
            return;

        chain.Add(comGr);
        DerivedChain<T>(chain);
    }

    public static ConcreteGroup<T> Derived<T>(ConcreteGroup<T> gr) where T : struct, IElt<T>
    {
        return Commutator($"D({gr.Name})", gr, gr, gr);
    }

    public static List<ConcreteGroup<T>> DerivedChain<T>(ConcreteGroup<T> gr) where T : struct, IElt<T>
    {
        var g0 = Create("G0", gr);
        List<ConcreteGroup<T>> chain = new() { g0 };
        DerivedChain(chain);
        return chain;
    }

    public static ConcreteGroup<T> Normalize<T>(ConcreteGroup<T> h, ConcreteGroup<T> s) where T : struct, IElt<T>
    {
        var n = h.Where(x => s.SetEquals(s.Select(s0 => h.Op(h.Invert(x), h.Op(s0, x))))).ToArray();
        return Group.Generate($"N[{s}]({h})", h, n);
    }

    public static bool AreConjugate<T>(ConcreteGroup<T> g, ConcreteGroup<T> sg1, ConcreteGroup<T> sg2) where T : struct, IElt<T>
    {
        return sg1.ElementsOrders.Count == sg2.ElementsOrders.Count &&
               g.Any(x => sg1.SetEquals(sg2.Select(s0 => sg1.Op(sg1.Invert(x), sg1.Op(s0, x)))));
    }

    public static ConcreteGroup<T> Zentrum<T>(ConcreteGroup<T> gr) where T : struct, IElt<T>
    {
        var set = new HashSet<T>();
        Func<T, bool> Conj(T s, T si) => x => gr.Op(gr.Op(s, x), si).Equals(x);

        foreach (var s in gr)
        {
            var si = gr.Invert(s);
            var conjByS = Conj(s, si);
            if (gr.All(conjByS))
                set.Add(s);
        }

        return Generate($"Z({gr})", gr, set.ToArray());
    }

    static void ZentrumsChain<T>(ConcreteGroup<T> g, List<ConcreteGroup<T>> chain) where T : struct, IElt<T>
    {
        var i = chain.Count;
        var zi = chain.Last();

        var quo = g.Over(zi);
        var zQuo = Zentrum(quo);
        var naturalMap = AllHomomorphisms(g, quo).First(h => zi.SetEquals(h.Kernel()));
        var preImag = naturalMap.HomMap.Where(kp => zQuo.Contains(kp.Value)).Select(kp => kp.Key).ToArray();
        var zi1 = Generate($"Z{i}", g, preImag);

        if (zi1.SetEquals(zi))
            return;

        chain.Add(zi1);
        ZentrumsChain(g, chain);
    }

    // H.E. Rose. A Course on Finite Groups. Problem 10.1 page 223
    static void ZentrumsChainFast<T>(ConcreteGroup<T> gr, List<ConcreteGroup<T>> chain) where T : struct, IElt<T>
    {
        var i = chain.Count;
        var zi = chain.Last();

        T Commute(T x, T y) => gr.Op(gr.Op(x, y), gr.Op(gr.Invert(x), gr.Invert(y)));

        HashSet<T> set = new(gr.Count());
        foreach (var g in gr)
        {
            if (gr.All(a => zi.Contains(Commute(a, g))))
                set.Add(g);
        }

        if (zi.SetEquals(set))
            return;

        var zi1 = Generate($"Z{i}", gr, set.ToArray());
        chain.Add(zi1);
        ZentrumsChainFast(gr, chain);
    }

    public static List<ConcreteGroup<T>> ZentrumsChain<T>(ConcreteGroup<T> g) where T : struct, IElt<T>
    {
        var z0 = Group.Generate("Z0", g, g.Neutral());
        List<ConcreteGroup<T>> chain = new() { z0 };
        ZentrumsChain(g, chain);
        return chain;
    }

    public static List<ConcreteGroup<T>> ZentrumsChainFast<T>(ConcreteGroup<T> g) where T : struct, IElt<T>
    {
        var z0 = Group.Generate("Z0", g, g.Neutral());
        List<ConcreteGroup<T>> chain = new() { z0 };
        ZentrumsChainFast(g, chain);
        return chain;
    }
}
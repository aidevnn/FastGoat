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
        return (T g, Coset<T> xH) => cosets[grG.Op(g, xH.X)];
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
        q.Enqueue(x);
        while (q.Count != 0)
        {
            var x0 = q.Dequeue();
            foreach (var g in gr.PseudoGenerators)
            {
                var x1 = act(g, x0);
                if (set.Add(x1))
                    q.Enqueue(x1);
            }
        }

        // return gr.Select(g => act(g, x)).ToHashSet();
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

    public static Dictionary<T, (HashSet<T> Stabx, HashSet<T> Orbx)> AllOrbits<T>(ConcreteGroup<T> gr,
        GroupAction<T, T> act)
        where T : struct, IElt<T>
    {
        return AllOrbits(gr, gr.ToArray(), act);
    }

    public static void DisplayOrbx<T1, T2>(Dictionary<T2, (HashSet<T1> Stabx, HashSet<T2> Orbx)> allClasses)
        where T1 : struct, IElt<T1>
        where T2 : struct, IElt<T2>
    {
        foreach (var kp in allClasses.OrderBy(p => p.Value.Orbx.Count))
        {
            var x = kp.Key;
            var (stabx, orbx) = kp.Value;
            Console.WriteLine($"x={x,-40} Stab(x):{stabx.Count,-4} Orb(x):{orbx.Count} {orbx.Glue(", ")}");
        }
    }

    public static void DisplayOrbx<T1, T2>(ConcreteGroup<T1> gr,
        T2[] set, GroupAction<T1, T2> act)
        where T1 : struct, IElt<T1>
        where T2 : struct, IElt<T2>
    {
        DisplayOrbx(AllOrbits(gr, set, act));
    }

    public static void DisplayOrbx<T>(ConcreteGroup<T> gr, GroupAction<T, T> act)
        where T : struct, IElt<T>
    {
        DisplayOrbx(gr, gr.ToArray(), act);
    }

    public static bool AreConjugate<T>(ConcreteGroup<T> g, T a, T b) where T : struct, IElt<T>
    {
        return g.Contains(a) && g.Contains(b) && Orbits(g, ByConjugate(g), a).Contains(b);
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

    public static List<ConcreteGroup<T>> DerivedChain<T>(ConcreteGroup<T> gr) where T : struct, IElt<T>
    {
        var g0 = Create("G0", gr);
        List<ConcreteGroup<T>> chain = new() { g0 };
        DerivedChain(chain);
        return chain;
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
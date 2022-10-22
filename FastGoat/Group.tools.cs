using FastGoat.UserGroup;

namespace FastGoat;

public delegate T2 GroupAction<in T1, T2>(T1 g, T2 x) where T1 : struct, IElt<T1> where T2 : struct, IElt<T2>;

public static partial class Group
{
    public static GroupAction<T, T> ByTranslate<T>(IGroup<T> gr) where T : struct, IElt<T> => gr.Op;

    public static GroupAction<T, T> ByConjugate<T>(IGroup<T> gr) where T : struct, IElt<T>
    {
        return (T g, T x) => gr.Op(g, gr.Op(x, gr.Invert(g)));
    }

    public static GroupAction<T1, T2> ByAutomorphism<T1, T2>(IDictionary<T1, IDictionary<T2, T2>> aut)
        where T1 : struct, IElt<T1> where T2 : struct, IElt<T2>
    {
        return (T1 g, T2 x) => aut[g][x];
    }

    public static GroupAction<T1, T2> ByAutomorphism<T1, T2>(IDictionary<T1, Automorphism<T2>> aut)
        where T1 : struct, IElt<T1> where T2 : struct, IElt<T2>
    {
        return (T1 g, T2 x) => aut[g][x];
    }

    public static GroupAction<T, T> ByLeftCoset<T>(ConcreteGroup<T> grG, ConcreteGroup<T> grH) where T : struct, IElt<T>
    {
        var coset = Cosets(grG, grH, Coset.Left);
        return (T g, T x) => coset[grG.Op(g, x)];
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

    public static ConcreteGroup<T> Commutator<T>(ConcreteGroup<T> gr) where T : struct, IElt<T>
    {
        var set = new HashSet<T>();
        foreach (var x in gr)
        {
            var xi = gr.Invert(x);
            foreach (var y in gr)
            {
                var yi = gr.Invert(y);
                var d = gr.Op(gr.Op(x, y), gr.Op(xi, yi));
                set.Add(d);
            }
        }

        var com = Group.Generate($"D({gr})", gr, set.ToArray());
        return com;
    }

    static void CommutatorsChain<T>(List<ConcreteGroup<T>> chain) where T : struct, IElt<T>
    {
        var gr = chain.Last();
        var nb0 = gr.Count();
        if (nb0 == 1)
            return;

        var comGr = Commutator(gr);
        var nb1 = comGr.Count();
        if (nb0 == nb1)
            return;

        chain.Add(comGr);
        CommutatorsChain<T>(chain);
    }

    public static List<ConcreteGroup<T>> CommutatorsChain<T>(ConcreteGroup<T> gr) where T : struct, IElt<T>
    {
        List<ConcreteGroup<T>> chain = new() { gr };
        CommutatorsChain(chain);
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

}
using System.Collections.ObjectModel;

namespace FastGoat;

public static partial class Group
{
    public static T Times<T>(this IGroup<T> g, T e, int p) where T : struct, IElt<T>
    {
        var acc = g.Neutral();
        if (p == 0)
            return acc;

        var e0 = p > 0 ? e : g.Invert(e);
        var p0 = p > 0 ? p : -p;
        for (var k = 0; k < p0; ++k)
            acc = g.Op(acc, e0);

        return acc;
    }

    private static IEnumerable<T> GenerateElements<T>(IGroup<T> bg, IEnumerable<T> elements) where T : struct, IElt<T>
    {
        var generators = elements.ToHashSet().ToArray();
        var bg0 = bg.Neutral().BaseGroup;
        if (generators.Any(e => !bg0.Equals(e.BaseGroup)))
            throw new GroupException(GroupExceptionType.BaseGroup);

        var q = new Queue<T>();
        q.Enqueue(bg.Neutral());
        HashSet<T> generatedElements = new HashSet<T>() { bg.Neutral() };
        while (q.Count != 0)
        {
            var e1 = q.Dequeue();
            foreach (var e2 in generators)
            {
                var e3 = bg.Op(e1, e2);
                if (generatedElements.Add(e3))
                    q.Enqueue(e3);
            }
        }

        return generatedElements;
    }

    public static IEnumerable<T> GenerateElements<T>(IGroup<T> bg, params T[] elements) where T : struct, IElt<T>
    {
        return GenerateElements(bg, elements.ToList());
    }

    public static ReadOnlyDictionary<T, int> Cycle<T>(IGroup<T> g, T e) where T : struct, IElt<T>
    {
        var e0 = e;
        var n = g.Neutral();
        var p = 1;
        var elements = new Dictionary<T, int> { [e] = p };
        while (!n.Equals(e0))
        {
            e0 = g.Op(e0, e);
            elements[e0] = ++p;
        }

        return new ReadOnlyDictionary<T, int>(elements);
    }

    public static ReadOnlyDictionary<T, ReadOnlyDictionary<T, int>> LongestCycles<T>(IGroup<T> g,
        IEnumerable<T> elements)
        where T : struct, IElt<T>
    {
        var set = elements.ToHashSet();
        var allCycles = new Dictionary<T, ReadOnlyDictionary<T, int>>(set.Count);

        while (set.Count != 0)
        {
            var e0 = set.First();
            var cycle0 = Cycle(g, e0);
            set.ExceptWith(cycle0.Keys);
            if (allCycles.Count == 0)
            {
                allCycles[e0] = cycle0;
                continue;
            }

            var tmpCycles = allCycles.Select(p => (p.Key, p.Value)).ToArray();
            allCycles.Clear();
            var done = false;
            foreach (var (e1, cycle1) in tmpCycles)
                if (cycle0.Count % cycle1.Count != 0)
                {
                    allCycles[e1] = cycle1;
                }
                else
                {
                    if (cycle1.ContainsKey(e0))
                    {
                        allCycles[e1] = cycle1;
                        done = true;
                    }
                    else if (!cycle0.ContainsKey(e1))
                    {
                        allCycles[e1] = cycle1;
                    }
                }

            if (!done) allCycles[e0] = cycle0;
        }

        return new ReadOnlyDictionary<T, ReadOnlyDictionary<T, int>>(allCycles);
    }

    public static ReadOnlyDictionary<T, int> ElementsOrders<T>(
        ReadOnlyDictionary<T, ReadOnlyDictionary<T, int>> longCycles)
        where T : struct, IElt<T>
    {
        var orders = new Dictionary<T, int>();
        foreach (var (_, cycle) in longCycles)
        {
            var n = cycle.Count;
            foreach (var (e1, p) in cycle)
            {
                var o = n / IntExt.Gcd(n, p);
                if (!orders.ContainsKey(e1))
                    orders[e1] = o;
                else if (orders[e1] != o)
                    throw new Exception("TODO"); // TODO unexpected case
            }
        }

        return new ReadOnlyDictionary<T, int>(orders);
    }

    public static bool IsCommutative<T>(IGroup<T> g, IEnumerable<T> ts) where T : struct, IElt<T>
    {
        var ts0 = ts.ToArray();
        foreach (var e1 in ts0)
        foreach (var e2 in ts0)
        {
            var e12 = g.Op(e1, e2);
            var e21 = g.Op(e2, e1);
            if (!e12.Equals(e21))
                return false;
        }

        return true;
    }

    public static Dictionary<T, ReadOnlyCollection<T>> Cosets<T>(IGroup<T> g, IGroup<T> n) where T : struct, IElt<T>
    {
        Dictionary<T, ReadOnlyCollection<T>> cosets = new();
        List<HashSet<T>> sets = new();
        var setH = n.ToHashSet();
        if (!setH.IsSubsetOf(g))
            throw new GroupException(GroupExceptionType.NotSubGroup);

        cosets[g.Neutral()] = new ReadOnlyCollection<T>(setH.ToList());
        sets.Add(setH);
        var xH = new HashSet<T>(setH.Count);
        foreach (var x in g)
        {
            if (x.Equals(g.Neutral()))
                continue;

            xH.Clear();
            var xi = g.Invert(x);
            xH.UnionWith(setH.Select(h => g.Op(x, h)));
            if (xH.Any(xh => !setH.Contains(g.Op(xh, xi))))
                throw new GroupException(GroupExceptionType.NotNormal);

            var x0 = xH.First();
            if (sets.All(set => !set.Contains(x0)))
            {
                sets.Add(xH.ToHashSet());
                cosets[x0] = new ReadOnlyCollection<T>(xH.ToList());
            }
        }

        return cosets;
    }

    public static ConcreteGroup<T> Create<T>(string name, IGroup<T> g) where T : struct, IElt<T>
    {
        return new ConcreteGroup<T>(name, g);
    }

    public static ConcreteGroup<T> Create<T>(IGroup<T> g) where T : struct, IElt<T>
    {
        return new ConcreteGroup<T>(g);
    }

    public static ConcreteGroup<T> Generate<T>(T e, params T[] others) where T : struct, IElt<T>
    {
        return Generate(e.BaseGroup.Name, e, others);
    }

    public static ConcreteGroup<T> Generate<T>(string name, T e, params T[] others) where T : struct, IElt<T>
    {
        return new ConcreteGroup<T>(name, e.BaseGroup, others.Prepend(e).ToArray());
    }

    public static ConcreteGroup<T> Generate<T>(IGroup<T> g, params T[] generators) where T : struct, IElt<T>
    {
        return new ConcreteGroup<T>(g, generators);
    }

    public static ConcreteGroup<T> Generate<T>(string name, IGroup<T> g, params T[] generators)
        where T : struct, IElt<T>
    {
        return new ConcreteGroup<T>(name, g, generators);
    }

    public static ConcreteGroup<T> DirectProduct<T>(ConcreteGroup<T> g1, ConcreteGroup<T> g2) where T : struct, IElt<T>
    {
        if (!g1.BaseGroup.Equals(g2.BaseGroup))
            throw new GroupException(GroupExceptionType.GroupDef);

        if (g1.SuperGroup is null || g2.SuperGroup is null || !g1.SuperGroup.Equals(g2.SuperGroup))
            throw new GroupException(GroupExceptionType.GroupDef);

        var generators = g1.LongestCycles.Keys.Concat(g2.LongestCycles.Keys).Distinct().ToArray();
        return new ConcreteGroup<T>(g1.SuperGroup, generators);
    }

    public static QuotientGroup<T> Over<T>(this ConcreteGroup<T> g, ConcreteGroup<T> h) where T : struct, IElt<T>
    {
        var gName = g.Name.Trim().Contains(' ') ? $"({g.Name})" : g.Name;
        var hName = h.Name.Trim().Contains(' ') ? $"({h.Name})" : h.Name;
        return new QuotientGroup<T>(g, h, $"{gName}/{hName}");
    }

    public static QuotientGroup<T> Over<T>(this ConcreteGroup<T> g, ConcreteGroup<T> h, string name)
        where T : struct, IElt<T>
    {
        return new QuotientGroup<T>(g, h, name);
    }

    public static SemiDirectProduct<T1, T2> SemiDirectProd<T1, T2>(string name, ConcreteGroup<T1> n,
        ConcreteGroup<T2> g)
        where T1 : struct, IElt<T1> where T2 : struct, IElt<T2>
    {
        return new SemiDirectProduct<T1, T2>(name, n, g);
    }

    public static SemiDirectProduct<T1, T2> SemiDirectProd<T1, T2>(ConcreteGroup<T1> n,
        ConcreteGroup<T2> g)
        where T1 : struct, IElt<T1> where T2 : struct, IElt<T2>
    {
        return new SemiDirectProduct<T1, T2>(n, g);
    }
}
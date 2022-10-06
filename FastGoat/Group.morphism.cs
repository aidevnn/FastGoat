using FastGoat.Examples;
using FastGoat.Gp;
using FastGoat.UserGroup;

namespace FastGoat;

public static partial class Group
{
    public static Dictionary<T1, T2> PartialMap<T1, T2>(params (T1, T2)[] elements)
        where T1 : struct, IElt<T1>
        where T2 : struct, IElt<T2>
    {
        return elements.ToDictionary(p => p.Item1, p => p.Item2);
    }

    public static Dictionary<T1, T2> HomomorphismMap<T1, T2>(ConcreteGroup<T1> g, ConcreteGroup<T2> h,
        Dictionary<T1, T2> pMap)
        where T1 : struct, IElt<T1>
        where T2 : struct, IElt<T2>
    {
        Dictionary<T1, T2> map = new();
        if (pMap.ContainsKey(g.Neutral()) && !pMap[g.Neutral()].Equals(h.Neutral()))
            return map;

        if (pMap.Any(kp => !g.Contains(kp.Key) || !h.Contains(kp.Value)))
            return map;

        map[g.Neutral()] = h.Neutral();
        var q = new Queue<T1>();
        q.Enqueue(g.Neutral());
        while (q.Count != 0)
        {
            var e1 = q.Dequeue();
            var f1 = map[e1];
            foreach (var kp in pMap)
            {
                var e2 = kp.Key;
                var f2 = kp.Value;
                var e3 = g.Op(e1, e2);
                var f3 = h.Op(f1, f2);
                if (!map.ContainsKey(e3))
                {
                    map[e3] = f3;
                    q.Enqueue(e3);
                }
                else if (!map[e3].Equals(f3))
                {
                    return new();
                }
            }
        }

        return map;
    }

    public static Dictionary<T, T> EndomorphismMap<T>(ConcreteGroup<T> g, Dictionary<T, T> pMap)
        where T : struct, IElt<T>
    {
        return HomomorphismMap(g, g, pMap);
    }

    public static Dictionary<T1, T2> IsomorphismMap<T1, T2>(ConcreteGroup<T1> g, ConcreteGroup<T2> h,
        Dictionary<T1, T2> pMap)
        where T1 : struct, IElt<T1>
        where T2 : struct, IElt<T2>
    {
        if (pMap.Any(kp => g.ElementsOrders[kp.Key] != h.ElementsOrders[kp.Value]))
            return new();

        var map = HomomorphismMap(g, h, pMap);
        if (map.Any(kp => g.ElementsOrders[kp.Key] != h.ElementsOrders[kp.Value]))
            return new();

        return map;
    }

    public static Dictionary<T, T> AutomorphismMap<T>(ConcreteGroup<T> g, Dictionary<T, T> pMap)
        where T : struct, IElt<T>
    {
        return IsomorphismMap(g, g, pMap);
    }

    public static AutomorphismGroup<T> Aut<T>(ConcreteGroup<T> g) where T : struct, IElt<T>
    {
        return new AutomorphismGroup<T>(g);
    }

    public static ConcreteGroup<Automorphism<T>> Aut<T>(string name, IGroup<T> bg, T[] generators)
        where T : struct, IElt<T>
    {
        var g = Generate(name, bg, generators);
        var bgAut = Aut(g);
        var gByOrders = g.ElementsOrders.GroupBy(kp => kp.Value)
            .ToDictionary(a => a.Key, a => a.Select(kp => kp.Key).ToArray());
        var gensPossibles = generators.ToDictionary(a => a, a => gByOrders[g.ElementsOrders[a]]);
        var allTuples = gensPossibles.Select(kp => kp.Value.Select(a => (kp.Key, a))).ToArray();
        List<Dictionary<T, T>> allMaps = new();
        foreach (var tuples in allTuples)
        {
            if (allMaps.Count == 0)
            {
                allMaps = tuples.Select(t => new Dictionary<T, T>() { [t.Key] = t.a }).ToList();
                continue;
            }

            List<Dictionary<T, T>> tmpMaps = new();
            foreach (var t in tuples)
            {
                foreach (var map in allMaps)
                {
                    var map0 = new Dictionary<T, T>(map) { [t.Key] = t.a };
                    tmpMaps.Add(map0);
                }
            }

            allMaps.Clear();
            allMaps.AddRange(tmpMaps);
        }

        var allAuts = allMaps.Select(pMap => AutomorphismMap(g, pMap)).Where(map => map.Count == g.Count())
            .Select(map => new Automorphism<T>(bgAut, map))
            .ToHashSet();

        var autG = Generate($"Aut[{name}]", allAuts.First(), allAuts.Skip(1).ToArray());
        return autG;
    }

    public static ConcreteGroup<Automorphism<T>> Aut<T>(T e, params T[] others) where T : struct, IElt<T>
    {
        return Aut(e.BaseGroup.Name, e, others);
    }

    public static ConcreteGroup<Automorphism<T>> Aut<T>(string name, T e, params T[] others) where T : struct, IElt<T>
    {
        return Aut(name, e.BaseGroup, others.Prepend(e).ToArray());
    }

    public static ConcreteGroup<Automorphism<T>> Aut<T>(IGroup<T> bg, T[] generators)
        where T : struct, IElt<T>
    {
        return Aut(bg.Name, bg, generators);
    }
}
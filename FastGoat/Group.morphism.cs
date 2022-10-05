using FastGoat.Examples;
using FastGoat.Gp;

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

        return HomomorphismMap(g, h, pMap);
    }

    public static Dictionary<T, T> AutomorphismMap<T>(ConcreteGroup<T> g, Dictionary<T, T> pMap)
        where T : struct, IElt<T>
    {
        return IsomorphismMap(g, g, pMap);
    }
}
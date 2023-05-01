using FastGoat.Commons;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.VecSpace;

namespace FastGoat.Structures;

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
        var gCount = g.Count();
        Dictionary<T1, T2> map = new(gCount);
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

    public static AutomorphismGroup<T> AutBase<T>(ConcreteGroup<T> g) where T : struct, IElt<T>
    {
        return new AutomorphismGroup<T>(g);
    }

    public enum MorphismType
    {
        Homomorphism,
        Isomorphism
    }

    public static IEnumerable<Homomorphism<T1, T2>> AllMorphisms<T1, T2>(IGroup<T1> bg, ConcreteGroup<T2> g2,
        MorphismType mType = MorphismType.Homomorphism)
        where T1 : struct, IElt<T1>
        where T2 : struct, IElt<T2>
    {
        var gGens = bg.GetGenerators().ToArray();
        ConcreteGroup<T1> g1;
        if (bg is ConcreteGroup<T1> gi)
            g1 = gi;
        else
            g1 = Group.Generate(bg, gGens);

        bool Filter(int e, int a) => mType == MorphismType.Homomorphism ? e % a == 0 : e == a;

        var g2ByOrders = g2.GroupBy(e => g2.ElementsOrders[e])
            .Select(e => (ord: e.Key, elt: e.ToArray()))
            .ToArray();
        var gGensOrders = gGens.Select(e => (g: e, ord: g1.ElementsOrders[e])).ToArray();
        var gpMap = gGensOrders.Select(e =>
                g2ByOrders.Where(a => Filter(e.ord, a.ord)).SelectMany(a => a.elt).Select(a => (e.g, a))
                    .ToArray())
            .ToArray();

        var ng = g1.Count();
        // Console.WriteLine(gpMap.Aggregate(1, (acc, b) => acc * b.Length));
        foreach (var arr in gpMap.MultiLoop())
        {
            var map = arr.ToDictionary(t => t.g, t => t.a);
            if (mType == MorphismType.Homomorphism)
            {
                var hom = Group.HomomorphismMap(g1, g2, map);
                if (hom.Count == ng)
                    yield return new(g1, hom);
            }
            else
            {
                var iso = Group.IsomorphismMap(g1, g2, map);
                if (iso.Count == ng && iso.Values.Count() == ng)
                    yield return new(g1, iso);
            }
        }
    }

    public static List<Homomorphism<T1, T2>> AllIsomorphisms<T1, T2>(IGroup<T1> bg, ConcreteGroup<T2> g2)
        where T1 : struct, IElt<T1>
        where T2 : struct, IElt<T2>
    {
        return AllMorphisms(bg, g2, MorphismType.Isomorphism).ToList();
    }

    public static List<Homomorphism<T, T>> AllAutomorphisms<T>(ConcreteGroup<T> g)
        where T : struct, IElt<T>
    {
        return AllIsomorphisms(g, g);
    }

    public static ConcreteGroup<Automorphism<T>> AutomorphismGroup<T>(ConcreteGroup<T> g)
        where T : struct, IElt<T>
    {
        var bgAut = new AutomorphismGroup<T>(g);
        var allAut = AllAutomorphisms(g);
        var autG = Generate($"Aut[{g.Name}]", bgAut,
            allAut.Select(aut => new Automorphism<T>(bgAut, aut.HomMap)).ToArray());
        return autG;
    }

    public static List<Homomorphism<T1, T2>> AllHomomorphisms<T1, T2>(IGroup<T1> bg, ConcreteGroup<T2> g2)
        where T1 : struct, IElt<T1>
        where T2 : struct, IElt<T2>
    {
        return AllMorphisms(bg, g2).ToList();
    }

    public static List<Homomorphism<T2, Automorphism<T1>>> AllOpsByAutomorphisms<T1, T2>(IGroup<T2> bg,
        ConcreteGroup<T1> bn)
        where T1 : struct, IElt<T1>
        where T2 : struct, IElt<T2>
    {
        var nGens = bn.GetGenerators().ToArray();
        var autN = Group.AutomorphismGroup(bn);
        return AllHomomorphisms(bg, autN);
    }

    public static Homomorphism<T1, T2> Hom<T1, T2>(ConcreteGroup<T1> g, Dictionary<T1, T2> map)
        where T1 : struct, IElt<T1>
        where T2 : struct, IElt<T2>
    {
        if (map.Count != g.Count())
            throw new GroupException(GroupExceptionType.GroupDef);

        return new(g, map);
    }

    public static HashSet<T2> Image<T1, T2>(Dictionary<T1, T2> map)
        where T1 : struct, IElt<T1>
        where T2 : struct, IElt<T2>
    {
        return map.Values.ToHashSet();
    }

    public static HashSet<T1> Kernel<T1, T2>(T2 neutral, Dictionary<T1, T2> map)
        where T1 : struct, IElt<T1>
        where T2 : struct, IElt<T2>
    {
        return map.Where(kp => kp.Value.Equals(neutral)).Select(kp => kp.Key).ToHashSet();
    }

    public static Dictionary<T1, T3> Compose<T1, T2, T3>(Dictionary<T1, T2> map1, Dictionary<T2, T3> map2)
        where T1 : struct, IElt<T1>
        where T2 : struct, IElt<T2>
    {
        return map1.Where(kp => map2.ContainsKey(kp.Value)).ToDictionary(a => a.Key, a => map2[a.Value]);
    }

    public static SemiDirectProduct<T1, T2> SemiDirectProd<T1, T2>(string name, ConcreteGroup<T1> n, ConcreteGroup<T2> g)
        where T1 : struct, IElt<T1> where T2 : struct, IElt<T2>
    {
        var allOps = AllOpsByAutomorphisms(g, n);
        var idx = allOps.FindIndex(kp => kp.Image().Count() > 1 && kp.Kernel().Count() == 1); // First faithfull action
        if (idx == -1)
        {
            idx = allOps.FindIndex(kp => kp.Image().Count() > 1); // First non faithfull action
            if (idx == -1)
                throw new GroupException(GroupExceptionType.SemiDirectProductDontExist);
        }

        var theta = allOps[idx];
        return new SemiDirectProduct<T1, T2>(name, n, theta, g);
    }

    public static SemiDirectProduct<T1, T2> SemiDirectProd<T1, T2>(ConcreteGroup<T1> n, ConcreteGroup<T2> g)
        where T1 : struct, IElt<T1> where T2 : struct, IElt<T2>
    {
        return SemiDirectProd("", n, g);
    }

    public static SemiDirectProduct<T1, T2> SemiDirectProd<T1, T2>(ConcreteGroup<T1> n, Homomorphism<T2, Automorphism<T1>> theta, ConcreteGroup<T2> g)
        where T1 : struct, IElt<T1> where T2 : struct, IElt<T2>
    {
        return new SemiDirectProduct<T1, T2>("", n, theta, g);
    }

    public static SemiDirectProduct<T1, T2> SemiDirectProd<T1, T2>(string name, ConcreteGroup<T1> n,
        Homomorphism<T2, Automorphism<T1>> theta, ConcreteGroup<T2> g)
        where T1 : struct, IElt<T1> where T2 : struct, IElt<T2>
    {
        return new SemiDirectProduct<T1, T2>(name, n, theta, g);
    }

    public static ConcreteGroup<EPoly<K>> KAut<K>(KPoly<K> poly) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        return new ConcreteGroup<EPoly<K>>(new KAut<K>(poly));
    }
    public static ConcreteGroup<EPoly<K>> KAut<K>(params EPoly<K>[] gens) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        if (gens.Length == 0 || gens.Select(g => g.F).Distinct().Count() != 1)
            throw new GroupException(GroupExceptionType.GroupDef);
        
        var poly = gens[0].F;
        return new ConcreteGroup<EPoly<K>>(new KAut<K>(poly), gens);
    }
}
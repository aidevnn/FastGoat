using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.GenericGroup;

namespace Craft.Craft;

public static partial class GroupCraft
{
    static void HomomorphismMapInplace<T1, T2>(IGroup<T1> g, IGroup<T2> h, Dictionary<T1, T2> pMap, Dictionary<T1, T2> map)
        where T1 : struct, IElt<T1>
        where T2 : struct, IElt<T2>
    {
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
                    map.Clear();
                    return;
                }
            }
        }
    }

    public static IEnumerable<Dictionary<T1, T2>> MorphismPruning<T1, T2>(IGroup<T1> G, IGroup<T2> H,
        Dictionary<T1, int> sizes,
        Dictionary<T1, T2> pMap, Dictionary<T1, T2> map, (T1 g, T2 a)[][] gens, int idx, Group.MorphismType mType)
        where T1 : struct, IElt<T1>
        where T2 : struct, IElt<T2>
    {
        foreach (var (g, a) in gens[idx])
        {
            var map2 = map.ToDictionary();
            var pMap2 = pMap.ToDictionary();
            pMap2[g] = a;
            if (mType == Group.MorphismType.Isomorphism && pMap2.Keys.Count != pMap2.Values.Distinct().Count())
                continue; 
            HomomorphismMapInplace(G, H, pMap2, map2);
            if (map2.Count == sizes[g])
            {
                if (idx == gens.Length - 1)
                    yield return map2;
                else
                {
                    foreach (var hom in MorphismPruning(G, H, sizes, pMap2, map2, gens, idx + 1, mType))
                        yield return hom;
                }
            }
        }
    }

    public static IEnumerable<Homomorphism<T1, T2>> AllMorphismsWithPruning<T1, T2>(ConcreteGroup<T1> G,
        ConcreteGroup<T2> H, Dictionary<T1,T2> mapInit, Group.MorphismType mType)
        where T1 : struct, IElt<T1> where T2 : struct, IElt<T2>
    {
        var gGens = G.GetGenerators().ToArray();
        bool Filter(int e, int a) => mType == Group.MorphismType.Homomorphism ? e % a == 0 : e == a;
        var hByOrders = H.GroupBy(e => H.ElementsOrders[e])
            .Select(e => (ord: e.Key, elt: e.ToArray()))
            .ToArray();
        var gGensOrders = gGens.Select(e => (g: e, ord: G.ElementsOrders[e])).ToArray();
        var gpMap = gGensOrders.Select(e => hByOrders.Where(a => Filter(e.ord, a.ord)).SelectMany(a => a.elt)
                .Select(a => (e.g, a))
                .ToArray())
            .ToArray();

        var ng = G.Order;
        var sizes = (gGens.Length - 1).SeqLazy(1)
            .ToDictionary(i => gGens[i - 1], i => Group.GenerateElements(G, gGens.Take(i).ToArray()).Count);
        sizes[gGens.Last()] = ng;
        foreach (var hom in MorphismPruning(G, H, sizes, new(), mapInit, gpMap, 0, mType))
        {
            if (mType == Group.MorphismType.Homomorphism)
                yield return new(G, hom);
            else if (hom.Values.Distinct().Count() == ng)
                yield return new(G, hom);
        }
    }


    public static IEnumerable<Homomorphism<T1, T2>> AllMorphismsWithPruning<T1, T2>(ConcreteGroup<T1> G,
        ConcreteGroup<T2> H, Group.MorphismType mType = Group.MorphismType.Homomorphism)
        where T1 : struct, IElt<T1> where T2 : struct, IElt<T2>
    {
        return AllMorphismsWithPruning(G, H, new() { [G.Neutral()] = H.Neutral() }, mType);
    }

    public static bool AreIsomorphic<T1, T2>(ConcreteGroup<T1> G, ConcreteGroup<T2> H) where T1 : struct, IElt<T1>
        where T2 : struct, IElt<T2>
    {
        if (G.Order != H.Order || G.GroupType != H.GroupType)
            return false;
        if (!G.ElementsOrdersList().SequenceEqual(H.ElementsOrdersList()))
            return false;

        Homomorphism<T1, T2> nullHom = new();
        var iso = AllMorphismsWithPruning(G, H, Group.MorphismType.Isomorphism).FirstOrDefault(nullHom);
        return !iso.IsNull;
    }

    public static void DisplayAreIsomorphic<T1, T2>(ConcreteGroup<T1> G, ConcreteGroup<T2> H) where T1 : struct, IElt<T1>
        where T2 : struct, IElt<T2>
    {
        var test = AreIsomorphic(G, H);
        Console.WriteLine($"{G} IsIsomorphic to {H}:{test}");
        Console.WriteLine();
    }

    public static ConcreteGroup<T1> IsomorphicSubgroup<T1, T2>(ConcreteGroup<T1> G, ConcreteGroup<T2> H)
        where T1 : struct, IElt<T1>
        where T2 : struct, IElt<T2>
    {
        Homomorphism<T2, T1> nullHom = new();
        var iso = AllMorphismsWithPruning(H, G, Group.MorphismType.Isomorphism).FirstOrDefault(nullHom);
        if (iso.IsNull)
            return G.C1;
        return Group.Generate(H.Name, G, H.GetGenerators().Select(e => iso[e]).ToArray());
    }

    public static IEnumerable<ConcreteGroup<T1>> AllIsomorphicSubGroups<T1, T2>(ConcreteGroup<T1> G, ConcreteGroup<T2> H)
        where T1 : struct, IElt<T1>
        where T2 : struct, IElt<T2>
    {
        foreach (var iso in AllMorphismsWithPruning(H, G, Group.MorphismType.Isomorphism))
            yield return Group.Generate(H.Name, G, H.GetGenerators().Select(e => iso[e]).ToArray());
    }
}
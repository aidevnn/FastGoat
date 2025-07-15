using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;

namespace FastGoat.UserGroup.GModuleN;

public static class TC
{
    public static IOrderedEnumerable<KeyValuePair<Ep<Tg>, Tn>> OrderKeys<Tn, Tg>(this IEnumerable<KeyValuePair<Ep<Tg>, Tn>> kvs,
        ConcreteGroup<Tg> G)
        where Tg : struct, IElt<Tg>
        where Tn : struct, IElt<Tn>
    {
        return kvs.OrderByDescending(kv => kv.Key.Ei.Any(e => e.Equals(G.Neutral())) ? 1 : 0).ThenBy(kv => kv.Key);
    }

    public static MapElt<Tg, Tn> ToMapElt<Tg, Tn>(this IMap<Tg, Tn> imap, ConcreteGroup<Tn> N)
        where Tg : struct, IElt<Tg>
        where Tn : struct, IElt<Tn>
    {
        return new(imap.Domain, N, new(imap.Domain.ToDictionary(e => e, e => imap[e])));
    }

    private static Tn[] ArrImages<Tn, Tg>(MapElt<Ep<Tg>, Tn> mapElt, ConcreteGroup<Tg> G)
        where Tg : struct, IElt<Tg>
        where Tn : struct, IElt<Tn>
    {
        return mapElt.map.OrderKeys(G).Select(m => m.Value).ToArray();
    }

    public static IOrderedEnumerable<MapElt<Ep<Tg>, Tn>> OrderMaps<Tn, Tg>(this IEnumerable<MapElt<Ep<Tg>, Tn>> maps,
        ConcreteGroup<Tg> G)
        where Tg : struct, IElt<Tg>
        where Tn : struct, IElt<Tn>
    {
        return maps.OrderBy(m => m,
            Comparer<MapElt<Ep<Tg>, Tn>>.Create((m0, m1) => ArrImages(m0, G).SequenceCompareTo(ArrImages(m1, G))));
    }

    public static Dictionary<MapElt<Tg, Automorphism<Tn>>, TwoCocyclesDFS<Tn, Tg>> TwoCocycles<Tn, Tg>(ConcreteGroup<Tn> N,
        ConcreteGroup<Tg> G, bool trivialActionOnly = true)
        where Tg : struct, IElt<Tg>
        where Tn : struct, IElt<Tn>
    {
        var autN = Group.AutomorphismGroup(N);
        var allOps = Group.AllHomomorphisms(G, autN);
        var trivL = new Homomorphism<Tg, Automorphism<Tn>>(G, G.ToDictionary(e => e, _ => autN.Neutral()));
        if (trivialActionOnly)
            allOps = new() { trivL };
        return allOps.Select(L => new MapElt<Tg, Automorphism<Tn>>(G, autN, new(L.HomMap)))
            .Select((L, lbl) => new TwoCocyclesDFS<Tn, Tg>(N, G, L, $"Lbl{lbl + 1}/{allOps.Count}"))
            .ToDictionary(e => e.L, e => e);
    }

    public static TwoCocyclesDFS<Tn, Tg> TwoCocycles<Tn, Tg>(ConcreteGroup<Tn> N, ConcreteGroup<Tg> G, MapElt<Tg, Automorphism<Tn>> L,
        string lbl = "test")
        where Tg : struct, IElt<Tg>
        where Tn : struct, IElt<Tn>
    {
        return new TwoCocyclesDFS<Tn, Tg>(N, G, L, lbl);
    }

    
}
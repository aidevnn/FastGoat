using System.Collections;
using FastGoat.Structures.CartesianProduct;

namespace FastGoat.Structures.GenericGroup;

public readonly struct ExtensionGroup<Tn, Tg> : IGroup<Ep2<Tn, Tg>> where Tg : struct, IElt<Tg> where Tn : struct, IElt<Tn>
{
    public MapGroups<Ep2<Tg, Tg>, Tn> Map { get; }
    public ConcreteGroup<Tg> G { get; }
    public ConcreteGroup<Tn> N { get; }
    private HashSet<Ep2<Tn, Tg>> Elements { get; }

    public ExtensionGroup(ConcreteGroup<Tn> n, MapGroups<Ep2<Tg, Tg>, Tn> map, ConcreteGroup<Tg> g)
    {
        G = g;
        N = n;
        Map = map;
        var nName = N.Name.Contains('x') ? $"({N})" : $"{N}";
        var gName = G.Name.Contains('x') ? $"({G})" : $"{G}";
        Name = $"{nName} . {gName}";
        Elements = Product.Generate(N, G).ToHashSet();
        Hash = (Name, Elements.Count).GetHashCode();
    }

    public IEnumerator<Ep2<Tn, Tg>> GetEnumerator() => Elements.GetEnumerator();

    IEnumerator IEnumerable.GetEnumerator()
    {
        return GetEnumerator();
    }

    public bool Equals(IGroup<Ep2<Tn, Tg>>? other) => other?.Hash == Hash;

    public int Hash { get; }
    public string Name { get; }

    public Ep2<Tn, Tg> this[params ValueType[] us] => new(N[us[0]], G[us[1]]);

    public IEnumerable<Ep2<Tn, Tg>> GetElements() => Elements;

    public IEnumerable<Ep2<Tn, Tg>> GetGenerators()
    {
        foreach (var g0 in G)
            yield return new(N.Neutral(), g0);

        foreach (var n0 in N)
            yield return new(n0, G.Neutral());
    }

    public Ep2<Tn, Tg> Neutral() => new(N.Neutral(), G.Neutral());

    public Ep2<Tn, Tg> Invert(Ep2<Tn, Tg> e)
    {
        var (n, t) = (e.E1, e.E2);
        var ni = N.Invert(n);
        var s = G.Invert(t);
        var wi = N.Invert(Map[new(s, t)]);
        var m = N.Op(wi, ni);
        return new(m, s);
    }

    public Ep2<Tn, Tg> Op(Ep2<Tn, Tg> e1, Ep2<Tn, Tg> e2)
    {
        // “twisted” multiplication ·ω on N × G by (m, s) ·ω (n, t) := (mnω(s, t), st);
        var (m, n, s, t) = (e1.E1, e2.E1, e1.E2, e2.E2);
        var st = G.Op(s, t);
        var mn = N.Op(m, n);
        return new(N.Op(mn, Map[new(s, t)]), st);
    }

    public override int GetHashCode() => Hash;
    public override string ToString() => Name;
}
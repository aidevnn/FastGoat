using System.Collections;
using FastGoat.Structures.CartesianProduct;

namespace FastGoat.Structures.GenericGroup;

public struct ExtensionGroupBase<Tn, Tg> : IGroup<Ep2<Tn, Tg>> where Tg : struct, IElt<Tg> where Tn : struct, IElt<Tn>
{
    public MapElt<Ep2<Tg, Tg>, Tn> Map { get; }
    public MapElt<Tg, Automorphism<Tn>> L { get; }
    public ConcreteGroup<Tg> G { get; }
    public ConcreteGroup<Tn> N { get; }

    public ExtensionGroupBase(ConcreteGroup<Tn> n, MapElt<Tg, Automorphism<Tn>> l, MapElt<Ep2<Tg, Tg>, Tn> map, ConcreteGroup<Tg> g)
    {
        G = g;
        N = n;
        Map = map;
        L = l.Clone();
        var nName = N.Name.Contains('x') ? $"({N})" : $"{N}";
        var gName = G.Name.Contains('x') ? $"({G})" : $"{G}";
        Name = $"{nName} . {gName}";

        Elements = Product.Generate(N, G).ToHashSet();

        PseudoGenerators = new List<Ep2<Tn, Tg>>();
        foreach (var e in G.GetGenerators())
            PseudoGenerators.Add(Product.Elt(N.Neutral(), e));

        foreach (var e in N.GetGenerators())
            PseudoGenerators.Add(Product.Elt(e, G.Neutral()));

        Hash = (Name, Elements.Count).GetHashCode();
        IsGroup = Group.IsGroup(this, Elements);
        // IsGroup = true;
    }
    
    public bool IsGroup { get; }

    public List<Ep2<Tn, Tg>> PseudoGenerators { get; }
    public HashSet<Ep2<Tn, Tg>> Elements { get; }

    public IEnumerator<Ep2<Tn, Tg>> GetEnumerator() => GetElements().GetEnumerator();

    IEnumerator IEnumerable.GetEnumerator()
    {
        return GetEnumerator();
    }

    public bool Equals(IGroup<Ep2<Tn, Tg>>? other) => other?.Hash == Hash;

    public int Hash { get; }
    public string Name { get; }

    public Ep2<Tn, Tg> this[params ValueType[] us]
    {
        get
        {
            var (us0, us1) = (us[0], us[1]);
            return new(N[us0], G[us1]);
        }
    }

    public IEnumerable<Ep2<Tn, Tg>> GetElements() => Elements;

    public IEnumerable<Ep2<Tn, Tg>> GetGenerators() => PseudoGenerators;

    public Ep2<Tn, Tg> Neutral() => new(N.Neutral(), G.Neutral());

    public Ep2<Tn, Tg> Invert(Ep2<Tn, Tg> e)
    {
        var (n, g) = (e.E1, e.E2);
        var ni = N.Invert(n);
        var gi = G.Invert(g);
        var nop = N.Op(Map[new(gi, g)], ni);
        var m = L[g][nop];
        return new(m, gi);
    }

    public Ep2<Tn, Tg> Op(Ep2<Tn, Tg> e1, Ep2<Tn, Tg> e2)
    {
        var (n1, n2, g1, g2) = (e1.E1, e2.E1, e1.E2, e2.E2);
        var nop = N.Op(L[G.Invert(g2)][n1], n2);
        var fi = N.Invert(Map[new(g1, g2)]);
        var g1g2 = G.Op(g1, g2);
        return new(N.Op(fi, nop), g1g2);
    }

    public override int GetHashCode() => Hash;
    public override string ToString() => Name;
}
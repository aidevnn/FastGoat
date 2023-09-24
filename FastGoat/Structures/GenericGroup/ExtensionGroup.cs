using System.Collections;
using FastGoat.Structures.CartesianProduct;

namespace FastGoat.Structures.GenericGroup;

public class ExtensionGroup<Tn, Tg> : ConcreteGroup<Ep2<Tn, Tg>> where Tg : struct, IElt<Tg> where Tn : struct, IElt<Tn>
{
    public MapGroups<Ep2<Tg, Tg>, Tn> Map { get; }
    public MapGroups<Tg, Automorphism<Tn>> L { get; }
    public ConcreteGroup<Tg> G { get; }
    public ConcreteGroup<Tn> N { get; }

    public ExtensionGroup(ConcreteGroup<Tn> n, MapGroups<Tg, Automorphism<Tn>> l, MapGroups<Ep2<Tg, Tg>, Tn> map, ConcreteGroup<Tg> g)
    :base("Ext",  Product.Group(n, g), true)
    {
        G = g;
        N = n;
        Map = map;
        L = l.Clone();
        var nName = N.Name.Contains('x') ? $"({N})" : $"{N}";
        var gName = G.Name.Contains('x') ? $"({G})" : $"{G}";
        Name = $"{nName} . {gName}";

        List<Ep2<Tn, Tg>> generators = new List<Ep2<Tn, Tg>>();
        foreach (var e in G.GetGenerators())
            generators.Add(Product.Elt(N.Neutral(), e));

        foreach (var e in N.GetGenerators())
            generators.Add(Product.Elt(e, G.Neutral()));
        
        var (tmpElements, uniqueGenerators) = Group.UniqueGenerators(this, generators.ToArray());
        PseudoGenerators = new(uniqueGenerators);
        Elements = tmpElements;
        ElementsOrders = Group.ElementsOrders(this, Elements);
        GroupType = Group.IsCommutative(this, PseudoGenerators)
            ? GroupType.AbelianGroup
            : GroupType.NonAbelianGroup;
        Hash = (Name, Elements.Count).GetHashCode();
    }

    public override Ep2<Tn, Tg> Neutral() => new(N.Neutral(), G.Neutral());

    public override Ep2<Tn, Tg> Invert(Ep2<Tn, Tg> e)
    {
        // var (n, t) = (e.E1, e.E2);
        // var ni = N.Invert(n);
        // var s = G.Invert(t);
        // var wi = N.Invert(Map[new(s, t)]);
        // var m = N.Op(wi, ni);
        // return new(m, s);
        var (n, g) = (e.E1, e.E2);
        var ni = N.Invert(n);
        var gi = G.Invert(g);
        var nop = N.Op(Map[new(gi, g)], ni);
        var m = L[g][nop];
        return new(m, gi);
    }

    public override Ep2<Tn, Tg> Op(Ep2<Tn, Tg> e1, Ep2<Tn, Tg> e2)
    {
        // “twisted” multiplication ·ω on N × G by (m, s) ·ω (n, t) := (mnω(s, t), st);
        // var (m, n, s, t) = (e1.E1, e2.E1, e1.E2, e2.E2);
        // var st = G.Op(s, t);
        // var mn = N.Op(m, n);
        // return new(N.Op(mn, Map[new(s, t)]), st);
        var (n1, n2, g1, g2) = (e1.E1, e2.E1, e1.E2, e2.E2);
        var nop = N.Op(L[g2].Invert()[n1], n2);
        var fi = N.Invert(Map[new(g1, g2)]);
        var g1g2 = G.Op(g1, g2);
        return new(N.Op(fi, nop), g1g2);
    }

    public override int GetHashCode() => Hash;
    public override string ToString() => Name;
}
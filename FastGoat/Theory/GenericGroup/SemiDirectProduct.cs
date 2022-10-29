using FastGoat.Theory.CartesianProduct;

namespace FastGoat.Theory.GenericGroup;

public class SemiDirectProduct<T1, T2> : ConcreteGroup<Ep2<T1, T2>>
    where T1 : struct, IElt<T1> where T2 : struct, IElt<T2>
{
    public Homomorphism<T2, Automorphism<T1>> Theta { get; }
    public ConcreteGroup<T1> N { get; }
    public ConcreteGroup<T2> G { get; }
    public ConcreteGroup<Ep2<T1, T2>> Ncan { get; }
    public ConcreteGroup<Ep2<T1, T2>> Gcan { get; }

    public SemiDirectProduct(string name, ConcreteGroup<T1> n, Homomorphism<T2, Automorphism<T1>> theta,
        ConcreteGroup<T2> g) : base(name, Product.Group(n, g), true)
    {
        G = g;
        N = n;
        if (string.IsNullOrEmpty(name))
        {
            var nName = n.Name.Contains(' ') ? $"({n.Name})" : n.Name;
            var gName = g.Name.Contains(' ') ? $"({g.Name})" : g.Name;
            Name = $"{nName} x: {gName}";
        }
        else
        {
            Name = name;
        }

        Theta = theta;
        List<Ep2<T1, T2>> generators = new List<Ep2<T1, T2>>();
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

        Ncan = Group.Generate(N.Name, this, n.Select(e => Product.Elt(e, G.Neutral())).ToArray());
        Gcan = Group.Generate(G.Name, this, g.Select(e => Product.Elt(N.Neutral(), e)).ToArray());
    }

    public Ep2<T1, T2> Act(T1 en, T2 eg)
    {
        if (!N.Contains(en) || !G.Contains(eg))
            throw new GroupException(GroupExceptionType.BaseGroup);

        return Product.Elt(Theta[eg][en], eg);
    }

    public override Ep2<T1, T2> Invert(Ep2<T1, T2> e)
    {
        var gi = G.Invert(e.E2);
        var xi = N.Invert(e.E1);
        return Product.Elt(Theta[gi][xi], gi);
    }

    public override Ep2<T1, T2> Op(Ep2<T1, T2> e1, Ep2<T1, T2> e2)
    {
        var n = N.Op(e1.E1, Theta[e1.E2][e2.E1]);
        var g = G.Op(e1.E2, e2.E2);
        return Product.Elt(n, g);
    }
}
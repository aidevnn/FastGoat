using FastGoat.Structures.CartesianProduct;

namespace FastGoat.Structures.GenericGroup;

public class SemiDirectProduct<T1, T2> : ConcreteGroup<Ep2<T1, T2>>
    where T1 : struct, IElt<T1> where T2 : struct, IElt<T2>
{
    public Homomorphism<T2, Automorphism<T1>> Theta { get; }
    public ConcreteGroup<T1> N { get; }
    public ConcreteGroup<T2> G { get; }
    public ConcreteGroup<Ep2<T1, T2>> Ncan { get; }
    public ConcreteGroup<Ep2<T1, T2>> Gcan { get; }
    public Dictionary<Ep2<T1, T2>, Ep2<T1, T2>> InvertTable { get; }
    public Dictionary<(Ep2<T1, T2>, Ep2<T1, T2>), Ep2<T1, T2>> OpTable { get; }

    public SemiDirectProduct(string name, ConcreteGroup<T1> n, Homomorphism<T2, Automorphism<T1>> theta,
        ConcreteGroup<T2> g) : base(name, Product.Group(n, g), true)
    {
        G = g;
        N = n;
        Name = !string.IsNullOrEmpty(name) ? name : $"{n.NameParenthesis()} x: {g.NameParenthesis()}";

        Theta = theta;
        List<Ep2<T1, T2>> generators = new List<Ep2<T1, T2>>();
        foreach (var e in G.GetGenerators())
            generators.Add(Product.Elt(N.Neutral(), e));

        foreach (var e in N.GetGenerators())
            generators.Add(Product.Elt(e, G.Neutral()));

        Ord = N.Count() * G.Count();
        if (Ord < Group.GetStorageCapacity())
        {
            InvertTable = new(2 * Ord);
            OpTable = new(2 * Ord * Ord);
        }
        else
        {
            InvertTable = new();
            OpTable = new();
        }
        
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

    private int Ord { get; }
    public bool IsFaithFull() => Theta.Kernel().Count() == 1;

    public Ep2<T1, T2> Act(T1 en, T2 eg)
    {
        if (!N.Contains(en) || !G.Contains(eg))
            throw new GroupException(GroupExceptionType.BaseGroup);

        return Product.Elt(Theta[eg][en], eg);
    }

    public override Ep2<T1, T2> Invert(Ep2<T1, T2> e)
    {
        if (Ord >= Group.GetStorageCapacity())
        {
            var gi = G.Invert(e.E2);
            var xi = N.Invert(e.E1);
            return Product.Elt(Theta[gi][xi], gi);
        }
        else
        {
            if (InvertTable.TryGetValue(e, out var r))
                return r;
        
            var gi = G.Invert(e.E2);
            var xi = N.Invert(e.E1);
            var ei = InvertTable[e] = Product.Elt(Theta[gi][xi], gi);
            return ei;
        }
    }

    public override Ep2<T1, T2> Op(Ep2<T1, T2> e1, Ep2<T1, T2> e2)
    {
        if (Ord >= Group.GetStorageCapacity())
        {
            var n = N.Op(e1.E1, Theta[e1.E2][e2.E1]);
            var g = G.Op(e1.E2, e2.E2);
            return Product.Elt(n, g);
        }
        else
        {
            var e12 = (e1, e2);
            if (OpTable.TryGetValue(e12, out var r))
                return r;
        
            var n = N.Op(e1.E1, Theta[e1.E2][e2.E1]);
            var g = G.Op(e1.E2, e2.E2);
            var e3 = OpTable[e12] = Product.Elt(n, g);
            return e3;
        }
    }
}

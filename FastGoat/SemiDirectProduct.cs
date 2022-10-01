using System.Collections.ObjectModel;
using FastGoat.Gp;

namespace FastGoat;

public class SemiDirectProduct<T1, T2> : ConcreteGroup<Ep2<T1, T2>>
    where T1 : struct, IElt<T1> where T2 : struct, IElt<T2>
{
    public ReadOnlyDictionary<T2, Func<T1, T1>> Actions { get; }
    public ReadOnlyDictionary<T2, string> ActionsStr { get; }
    public ConcreteGroup<T1> N { get; }
    public ConcreteGroup<T2> G { get; }

    public SemiDirectProduct(string name, ConcreteGroup<T1> n, ConcreteGroup<T2> g) : base(name,
        Product.Elt(n.Neutral(), g.Neutral()).BaseGroup, true)
    {
        if (g.LongestCycles.Count != 1 || n.LongestCycles.Count != 1)
            throw new GroupException(GroupExceptionType.OnlyCyclicGroups);

        G = g;
        N = n;
        var cg = G.Count();
        var cn = N.Count();
        var pow = IntExt.Solve_x_pow_n_equal_one_mod_m(cn, cg);
        if (pow == 0)
            throw new GroupException(GroupExceptionType.SemiDirectProductDontExist);

        var actions = new Dictionary<T2, Func<T1, T1>>();
        var actionsStr = new Dictionary<T2, string>();
        actions[g.Neutral()] = x => x;
        actionsStr[g.Neutral()] = "x->x";
        foreach (var (_, set) in G.LongestCycles)
        {
            foreach (var (e, p) in set.SkipLast(1))
            {
                var p0 = IntExt.PowMod(pow, p, cn);
                actions[e] = x => N.Times(x, p0);
                actionsStr[e] = p0 == 1 ? "x->x" : $"x->x^{p0}";
            }
        }

        Actions = new ReadOnlyDictionary<T2, Func<T1, T1>>(actions);
        ActionsStr = new ReadOnlyDictionary<T2, string>(actionsStr);

        List<Ep2<T1, T2>> generators = new List<Ep2<T1, T2>>();
        foreach (var (e, _) in G.LongestCycles)
            generators.Add(Product.Elt(N.Neutral(), e));

        foreach (var (e, _) in N.LongestCycles)
            generators.Add(Product.Elt(e, G.Neutral()));

        Elements = Group.GenerateElements(this, generators.ToArray()).ToHashSet();
        LongestCycles = Group.LongestCycles(this, Elements);
        ElementsOrders = Group.ElementsOrders(LongestCycles);
        GroupType = Group.IsCommutative(this, LongestCycles.Keys)
            ? GroupType.AbelianGroup
            : GroupType.NonAbelianGroup;
    }

    public SemiDirectProduct(ConcreteGroup<T1> n, ConcreteGroup<T2> g) : this($"{n} x: {g}", n, g)
    {
    }

    public override Ep2<T1, T2> Invert(Ep2<T1, T2> e)
    {
        var gi = G.Invert(e.E2);
        var xi = N.Invert(e.E1);
        return Product.Elt(Actions[gi](xi), gi);
    }

    public override Ep2<T1, T2> Op(Ep2<T1, T2> e1, Ep2<T1, T2> e2)
    {
        var n = N.Op(e1.E1, Actions[e1.E2](e2.E1));
        var g = G.Op(e1.E2, e2.E2);
        return Product.Elt(n, g);
    }
}
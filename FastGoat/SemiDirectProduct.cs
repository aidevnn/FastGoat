namespace FastGoat;

public class SemiDirectProduct<T1, T2> : ConcreteGroup<Ep<T1, T2>> where T1 : struct, IElt<T1> where T2 : struct, IElt<T2>
{
    public SemiDirectProduct(WorkGroup<T1> n, WorkGroup<T2> g, Func<T2, T1, T1> action) : base(new Gp<T1, T2>(n.BaseGroup, g.BaseGroup))
    {
        ControlGroup = new ConcreteGroup<Ep<T1, T2>>(BaseGroup);
        N = n;
        G = g;
        this.action = action;
        var nGens = n.Select(n0 => new Ep<T1, T2>(n0, g.Neutral()));
        var gGens = g.Select(g0 => new Ep<T1, T2>(n.Neutral(), g0));
        var tmpElements = Generate(nGens, gGens);
        (groupType, elementOrder, monogenics) = ComputeDetails(tmpElements);
        elements = new(tmpElements);
    }
    WorkGroup<T1> N { get; }
    WorkGroup<T2> G { get; }
    Func<T2, T1, T1> action { get; }
    public override Ep<T1, T2> Neutral() => (N.Neutral(), G.Neutral());
    public override Ep<T1, T2> Invert(Ep<T1, T2> a)
    {
        var g = a.e2;
        var gi = G.Invert(g);
        var xi = N.Invert(a.e1);
        return (action(gi, xi), gi);
    }
    public override Ep<T1, T2> Op(Ep<T1, T2> a, Ep<T1, T2> b)
    {
        var g = a.e2;
        var x = a.e1;
        var h = b.e2;
        var y = b.e1;
        return (N.Op(x, action(g, y)), G.Op(g, h));
    }
}

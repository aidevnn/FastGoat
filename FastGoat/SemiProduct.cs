namespace FastGoat;

public class SemiProduct<T1, T2> : WorkGroup<Ep<T1, T2>> where T1 : struct, IElt<T1> where T2 : struct, IElt<T2>
{
    public SemiProduct(WorkGroup<T1> n, WorkGroup<T2> g) : base(new Gp<T1, T2>(n.BaseGroup, g.BaseGroup))
    {
        if (n.GetMonogenics().Count() != 1 || g.GetMonogenics().Count() != 1)
            throw new Exception("Only cyclic groups are supported");

        ControlGroup = new ConcreteGroup<Ep<T1, T2>>(BaseGroup);
        N = n;
        G = g;

        action = new();
        var arrG = G.Ascending().ToArray();
        var cn = N.Count();
        var cg = G.Count();

        var seq = Enumerable.Range(2, cn - 2);
        var criteria = seq.Where(i => IntExt.GCD(i, cn) == 1 && IntExt.PowMod(i, cg, cn) == 1);
        if (criteria.Count() == 0)
            throw new Exception("Semiproduct doesnt exist");

        int pow = criteria.First();
        action[G.Neutral()] = x => x;
        for (int k = 1; k < arrG.Length; ++k)
        {
            int k0 = IntExt.PowMod(pow, k, cn);
            action[arrG[k]] = x => N.Pow(x, k0);
        }

        var nGens = n.Select(n0 => new Ep<T1, T2>(n0, g.Neutral()));
        var gGens = g.Select(g0 => new Ep<T1, T2>(n.Neutral(), g0));
        Ncan = Group.Generate(nGens.ToArray());
        Gcan = Group.Generate(gGens.ToArray());
        var tmpElements = Generate(nGens, gGens);
        (groupType, elementOrder, monogenics) = ComputeDetails(tmpElements);
        elements = new(tmpElements);
    }

    public WorkGroup<Ep<T1, T2>> Ncan { get; }
    public WorkGroup<Ep<T1, T2>> Gcan { get; }
    WorkGroup<T1> N { get; }
    WorkGroup<T2> G { get; }
    Dictionary<T2, Func<T1, T1>> action { get; }
    public override Ep<T1, T2> Neutral() => (N.Neutral(), G.Neutral());
    public override Ep<T1, T2> Invert(Ep<T1, T2> a)
    {
        var g = a.e2;
        var gi = G.Invert(g);
        var xi = N.Invert(a.e1);
        return (action[gi](xi), gi);
    }
    public override Ep<T1, T2> Op(Ep<T1, T2> a, Ep<T1, T2> b)
    {
        var g = a.e2;
        var x = a.e1;
        var h = b.e2;
        var y = b.e1;
        return (N.Op(x, action[g](y)), G.Op(g, h));
    }

    public void GNGi_it()
    {
        Console.WriteLine();
        Console.WriteLine("GNGi_it !!!");
        foreach (var g0 in G)
        {
            Ep<T1, T2> g = (N.Neutral(), g0);
            Ep<T1, T2> gi = (N.Neutral(), G.Invert(g0));
            foreach (var n0 in N)
            {
                Ep<T1, T2> n = (n0, G.Neutral());
                var gngi = this.Op(this.Op(g, n), gi);
                var gngi_it = action[g0](n0);
                Console.WriteLine("g={0} n={1}; y(g)(n) = {2} ? {3} = gngi", g0, n0, gngi_it, gngi);
            }
        }

        Console.WriteLine();
    }
}

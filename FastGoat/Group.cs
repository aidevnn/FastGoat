namespace FastGoat;

public static class Group
{
    public static Gp<U1, U2> CartesianProduct<U1, U2>(IGroup<U1> g1, IGroup<U2> g2)
        where U1 : struct, IElt<U1>
        where U2 : struct, IElt<U2>
    {
        return new Gp<U1, U2>(g1, g2);
    }
    public static Gp<U1, U2, U3> CartesianProduct<U1, U2, U3>(IGroup<U1> g1, IGroup<U2> g2, IGroup<U3> g3)
        where U1 : struct, IElt<U1>
        where U2 : struct, IElt<U2>
        where U3 : struct, IElt<U3>
    {
        return new Gp<U1, U2, U3>(g1, g2, g3);
    }
    public static Gp<U1, U2, U3, U4> CartesianProduct<U1, U2, U3, U4>(IGroup<U1> g1, IGroup<U2> g2, IGroup<U3> g3, IGroup<U4> g4)
        where U1 : struct, IElt<U1>
        where U2 : struct, IElt<U2>
        where U3 : struct, IElt<U3>
        where U4 : struct, IElt<U4>
    {
        return new Gp<U1, U2, U3, U4>(g1, g2, g3, g4);
    }
    public static WorkGroup<T> Generate<T>(params T[] te) where T : struct, IElt<T>
    {
        if (te.Length == 0)
            throw new Exception("Elements must be specified");

        return new WorkGroup<T>(te);
    }
    public static WorkGroup<T> Generate<T>(this IGroup<T> gr, params T[] te) where T : struct, IElt<T>
    {
        if (te.Length == 0)
            return new WorkGroup<T>(gr.Neutral());

        if (!gr.Equals(te.First().Group))
            throw new BaseGroupException();

        return new WorkGroup<T>(te);
    }
    public static WorkGroup<Ep<U1, U2>> Generate<U1, U2>(this Gp<U1, U2> gr, params Ep<U1, U2>[] te)
        where U1 : struct, IElt<U1>
        where U2 : struct, IElt<U2>
    {
        if (te.Length == 0)
            return new WorkGroup<Ep<U1, U2>>(gr.Neutral());

        if (!gr.Equals(te.First().Group))
            throw new BaseGroupException();

        return new WorkGroup<Ep<U1, U2>>(te);
    }
    public static WorkGroup<Ep<U1, U2, U3>> Generate<U1, U2, U3>(this Gp<U1, U2, U3> gr, params Ep<U1, U2, U3>[] te)
        where U1 : struct, IElt<U1>
        where U2 : struct, IElt<U2>
        where U3 : struct, IElt<U3>
    {
        if (te.Length == 0)
            return new WorkGroup<Ep<U1, U2, U3>>(gr.Neutral());

        if (!gr.Equals(te.First().Group))
            throw new BaseGroupException();

        return new WorkGroup<Ep<U1, U2, U3>>(te);
    }
    public static WorkGroup<Ep<U1, U2, U3, U4>> Generate<U1, U2, U3, U4>(this Gp<U1, U2, U3, U4> gr, params Ep<U1, U2, U3, U4>[] te)
        where U1 : struct, IElt<U1>
        where U2 : struct, IElt<U2>
        where U3 : struct, IElt<U3>
        where U4 : struct, IElt<U4>
    {
        if (te.Length == 0)
            return new WorkGroup<Ep<U1, U2, U3, U4>>(gr.Neutral());

        if (!gr.Equals(te.First().Group))
            throw new BaseGroupException();

        return new WorkGroup<Ep<U1, U2, U3, U4>>(te);
    }
    public static QuotientGroup<T> Over<T>(this WorkGroup<T> g, WorkGroup<T> h) where T : struct, IElt<T>
    {
        try
        {
            return new QuotientGroup<T>(g, h);
        }
        catch (Exception e)
        {
            Console.WriteLine("{0} : {1}", e.GetType(), e.Message);
        }

        return new QuotientGroup<T>(g, g);
    }
    static (T factor, int order, WorkGroup<T> subGroup) InvariantFactor<T>(this WorkGroup<T> group) where T : struct, IElt<T>
    {
        if (group.groupType == GroupType.NotAbelianGroup)
            throw new Exception("Not Abelian group");

        var eMax = group.OrderByDescending(group.GetOrderOf).ThenAscending().First();
        var H = group.GenerateSubgroup(eMax);
        return (eMax, group.GetOrderOf(eMax), group.Over(H));
    }
    public static void InvariantFactors<T>(this WorkGroup<T> group) where T : struct, IElt<T>
    {
        var g0 = new WorkGroup<T>(group);
        var e0 = group.Last();
        var o = 1;
        List<string> decom = new();
        var name = "G";
        while (!e0.Equals(g0.Neutral()))
        {
            var name0 = name;
            var count = g0.Count();
            (e0, o, g0) = g0.InvariantFactor();
            var co = $"C{o}";
            decom.Add(co);
            name = $"{name}/{co}";

            if (!e0.Equals(g0.Neutral()))
                Console.WriteLine("{0,18} = {1,-5} max order element is {2,-5} = {3}", name0, count, co, e0);

        }

        Console.WriteLine();
        Console.WriteLine("G[{0}] ~ {1}", group.Count(), decom.SkipLast(1).Glue(" x "));
        Console.WriteLine($"in {group.BaseGroup}");
        Console.WriteLine();
    }

    public static SemiDirectProduct<U1, U2> SemiDirectProd<U1, U2>(WorkGroup<U1> n, WorkGroup<U2> g, Func<U2, U1, U1> action)
        where U1 : struct, IElt<U1>
        where U2 : struct, IElt<U2>
    {
        return new SemiDirectProduct<U1, U2>(n, g, action);
    }
}

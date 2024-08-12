using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.GenericGroup;
using FastGoat.UserGroup;
using FastGoat.UserGroup.GModuleN;

namespace FastGoat.Examples;

public static class Cohomology
{
    static Cohomology()
    {
        Logger.Level = LogLevel.Level1;
    }

    static void SolveRCohomologies<Tn, Tg>(ConcreteGroup<Tn> N, ConcreteGroup<Tg> G, int r = 2)
        where Tg : struct, IElt<Tg>
        where Tn : struct, IElt<Tn>
    {
        var CN = Group.Zentrum(N);
        var (autN, innN, outN) = Group.OuterAutomorphismGroup(N);
        var ops = Group.AllHomomorphisms(G, outN);
        foreach (var (op, i) in ops.Select((l, i) => (l, i + 1)))
        {
            var L = op.ToMapElt(outN);
            var lbl = $"H{r}(G,N) Lbl{i}/{ops.Count}";
            var (hr, (_, _, br), (_, _, zr)) = ZNSolver.ReduceCohomologies(CN, G, L, r, lbl);
        }
    }

    public static void H4_C2_C2()
    {
        var (N, G) = (FG.Abelian(2), FG.Abelian(2));
        SolveRCohomologies(N, G, r: 1);
        SolveRCohomologies(N, G, r: 2);
        SolveRCohomologies(N, G, r: 3);
        SolveRCohomologies(N, G, r: 4);
    }

    public static void H4_C2_C4()
    {
        var (N, G) = (FG.Abelian(4), FG.Abelian(2));
        SolveRCohomologies(N, G, r: 1);
        SolveRCohomologies(N, G, r: 2);
        SolveRCohomologies(N, G, r: 3);
        SolveRCohomologies(N, G, r: 4);
    }

    public static void H4_C2_C2C2()
    {
        var (N, G) = (FG.Abelian(2, 2), FG.Abelian(2));
        SolveRCohomologies(N, G, r: 1);
        SolveRCohomologies(N, G, r: 2);
        SolveRCohomologies(N, G, r: 3);
        SolveRCohomologies(N, G, r: 4);
    }

    public static void H3_C2C2_C2()
    {
        var (N, G) = (FG.Abelian(2), FG.Abelian(2, 2));
        SolveRCohomologies(N, G, r: 1);
        SolveRCohomologies(N, G, r: 2);
        SolveRCohomologies(N, G, r: 3);
    }

    public static void H3_C2C2_C2C2()
    {
        var (N, G) = (FG.Abelian(2, 2), FG.Abelian(2, 2));
        SolveRCohomologies(N, G, r: 1);
        SolveRCohomologies(N, G, r: 2);
        SolveRCohomologies(N, G, r: 3);
    }

    public static void H3_C2C2_C4()
    {
        var (N, G) = (FG.Abelian(4), FG.Abelian(2, 2));
        SolveRCohomologies(N, G, r: 1);
        SolveRCohomologies(N, G, r: 2);
        SolveRCohomologies(N, G, r: 3);
    }

    public static void H4_C3_C3()
    {
        var (N, G) = (FG.Abelian(3), FG.Abelian(3));
        SolveRCohomologies(N, G, r: 1);
        SolveRCohomologies(N, G, r: 2);
        SolveRCohomologies(N, G, r: 3);
        SolveRCohomologies(N, G, r: 4);
    }

    public static void H4_C3_C9()
    {
        var (N, G) = (FG.Abelian(9), FG.Abelian(3));
        SolveRCohomologies(N, G, r: 1);
        SolveRCohomologies(N, G, r: 2);
        SolveRCohomologies(N, G, r: 3);
        SolveRCohomologies(N, G, r: 4);
    }

    public static void H4_C3_C3C3()
    {
        var (N, G) = (FG.Abelian(9), FG.Abelian(3));
        SolveRCohomologies(N, G, r: 1);
        SolveRCohomologies(N, G, r: 2);
        SolveRCohomologies(N, G, r: 3);
        SolveRCohomologies(N, G, r: 4);
    }

}

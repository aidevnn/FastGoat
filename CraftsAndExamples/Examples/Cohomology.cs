using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.UserGroup;
using FastGoat.UserGroup.GModuleN;

namespace CraftsAndExamples.Examples;

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

    static void DiagAction<Tg, Tn>(ConcreteGroup<Tn> N, ConcreteGroup<Tg> G, int r)
        where Tn : struct, IElt<Tn>
        where Tg : struct, IElt<Tg>
    {
        var nab = N.AbelianDirectSum();
        var subgs = nab.ElementarySubgroups();

        var autG = Group.AutomorphismGroup(G);
        var autN = Group.AutomorphismGroup(N);
        Console.WriteLine($"{autG.ShortName}    {autN.ShortName}");
        var ops = Group.AllHomomorphisms(G, autN).ToHashSet(FG.EqOpByAut(G, autG, autN));
        var diagOps = ops.Where(op => G.All(g => subgs.All(sg => sg.SetEquals(sg.Select(e => op[g][e]))))).ToArray();
        foreach (var (L, k) in diagOps.Select((op, k) => (op.ToMapElt(autN), k + 1)))
        {
            var lbl = $"Lbl{k}/{diagOps.Length}/{ops.Count}";
            var decomp = subgs.Select(sg => ZNSolver.ReduceCohomologies(sg, G, L, r, lbl)).ToArray();
            var sol = ZNSolver.ReduceCohomologies(N, G, L, r, lbl);

            var nbH2 = sol.solsCohs.Length;
            var prod = decomp.Select(s => s.solsCohs.Length).Aggregate((a0, a1) => a0 * a1);
            if (nbH2 == prod)
            {
                Console.WriteLine("PASS");
                if (sol.solsCocs.total > 1000)
                {
                    Console.WriteLine(" *** TOO BIG *** ");
                    Console.WriteLine();
                    continue;
                }

                var cobsExpected = sol.solsCobs.allMaps.ToGroupMapElt();
                var cocsExpected = sol.solsCocs.allMaps.ToGroupMapElt();
                var cohsExpected = sol.solsCohs.Select(c => c.ToMapElt).ToHashSet();

                var gb = (MapGroupBase<Ep<Tg>, Tn>)cobsExpected.BaseGroup;
                
                var cobsSeq = decomp.Select(e => e.solsCobs.allMaps.ToGroupMapElt()).ToArray();
                var cobsActual = cobsSeq.MultiLoop().Select(l => gb.OpSeq(l)).ToHashSet();
                if (!cobsExpected.SetEquals(cobsActual))
                    throw new("#Coboundaries");

                var cocsSeq = decomp.Select(e => e.solsCocs.allMaps.ToGroupMapElt()).ToArray();
                var cocsActual = cocsSeq.MultiLoop().Select(l => gb.OpSeq(l)).ToHashSet();
                if (!cocsExpected.SetEquals(cocsActual))
                    throw new("#Cocycles");

                var cohsSeq = decomp.Select(e => e.solsCohs.Select(c => c.ToMapElt)).ToArray();
                var cohsActual = cohsSeq.MultiLoop().Select(l => gb.OpSeq(l)).ToHashSet();
                if (!cohsExpected.SetEquals(cohsActual))
                    throw new("#Cohomologies");
            }
            else
                throw new("FAIL");

            Console.WriteLine();
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

    public static void Example_DiagonalAction()
    {
        GlobalStopWatch.Restart();
        for (int r = 1; r <= 3; r++)
        {
            DiagAction(N: FG.Abelian(2, 2), G: FG.Abelian(2), r);
            DiagAction(N: FG.Abelian(2, 4), G: FG.Abelian(2), r);
            DiagAction(N: FG.Abelian(2, 2, 2), G: FG.Abelian(2), r);
            DiagAction(N: FG.Abelian(2, 2, 4), G: FG.Abelian(2), r);
            DiagAction(N: FG.Abelian(4, 4), G: FG.Abelian(2), r);
            
            DiagAction(N: FG.Abelian(3, 3), G: FG.Abelian(3), r);
            DiagAction(N: FG.Abelian(3, 9), G: FG.Abelian(3), r);
            
            DiagAction(N: FG.Abelian(2, 2), G: FG.Abelian(4), r);
            DiagAction(N: FG.Abelian(2, 4), G: FG.Abelian(4), r);
            DiagAction(N: FG.Abelian(2, 2, 2), G: FG.Abelian(4), r);
            
            Console.WriteLine();
        }
        
        GlobalStopWatch.Show("END");
    }
}
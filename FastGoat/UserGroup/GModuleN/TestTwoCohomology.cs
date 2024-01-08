using FastGoat.Commons;
using FastGoat.Examples;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.UserGroup.Integers;

namespace FastGoat.UserGroup.GModuleN;

public static class TestTwoCohomology
{
    private static void SolveTwoCohomOrder<Tn, Tg>(ConcreteGroup<Tn> N, ConcreteGroup<Tg> G)
        where Tg : struct, IElt<Tg>
        where Tn : struct, IElt<Tn>
    {
        var autN = Group.AutomorphismGroup(N);
        var ops = Group.AllHomomorphisms(G, autN);
        foreach (var (op, i) in ops.Select((l, i) => (l, i + 1)))
        {
            var L = op.ToMapElt(autN);
            var lbl = $"Lbl{i}/{ops.Count}";
            ZNSolver.RCohomologyOrder(N, G, L, lbl: lbl);
        }
    }

    private static void SolveTwoCohom<Tn, Tg>(ConcreteGroup<Tn> N, ConcreteGroup<Tg> G)
        where Tg : struct, IElt<Tg>
        where Tn : struct, IElt<Tn>
    {
        var autN = Group.AutomorphismGroup(N);
        var ops = Group.AllHomomorphisms(G, autN);
        foreach (var (op, i) in ops.Select((l, i) => (l, i + 1)))
        {
            var L = op.ToMapElt(autN);
            var lbl = $"Lbl{i}/{ops.Count}";
            ZNSolver.ReduceCohomologies(N, G, L, lbl: lbl);
        }
    }

    private static void SolveTwoCohomExplicit<Tn, Tg>(ConcreteGroup<Tn> N, ConcreteGroup<Tg> G)
        where Tg : struct, IElt<Tg>
        where Tn : struct, IElt<Tn>
    {
        var autN = Group.AutomorphismGroup(N);
        var ops = Group.AllHomomorphisms(G, autN);
        foreach (var (op, i) in ops.Select((l, i) => (l, i + 1)))
        {
            var L = op.ToMapElt(autN);
            var lbl = $"Lbl{i}/{ops.Count}";
            var (cosets, (nb2Cobs, _, maps2Cobs), (nb2Cocs, _, maps2Cocs)) = ZNSolver.ReduceCohomologies(N, G, L, lbl: lbl);
            if (nb2Cocs > 1100)
                continue;

            var cosetsMapElts = cosets.Select(e => e.ToMapElt).ToList();
            var g2Cobs = maps2Cobs.ToGroupMapElt("B2");
            var g2Cocs = maps2Cocs.ToGroupMapElt("Z2");
            var g2Cohs = g2Cocs.Over(g2Cobs, "H2");
            DisplayGroup.Head(g2Cohs);
            // CocyclesDFS.DisplayMapElt("H2 Actual", cosetsMapElts.ToArray());
            // CocyclesDFS.DisplayMapElt("H2 Expected", g2Cohs.Select(e => e.X).ToArray());
            var cosetsExpected = g2Cohs.Select(e => e.ToHashSet()).ToList();
            foreach (var elt in cosetsMapElts)
            {
                var nb = cosetsExpected.RemoveAll(set => set.Contains(elt));
                if (nb != 1)
                    throw new("#############");
            }

            if (cosetsExpected.Count != 0)
                throw new("#############");

            Console.WriteLine($"#### {lbl} COMPLETE");
            Console.WriteLine();
        }
    }

    public static void TestAllAbelianGroup()
    {
        FG.AllAbelianGroupsOfOrder(1).Println();
        FG.AllAbelianGroupsOfOrder(2).Println();
        FG.AllAbelianGroupsOfOrder(3).Println();
        FG.AllAbelianGroupsOfOrder(4).Println();
        FG.AllAbelianGroupsOfOrder(6).Println();
        FG.AllAbelianGroupsOfOrder(8).Println();
        FG.AllAbelianGroupsOfOrder(12).Println();
        FG.AllAbelianGroupsOfOrder(16).Println();
        FG.AllAbelianGroupsOfOrder(18).Println();
    }

    public static void TestOrder()
    {
        var o = new[] { 2, 3, 4, 6, 8, 9, 12, 16 };
        var allProducts = o.Take(3).Grid2D(o).Where(e => e.t1 * e.t2 < 30)
            .SelectMany(e => FG.AllAbelianGroupsOfOrder(e.t2).Grid2D(FG.AllAbelianGroupsOfOrder(e.t1)))
            .ToList();

        allProducts.ForEach(e => SolveTwoCohomOrder(e.t1, e.t2));
    }

    public static void Test0Coset()
    {
        var o = new[] { 2, 3, 4, 6, 8, 9, 12, 16 };
        var allProducts = o.Take(4).Grid2D(o).Where(e => e.t1 * e.t2 <= 40)
            .SelectMany(e => FG.AllAbelianGroupsOfOrder(e.t2).Grid2D(FG.AllAbelianGroupsOfOrder(e.t1)))
            .ToList();

        allProducts.ForEach(e => SolveTwoCohom(e.t1, e.t2));
    }

    public static void Test1Coset()
    {
        var (ord4, ord8) = (FG.AllAbelianGroupsOfOrder(4), FG.AllAbelianGroupsOfOrder(8));
        var ord2_4 = ord4.Prepend(FG.Abelian(2)).ToList();
        var ord2_4_8 = ord2_4.Concat(ord8).ToList();

        {
            var allProducts = ord2_4.Grid2D(ord2_4_8).ToList();
            allProducts.ForEach(e => SolveTwoCohom(e.t1, e.t2));
        }

        {
            var allProducts = ord2_4_8.Grid2D(ord2_4).ToList();
            allProducts.ForEach(e => SolveTwoCohom(e.t1, e.t2));
        }
    }

    public static void Test2Coset()
    {
        var (ord9, ord27) = (FG.AllAbelianGroupsOfOrder(9), FG.AllAbelianGroupsOfOrder(27));
        var ord3_9_27 = ord9.Concat(ord27).Prepend(FG.Abelian(3)).ToList();
        var allProducts = ord3_9_27.Select(e => (e, FG.Abelian(3))).ToList();
        allProducts.ForEach(e => SolveTwoCohom(e.Item1, e.Item2));
    }

    public static void Test3Coset()
    {
        var ord25 = FG.AllAbelianGroupsOfOrder(25);
        var ord5_25 = ord25.Prepend(FG.Abelian(5)).ToList();
        var allProducts = ord5_25.Select(e => (e, FG.Abelian(5))).ToList();
        allProducts.ForEach(e => SolveTwoCohom(e.Item1, e.Item2));
        
        SolveTwoCohom(FG.Abelian(7), FG.Abelian(7));
        SolveTwoCohom(FG.Abelian(14), FG.Abelian(7));
        SolveTwoCohom(FG.Abelian(21), FG.Abelian(7));
    }

    public static void TestExplicit()
    {
        var o = new[] { 2, 3, 4, 6, 8, 9, 12, 16 };
        var allProducts = o.Take(3).Grid2D(o).Where(e => e.t1 * e.t2 < 33)
            .SelectMany(e => FG.AllAbelianGroupsOfOrder(e.t2).Grid2D(FG.AllAbelianGroupsOfOrder(e.t1)))
            .ToList();

        allProducts.ForEach(e => SolveTwoCohomExplicit(e.t1, e.t2));
    }
}

using FastGoat.Commons;
using FastGoat.Examples;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.UserGroup.Integers;

namespace FastGoat.UserGroup.GModuleN;

public static class TestH2Order
{
    private static void SolveTwoCohom<Tn, Tg>(ConcreteGroup<Tn> N, ConcreteGroup<Tg> G, bool details = false)
        where Tg : struct, IElt<Tg>
        where Tn : struct, IElt<Tn>
    {
        var autN = Group.AutomorphismGroup(N);
        var ops = Group.AllHomomorphisms(G, autN);
        foreach (var (op, i) in ops.Select((l, i) => (l, i + 1)))
        {
            var L = op.ToMapElt(autN);
            var lbl = $"Lbl{i}/{ops.Count}";
            ZNSolver.TwoCohomologyOrder(N, G, L, lbl);
        }
    }

    private static List<ConcreteGroup<Ep<ZnInt>>> AllAbelianGroupsOrder(int k)
    {
        var dec = IntExt.PrimesDec(k);
        return dec.Select(e => IntExt.Partitions32[e.Value].Select(l => l.Select(i => e.Key.Pow(i)).ToArray())).MultiLoop()
            .Select(l => FG.Abelian(l.SelectMany(i => i).ToArray())).ToList();
    }

    public static void TestAllAbelianGroup()
    {
        AllAbelianGroupsOrder(2).Println();
        AllAbelianGroupsOrder(3).Println();
        AllAbelianGroupsOrder(4).Println();
        AllAbelianGroupsOrder(6).Println();
        AllAbelianGroupsOrder(8).Println();
        AllAbelianGroupsOrder(12).Println();
        AllAbelianGroupsOrder(16).Println();
        AllAbelianGroupsOrder(18).Println();
    }

    public static void TestCohom()
    {
        var o = new[] { 2, 3, 4, 6, 8, 9, 12, 16 };
        var allProducts = o.Take(3).Grid2D(o).Where(e => e.t1 * e.t2 < 30)
            .SelectMany(e => AllAbelianGroupsOrder(e.t2).Grid2D(AllAbelianGroupsOrder(e.t1)))
            .ToList();

        allProducts.ForEach(e => SolveTwoCohom(e.t1, e.t2));
    }
}
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Integers;

namespace FastGoat.Examples;

public static class AbelianInvariantsFactorsPart2
{
    static int AllAbCount(int k)
    {
        if (k == 1)
            return 1;

        var dec = IntExt.PrimesDec(k);
        return dec.Select(e => IntExt.Partitions32[e.Value].Count).Aggregate(1, (acc, a) => acc * a);
    }

    static int[][] AllAbOfOrder(int n)
    {
        if (n == 1)
            return new[] { new[] { 1 } };

        var dec = IntExt.PrimesDec(n);
        return dec.Select(e => IntExt.Partitions32[e.Value].Select(l => l.Select(i => e.Key.Pow(i)).ToArray()))
            .MultiLoop()
            .Select(e => e.SelectMany(c => c).Order().ToArray())
            .ToArray();
    }

    static int[] ToFactors(int n)
    {
        if (n < 1)
            throw new();

        if (n == 1)
            return new[] { 1 };

        return IntExt.PrimesDec(n).Select(e => e.Key.Pow(e.Value)).Order().ToArray();
    }

    static int[] SeqToFactors(int[] seq)
    {
        return seq.SelectMany(e => ToFactors(e)).Order().ToArray();
    }

    static int[] FactorsToType(int[] seq)
    {
        if (seq.Length == 1 && seq[0] == 1)
            return new[] { 1 };

        var set = seq.ToList();
        while (true)
        {
            var arr = set.Descending().ToArray();
            var (e, l) = set.Descending().Select(e => (e, arr.FirstOrDefault(c => IntExt.Gcd(e, c) == 1, -1)))
                .FirstOrDefault(e => e.Item2 != -1, (-1, -1));
            if (e == -1)
                break;

            set.Remove(e);
            set.Remove(l);
            set.Add(e * l);
        }

        return set.Descending().ToArray();
    }

    static (int[] can, int[] facts) AbType(params int[] seq)
    {
        var facts = SeqToFactors(seq);
        var can = FactorsToType(facts);
        return (can, facts);
    }

    static (int o, int nb)[] CyclicElemsOrders(int n)
    {
        return n.Range()
            .Select(k => k == 0 ? (0, 1) : (k, IntExt.Lcm(n, k) / k))
            .GroupBy(e => e.Item2)
            .Select(e => (e.Key, e.Count()))
            .ToArray();
    }

    static (int o, int nb)[] AbGroupElemsOrders(int[] seq)
    {
        var ords = Array.Empty<(int o, int nb)>();
        foreach (var n in seq)
        {
            var ordsN = CyclicElemsOrders(n);
            if (ords.Length == 0)
                ords = ordsN;
            else
                ords = ords.Grid2D(ordsN)
                    .Select(e => (IntExt.Lcm(e.t1.o, e.t2.o), e.t1.nb * e.t2.nb))
                    .GroupBy(e => e.Item1)
                    .Select(e => (e.Key, e.Sum(a => a.Item2)))
                    .ToArray();
        }

        return ords;
    }

    static (T, int)[] GensToFactors<T>(ConcreteGroup<T> g, T e) where T : struct, IElt<T>
    {
        var elemsOrds = g.ElementsOrders;
        var g0 = Group.Cycle(g, e);
        var facts = ToFactors(elemsOrds[e]);
        return facts.Select(k => g0.First(a => elemsOrds[a.Key] == k)).Select(a => (a.Key, elemsOrds[a.Key])).ToArray();
    }

    static ((T, int)[] can, (T, int)[] facts) AbInvariants<T>(ConcreteGroup<T> g) where T : struct, IElt<T>
    {
        if (g.GroupType == GroupType.NonAbelianGroup)
            throw new GroupException(GroupExceptionType.OnlyAbelianGroups);

        var ne = g.Neutral();
        if (g.Count() == 1)
            return (new[] { (ne, 1) }, new[] { (ne, 1) });

        var set = g.GetGenerators().SelectMany(e => GensToFactors(g, e)).OrderBy(e => e.Item2).ToList();
        var facts = set.OrderBy(e => e.Item2).ToArray();
        while (true)
        {
            var arr = set.OrderByDescending(a => a.Item2).ToArray();
            var (e0, e1) = set.OrderByDescending(a => a.Item2)
                .Select(e => (e, arr.FirstOrDefault(c => IntExt.Gcd(e.Item2, c.Item2) == 1, (ne, -1))))
                .FirstOrDefault(a => a.Item2.Item2 != -1, ((ne, -1), (ne, -1)));
            if (e0.Item2 == -1)
                break;

            set.Remove(e0);
            set.Remove(e1);
            var e2 = g.Op(e0.Item1, e1.Item1);
            var o2 = e0.Item2 * e1.Item2;
            if (g.ElementsOrders[e2] != o2)
                throw new($"{e0} & {e1} => {(e2, o2)} | {g.ElementsOrders[e2]}");

            set.Add((e2, o2));
        }

        return (set.OrderByDescending(a => a.Item2).ToArray(), facts);
    }

    static IEnumerable<T> CycleExceptNeutral<T>(ConcreteGroup<T> g, T e) where T : struct, IElt<T>
    {
        var a = g.Neutral();
        var o = g.ElementsOrders[e];
        for (int i = 0; i < o - 1; i++)
        {
            a = g.Op(a, e);
            yield return a;
        }
    }

    // Give the cardinal of elements with order p^i for Group (Z/(p^i)Z )^k of order p^(ik) 
    static int PMagik(int p, int i, int k) => p.Pow((i - 1) * k) * (p.Pow(k) - 1); // TODO proof

    static int[] MagikFactors<T>(ConcreteGroup<T> g) where T : struct, IElt<T>
    {
        if (g.GroupType == GroupType.NonAbelianGroup)
            throw new GroupException(GroupExceptionType.OnlyAbelianGroups);
        
        var ord = g.Count();
        if (ord == 1)
            return new [] { 1 };
        
        var dec = IntExt.PrimesDec(ord);
        var lt = new List<(int, int, int)>();
        var elemsOrds = g.ElementsOrders.Values.GroupBy(e => e).ToDictionary(e => e.Key, e => e.Count());
        foreach (var (p, n) in dec)
        {
            var rg = n.Range(1).OrderDescending().ToArray();
            foreach (var i in rg)
            {
                var q = p.Pow(i);
                if (!elemsOrds.ContainsKey(q))
                    continue;

                foreach (var k in rg.Where(k => k <= n - i + 1))
                {
                    var mgk = PMagik(p, i, k);
                    if (elemsOrds[q] % mgk == 0)
                    {
                        lt.Add((p, i, k));
                        break;
                    }
                }
            }
        }

        var dim = lt.Max(e => e.Item3);
        var can = Enumerable.Repeat(1, dim).ToArray();
        foreach (var (p, i, k) in lt)
        {
            var q = p.Pow(i);
            for (int j = 0; j < k; j++)
            {
                if (can[j] % q != 0)
                    can[j] *= q;
            }
        }
        
        return can;
    }

    static List<(T, int)> MagikDecomposition<T>(ConcreteGroup<T> g, int[] abType) where T : struct, IElt<T>
    {
        var elemsOrds = g.ElementsOrders;
        var set = g.ToHashSet();
        var tmpSet = new HashSet<T>() { g.Neutral() };
        var gens = new List<(T, int)>();
        foreach (var o in abType)
        {
            var e = set.First(e => elemsOrds[e] == o && !tmpSet.Overlaps(CycleExceptNeutral(g, e)));
            tmpSet = Group.GenerateElements(g, tmpSet, new() { e });
            set.ExceptWith(tmpSet);
            gens.Add((e, o));
        }

        return gens;
    }

    static void ShowMagikDecomposition<T>(ConcreteGroup<T> g) where T : struct, IElt<T>
    {
        GlobalStopWatch.Restart();
        GlobalStopWatch.AddLap();
        var lt = MagikFactors(g);
        GlobalStopWatch.Show($"Magik Factors {g}");

        GlobalStopWatch.AddLap();
        var gens = MagikDecomposition(g, lt);
        var can0 = gens.Select(e => e.Item2).Descending().ToArray();
        GlobalStopWatch.Show($"Canonic Generators {g}");

        var can = AbelianInvariantsFactors.Reduce(g).OrderDescending().ToArray();
        lt.Println($"{g}  ~  {can.Glue(" x ", "C{0}")}");
        if (!can0.SequenceEqual(can))
            throw new();

        gens.Println("Gens");
        var g1 = Group.Generate(can0.Glue(" x ", "C{0}"), g, gens.Select(e => e.Item1).ToArray());
        Console.WriteLine($"{g.ShortName}  =>  {g1.ShortName}");
        if (g.Count() != g1.Count())
            throw new();

        Console.WriteLine();
    }

    static string ArrToString(int[] seq) => seq.Glue(" x ", "C{0}");

    static void ShowAbTypes(params int[] seq)
    {
        var (can, facts) = AbType(seq);
        Console.WriteLine("{0}  ~  {1}  ~  {2}", ArrToString(seq), ArrToString(can), ArrToString(facts));
        Console.WriteLine();
    }

    static void AbGroupOrders(params int[] seq)
    {
        var s = ArrToString(seq);

        GlobalStopWatch.Restart();
        var g = FG.Abelian(seq);
        var ords1 = g.ElementsOrders.Values.GroupBy(e => e).Select(e => (o: e.Key, nb: e.Count())).Order().ToArray();
        GlobalStopWatch.Show($"{s} Old");

        GlobalStopWatch.Restart();
        var facts = SeqToFactors(seq);
        var ords2 = AbGroupElemsOrders(facts).Order().ToArray();
        GlobalStopWatch.Show($"{s} New");

        var check = ords1.SequenceEqual(ords2);
        Console.WriteLine(check ? "PASS" : "FAIL");
        if (!check)
        {
            ords1.Println("Expected");
            ords2.Println("Actual");
        }

        Console.WriteLine();
    }

    static void ShowAllAb(int n)
    {
        var comp = Comparer<int[]>.Create((a, b) => a.SequenceCompareTo(b));
        var tab0 = AllAbOfOrder(n).Select(e => (e, FactorsToType(e)))
            .OrderBy(e => e.Item2.Length)
            .ThenBy(e => e.Item2, comp)
            .ToArray();
        var tab = tab0.Select(e => new[] { ArrToString(e.e), ArrToString(e.Item2) }).ToArray();

        var title = $"All Abelian Groups of Order:{n}";
        var max0 = int.Max(title.Length, tab.Max(e => e[0].Length));
        var max1 = tab.Max(e => e[1].Length);
        var fmt = $"{{0,{max0}}}  ~  {{1,-{max1}}}";

        Console.WriteLine($"{{0,{max0}}}     {{1,-{max1}}}", title, $"Count:{tab.Length} / {AllAbCount(n)}");
        foreach (var s in tab)
            Console.WriteLine(fmt, s[0], s[1]);

        Console.WriteLine();
    }

    static void ShowAbInvariants<T>(ConcreteGroup<T> g) where T : struct, IElt<T>
    {
        if (g.GroupType == GroupType.NonAbelianGroup)
            throw new GroupException(GroupExceptionType.OnlyAbelianGroups);

        DisplayGroup.Head(g);
        var gens = g.GetGenerators().ToDictionary(e => e, e => g.ElementsOrders[e]);

        Console.WriteLine("Gens: {0}", ArrToString(gens.Values.Order().ToArray()));
        Console.WriteLine("Expected");
        Console.WriteLine(ArrToString(SeqToFactors(gens.Values.ToArray())));
        var abType = AbType(gens.Values.ToArray());
        Console.WriteLine(ArrToString(abType.can));

        Console.WriteLine("Actual");
        GlobalStopWatch.Restart();
        var abInvs = AbInvariants(g);
        GlobalStopWatch.Show();
        Console.WriteLine(ArrToString(abInvs.facts.Select(e => e.Item2).ToArray()));
        Console.WriteLine(ArrToString(abInvs.can.Select(e => e.Item2).ToArray()));
        if (!abType.can.Order().SequenceEqual(abInvs.can.Select(e => e.Item2).Order()))
            throw new();

        Console.WriteLine();
    }

    static void ShowPMagik(int p, int maxOrder = 10000)
    {
        var n = (int)(double.Log10(maxOrder) / double.Log10(p));
        for (int i = 1; i <= n; i++)
        {
            var q = p.Pow(i);
            var kmax = n / i;
            for (int k = 1; k <= kmax; k++)
            {
                var seq = Enumerable.Repeat(q, k).ToArray();
                var e = AbGroupElemsOrders(seq).First(e => e.o == q);
                var nb = p.Pow((i - 1) * k) * (p.Pow(k) - 1);
                var name = seq.Glue(" x ", "C{0}");
                Console.WriteLine($"G={name} ord(G)={q.Pow(k)}    Card({{e in G, ord(e)={q}}})={e.nb}={nb}");
                if (e.nb != nb)
                    throw new($"Expected:{e.nb} Actual:{nb}");
            }

            Console.WriteLine();
        }
    }

    public static void Example1Integers()
    {
        ShowAbTypes(1);
        ShowAbTypes(41);
        ShowAbTypes(299);
        ShowAbTypes(14, 21);
        ShowAbTypes(20, 30);
        ShowAbTypes(8, 18, 30);
        ShowAbTypes(24, 18, 30);
    }

    public static void Example2Integers()
    {
        ShowAllAb(1);
        ShowAllAb(13);
        ShowAllAb(18);
        ShowAllAb(20);
        ShowAllAb(16);
        ShowAllAb(24);
        ShowAllAb(32);
        ShowAllAb(36);
        ShowAllAb(64);
        ShowAllAb(299);
        ShowAllAb(360);
        ShowAllAb(1024);
        ShowAllAb(4320);
        ShowAllAb(4096);
    }

    public static void Example3Integers()
    {
        AbGroupOrders(97);
        AbGroupOrders(143);
        AbGroupOrders(196);
        AbGroupOrders(196);
        AbGroupOrders(256);
        AbGroupOrders(14, 21);
        AbGroupOrders(20, 30);
        AbGroupOrders(12, 36);
        AbGroupOrders(4, 12, 36);
        AbGroupOrders(8, 16, 4, 16);
        AbGroupOrders(8, 18, 30);
        AbGroupOrders(2, 6, 360);
        AbGroupOrders(6, 12, 360);
    }

    public static void FastAbelianInvariants()
    {
        ShowAbInvariants(FG.Abelian(1));
        ShowAbInvariants(FG.Abelian(17));
        ShowAbInvariants(FG.Abelian(7, 11));
        ShowAbInvariants(FG.Abelian(14, 21));
        ShowAbInvariants(FG.Abelian(20, 30));
        ShowAbInvariants(Product.Generate(FG.Abelian(2, 9), FG.Abelian(6, 4)));
        ShowAbInvariants(FG.Abelian(2, 9, 6, 4));
        ShowAbInvariants(FG.Abelian(8, 18, 30));
        ShowAbInvariants(Product.Generate(FG.Abelian(24), FG.Abelian(18, 30)));
        ShowAbInvariants(FG.Abelian(48, 64));
        ShowAbInvariants(Product.Generate(FG.Abelian(6, 8), FG.Abelian(8, 8)));

    }

    public static void FastAbelianInvariantsErrors()
    {
        var g0 = FG.Abelian(2, 16);
        var g1 = Group.Generate("G", g0, g0[0, 1], g0[1, 2]);
        ShowAbInvariants(g1);
    }
    
    public static void TestPMagik()
    {
        ShowPMagik(2, maxOrder: 20000);
        ShowPMagik(3, maxOrder: 20000);
        ShowPMagik(5, maxOrder: 20000);
        ShowPMagik(7, maxOrder: 20000);
    }

    public static void ExamplePMagik()
    {
        ShowMagikDecomposition(FG.Abelian(14, 21));
        ShowMagikDecomposition(FG.Abelian(20, 30));
        ShowMagikDecomposition(FG.Abelian(8, 18, 30));
        ShowMagikDecomposition(FG.Abelian(8, 18, 30));
        ShowMagikDecomposition(FG.Abelian(8, 18, 30));

        foreach (var g in FG.AllAbelianGroupsOfOrder(324))
            ShowMagikDecomposition(g);
        
        foreach (var g in FG.AllAbelianGroupsOfOrder(576))
            ShowMagikDecomposition(g);

        foreach (var g in FG.AllAbelianGroupsOfOrder(750))
            ShowMagikDecomposition(g);

        foreach (var g in FG.AllAbelianGroupsOfOrder(810))
            ShowMagikDecomposition(g);
        
        foreach (var g in FG.AllAbelianGroupsOfOrder(1024))
            ShowMagikDecomposition(g);

        foreach (var g in FG.AllAbelianGroupsOfOrder(1458))
            ShowMagikDecomposition(g);

        var g0 = FG.Abelian(2, 16);
        var g1 = Group.Generate("G", g0, g0[0, 1], g0[1, 2]);
        ShowMagikDecomposition(g1);
    }

    public static void UnIsomorphisms()
    {
        var seqUn = 120.Range(3).Select(i => (i, un: new Un(i))).ToArray();
        foreach (var (i, un) in seqUn)
            Console.WriteLine("U{0,-3} ~ {1}", i, MagikFactors(un).Glue(" x ", "C{0}"));

        void NewMeth((int i, Un un)[] seq)
        {
            var arr = seq.Select(e => MagikFactors(e.un)).ToArray();
        }
        void OldMeth((int i, Un un)[] seq)
        {
            var arr = seq.Select(e => AbelianInvariantsFactors.Reduce(e.un)).ToArray();
        }
        
        GlobalStopWatch.Bench(5, "NewDecomp", () => NewMeth(seqUn));
        GlobalStopWatch.Bench(5, "OldDecomp", () => OldMeth(seqUn));
    }
}

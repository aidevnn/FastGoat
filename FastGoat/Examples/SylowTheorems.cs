using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Integers;

namespace FastGoat.Examples;

public static class SylowTheorems
{
    static void SylowFirstTheorem<T>(ConcreteGroup<T> g) where T : struct, IElt<T>
    {
        DisplayGroup.Head(g);
        var og = g.Count();
        var dec = IntExt.PrimesDec(og);
        var allPr = dec.SelectMany(e => e.Value.Range(1).Select(i => e.Key.Pow(i))).ToHashSet();
        var tablePSubGr = Group.AllPSubGroups(g);
        var allPSubGr = tablePSubGr.Values.SelectMany(e0 => e0).ToHashSet(new GroupSetEquality<T>());
        foreach (var (g0, conjs) in tablePSubGr.OrderBy(e => e.Key.ShortName))
            Console.WriteLine(
                $"{g0.ShortName} {g0.GroupType} Elements Orders {{{g0.ElementsOrders.Values.ToHashSet().Order().Glue(" ")}}} Conjugates {conjs.Count}");

        var allSubGr = Group.AllSubGroups(g);
        var tablePSubGr0 = allSubGr.Where(e => IntExt.PrimesDec(e.Key.Count()).Count == 1).ToDictionary(e => e.Key, e => e.Value);
        var allPr0 = tablePSubGr.Keys.Select(e => e.Count()).ToHashSet();
        var allPSubGr0 = tablePSubGr0.Values.SelectMany(e0 => e0).ToArray();
        var check = allPr0.SetEquals(allPr) && tablePSubGr.Count == tablePSubGr0.Count && allPSubGr.SetEquals(allPSubGr0);

        Console.WriteLine($"All p^r : expected [{allPr.Glue(" ")}] actual [{allPr0.Glue(" ")}]");
        Console.WriteLine($"########### Check {g} : {(check ? "PASS" : "FAIL")} ###########");
        Console.WriteLine();
    }

    static void SylowSecondTheorem<T>(ConcreteGroup<T> g) where T : struct, IElt<T>
    {
        DisplayGroup.Head(g);
        var og = g.Count();
        var tablePSubGr = Group.AllPSubGroups(g);
        var tablePSubGr0 = tablePSubGr.GroupBy(e => IntExt.PrimesDec(e.Key.Count()).First().Key)
            .ToDictionary(e => e.Key, e => e.ToDictionary(f => f.Key, f => f.Value));

        foreach (var (p, table) in tablePSubGr0)
        {
            if (!IntExt.Primes10000.Contains(p))
                throw new();
            var (sylow, conjs) = table.MaxBy(e => e.Key.Count());
            table.Remove(sylow);
            var np = conjs.Count;
            Console.WriteLine($"{sylow.ShortName} {sylow.GroupType} Conjugates np = {np}");
            Console.WriteLine($"np = {1} mod {p} {np % p == 1} and np | {og} {og % np == 0}");
            if (np == 1)
                DisplayGroup.Head(g.Over(sylow));

            foreach (var (sg, _) in table)
                Console.WriteLine($"    {sg.ShortName} {sg.GroupType} is subset of {sylow} {sg.SubSetOf(sylow)}");

            Console.WriteLine();
        }
    }

    static ConcreteGroup<Ep<ZnInt>> AbelianElementariesFactors<T>(ConcreteGroup<T> g) where T : struct, IElt<T>
    {
        if (g.GroupType == GroupType.NonAbelianGroup)
            throw new();

        var Sylows = Group.AllSylowPSubgroups(g).Keys;
        var elems = new List<int>();
        Console.WriteLine($"########################### {g} ###########################");
        DisplayGroup.Head(g);
        Console.WriteLine("Sylows subgroups");
        foreach (var sgSylow in Sylows)
        {
            DisplayGroup.Head(sgSylow);
            var lt = AbelianInvariantsFactors.Reduce(sgSylow);
            var ab0 = FG.Abelian(lt.Order().ToArray());
            DisplayGroup.AreIsomorphics(ab0, sgSylow);
            elems.AddRange(lt);
            Console.WriteLine();
        }

        var abElems = FG.Abelian(elems.Order().ToArray());
        abElems.SetName($"({abElems.Name})elem");
        var invFacts = AbelianInvariantsFactors.Reduce(g).Order().ToArray();
        var abInv = FG.Abelian(invFacts);
        abInv.SetName($"({abInv.Name})inv");
        
        DisplayGroup.AreIsomorphics(abElems, g);
        DisplayGroup.AreIsomorphics(abInv, g);
        var elems2 = invFacts.SelectMany(f => IntExt.PrimesDec(f).Select(kv => kv.Key.Pow(kv.Value))).Order().ToArray();
        var abElems2 = FG.Abelian(elems2.Order().ToArray());
        abElems2.SetName($"({abElems2.Name})elem2");
        Console.WriteLine($"{abElems} set equal {abElems2} : {abElems2.SetEquals(abElems)}");
        Console.WriteLine();
        return abElems;
    }

    // H.E.Rose, 10.2 Frattini and Fitting Subgroups page221
    static void FittingSubgroupProperties<T>(ConcreteGroup<T> g) where T : struct, IElt<T>
    {
        var tableSubgroups = Group.AllSubGroups(g);
        var frat = Group.FrattiniSubGroup(tableSubgroups);
        var tableSylows = Group.AllSylowPSubgroups(tableSubgroups);
        var allNormalSubgroups = Group.AllNormalSubgroups(tableSubgroups);

        ConcreteGroup<T> PRadical(List<ConcreteGroup<T>> sylows)
        {
            return Group.Generate(g, sylows.Aggregate(g.ToArray(), (acc, g0) => acc.Intersect(g0).ToArray()));
        }

        // Nilpotency definition
        var fitting1 = Group.Generate($"F({g})", g,
            allNormalSubgroups.Where(n => Group.IsNilpotent(n)).SelectMany(g0 => g0.GetGenerators()).ToArray());
        
        // PSylows property
        var fitting2 = Group.Generate($"F({g})2", g, tableSylows.SelectMany(kv => PRadical(kv.Value).GetGenerators()).ToArray());
        
        DisplayGroup.Head(g);
        DisplayGroup.Head(frat);
        DisplayGroup.Head(fitting1);
        if (fitting1.GroupType == GroupType.AbelianGroup)
        {
            var ab = FG.Abelian(AbelianInvariantsFactors.Reduce(fitting1).OrderDescending().ToArray());
            if (ab.Count() > 1)
                DisplayGroup.AreIsomorphics(fitting1, ab);
        }
        else if (Group.IsNilpotent(g))
        {
            DisplayGroup.AreIsomorphics(fitting1, g);
        }
        
        DisplayGroup.AreIsomorphics(fitting1, fitting2);
        Console.WriteLine($"{frat} is subgroup of {fitting1} {frat.SubSetOf(fitting1)}");
        Console.WriteLine();
    }

    public static void FirstTheoremExamples()
    {
        SylowFirstTheorem(FG.Dihedral(4));
        SylowFirstTheorem(FG.Symmetric(4));
        SylowFirstTheorem(FG.Alternate(5));
        SylowFirstTheorem(FG.Abelian(4, 9));
        SylowFirstTheorem(FG.Abelian(6, 8));
        SylowFirstTheorem(Group.SemiDirectProd(new Cn(5), new Cn(4)));
    }

    public static void SecondTheoremExamples()
    {
        SylowSecondTheorem(FG.Dihedral(4));
        SylowSecondTheorem(FG.Symmetric(4));
        SylowSecondTheorem(FG.Alternate(5));
        SylowSecondTheorem(FG.Abelian(4, 9));
        SylowSecondTheorem(FG.Abelian(6, 8));
        SylowSecondTheorem(Group.SemiDirectProd(new Cn(5), new Cn(4)));
    }

    public static void AbelianElementariesFactorsExamples()
    {
        AbelianElementariesFactors(FG.Abelian(2, 4, 6));
        AbelianElementariesFactors(FG.Abelian(2, 3, 8));
        AbelianElementariesFactors(FG.Abelian(3, 4, 4));
        AbelianElementariesFactors(FG.Abelian(6, 8, 18));
        AbelianElementariesFactors(FG.Abelian(8, 18, 30));
    }

    public static void PQGroups()
    {
        var p20 = IntExt.Primes10000.Take(10).ToArray();
        var allPQ = p20.Grid2D(p20).Select(e => (p: e.t1, q: e.t2)).Where(e => e.p != e.q && e.q % e.p != 1).ToArray();
        var abelianGroups = allPQ.Select(e => (e, FG.Abelian(e.p, e.q))).OrderBy(e => e.Item2.Count()).ToArray();

        foreach (var ((p, q), g0) in abelianGroups)
        {
            var Sylows = Group.AllSylowPSubgroups(g0);
            var (s0, conjs) = Sylows.First(e => IntExt.PrimesDecomposition(e.Key.Count()).All(k => k == p));
            if (conjs.Count != 1)
                throw new();

            Console.WriteLine($"{g0} Sylow {s0.ShortName}");
        }

        var sdpGroups = allPQ
            .SelectMany(e => Group.AllSemiDirectProd($"C{e.p} x: C{e.q}", new Cn(e.p), new Cn(e.q)).Select(g0 => (e, g0)))
            .ToArray();
        foreach (var ((p, q), g0) in sdpGroups)
        {
            var Sylows = Group.AllSylowPSubgroups(g0);
            var (s0, conjs) = Sylows.First(e => IntExt.PrimesDecomposition(e.Key.Count()).All(k => k == p));
            if (conjs.Count != 1)
                throw new();

            Console.WriteLine($"{g0} Sylow {s0.ShortName}");
        }

        Console.WriteLine("############## PASS ##############");
    }

    public static void FittingSubgroupExample()
    {
        FittingSubgroupProperties(FG.Symmetric(3));
        FittingSubgroupProperties(FG.Symmetric(4));
        FittingSubgroupProperties(FG.Symmetric(5));
        FittingSubgroupProperties(FG.Alternate(5));
        FittingSubgroupProperties(new Cn(32));
        FittingSubgroupProperties(FG.Dihedral(6));
        
        FittingSubgroupProperties(FG.Dihedral(4));
        FittingSubgroupProperties(FG.Dihedral(5));
        FittingSubgroupProperties(FG.DiCyclic(4));
        FittingSubgroupProperties(FG.DiCyclic(5));
    }
}

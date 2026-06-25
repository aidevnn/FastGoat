using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.Subgroups;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Perms;
using FastGoat.UserGroup.Words;

namespace Examples;

public static class HolomorphGroup
{
    public static SemiDirectProduct<T, Automorphism<T>> Holomorph<T>(ConcreteGroup<T> G) where T : struct, IElt<T>
    {
        var autG = Group.AutomorphismGroup(G);
        var theta = Group.Hom(autG,
            Group.HomomorphismMap(autG, autG, autG.GetGenerators().ToDictionary(e => e, e => e)));
        return Group.SemiDirectProd($"Hol[{G}]", G, theta, autG);
    }

    public static (ConcreteGroup<Perm> G, ConcreteGroup<Perm> AutG, ConcreteGroup<Perm> HolG)
        HolomorphPerm<T>(ConcreteGroup<T> G) where T : struct, IElt<T>
    {
        var (G1, AutG, _) = FG.RegPermAutGroup(G);
        var HolG = Group.DirectProduct($"Hol[{G}]", G1, AutG);
        return (G1, AutG, HolG);
    }

    static (Perm cnD2n, Perm c2D2n, Perm cnAutD2n, Perm[] phiAutD2n, Perm c2nHolD2n) HolD2nPerm(int n)
    {
        var pType = IntExt.PrimesDec(n).Select(e => e.Key.Pow(e.Value)).ToArray();
        var a0 = IntExt.PermAndCyclesFromType(pType);
        var n1 = a0.perm.Length;
        var sn = new Sn(2 * n1);
        var bCycles = a0.cycles.SelectMany(e => e.Zip(e.Select(i => i + n1).Reverse().ToArray())).ToArray();
        var cCycles = a0.cycles.Select(o => o.SelectMany(i => new[] { i, i + n1 }).ToArray()).ToArray();

        var cnD2n = FG.ConcatPerm(FG.Cycles(n), FG.Cycles(n));
        var c2D2n = sn.OpSeq(bCycles.Select(c => sn.CycleP1([c.First, c.Second])));
        var c2nHolD2n = sn.OpSeq(cCycles.Select(c => sn.CycleP1(c)));

        var cnAutD2n = FG.PaddingRight(FG.Cycles(n), n1);
        var phiAutD2n = FG.AutomorphismDihedralGens(n).gensUn.Select(e => FG.ConcatPerm(e, e)).ToArray();

        return (cnD2n, c2D2n, cnAutD2n, phiAutD2n, c2nHolD2n);
    }

    public static void Example1HolC7()
    {
        // Saunders MacLane, Garrett Birkhoff. Algebra (3rd ed.) page 416
        var c7 = new Cn(7);
        var c6 = new Cn(6);
        var thetas = Group.AllOpsByAutomorphisms(c6, c7);
        List<SemiDirectProduct<ZnInt, ZnInt>> nonIsomorphics = new();
        char k = 'a';
        var nm = "(C7 x: C6)";
        foreach (var theta in thetas)
        {
            var sdp = Group.SemiDirectProd($"{nm}{k++}", c7, theta, c6);
            if (nonIsomorphics.All(g => !g.IsIsomorphicTo(sdp)))
                nonIsomorphics.Add(sdp);
        }

        nonIsomorphics.Sort((a, b) => -a.ElementsOrdersList().Max().CompareTo(b.ElementsOrdersList().Max()));
        foreach (var sdp in nonIsomorphics)
        {
            DisplayGroup.HeadSdp(sdp);
            Console.WriteLine(sdp.ElementsOrdersList().Glue(", "));
            Console.WriteLine();
        }

        var c7c6 = Product.Generate(c6, c7);
        DisplayGroup.Head(c7c6);
        Console.WriteLine(c7c6.ElementsOrdersList().Glue(", "));
        Console.WriteLine();

        var wgC3D14 = new WordGroup("C3 x D14", "a7, b2, c3, abab, ac = ca, bc = cb");
        DisplayGroup.Head(wgC3D14);
        Console.WriteLine(wgC3D14.ElementsOrdersList().Glue(", "));
        Console.WriteLine();

        var wgC2H21 = new WordGroup("C2 x H21", "a7, b3, c2, a2 = bab-1, ac = ca, bc = cb");
        DisplayGroup.Head(wgC2H21);
        Console.WriteLine(wgC2H21.ElementsOrdersList().Glue(", "));
        Console.WriteLine();

        var u7 = new Un(7);
        var a = u7[(c7[1], c7[3])];
        var theta1 = Group.HomomorphismMap(u7, u7, new() { [a] = a });
        var homTheta1 = Group.Hom(u7, theta1);
        var hol7 = Group.SemiDirectProd("Hol7", c7, homTheta1, u7);
        DisplayGroup.HeadSdp(hol7);
        Console.WriteLine(hol7.ElementsOrdersList().Glue(", "));
        Console.WriteLine();

        var g0 = nonIsomorphics[0];
        Console.WriteLine("{0} IsIsomorphicTo {1} : {2}", c7c6, g0, c7c6.IsIsomorphicTo(g0));
        var g1 = nonIsomorphics[1];
        Console.WriteLine("{0} IsIsomorphicTo {1} : {2}", wgC3D14, g1, wgC3D14.IsIsomorphicTo(g1));
        var g2 = nonIsomorphics[2];
        Console.WriteLine("{0} IsIsomorphicTo {1} : {2}", wgC2H21, g2, wgC2H21.IsIsomorphicTo(g2));
        var g3 = nonIsomorphics[3];
        Console.WriteLine("{0} IsIsomorphicTo {1} : {2}", hol7, g3, hol7.IsIsomorphicTo(g3));
    }

    public static void Example2CyclicGroup()
    {
        DisplayGroup.HeadOrdersNames(Holomorph(new Cn(4)), setName: false);
        DisplayGroup.HeadOrdersNames(Holomorph(FG.Abelian(2, 2)), setName: false);
        DisplayGroup.HeadOrdersNames(Holomorph(new Cn(5)), setName: false);
        DisplayGroup.HeadOrdersNames(Holomorph(new Cn(6)), setName: false);
        DisplayGroup.HeadOrdersNames(Holomorph(new Cn(7)), setName: false);
        DisplayGroup.HeadOrdersNames(Holomorph(new Cn(8)), setName: false);
        DisplayGroup.HeadOrdersNames(Holomorph(new Cn(9)), setName: false);
        DisplayGroup.HeadOrdersNames(Holomorph(FG.Abelian(3, 3)), setName: false);
        DisplayGroup.HeadOrdersNames(Holomorph(new Cn(11)), setName: false);
    }

    public static void Example3DihedralAndDicyclic()
    {
        DisplayGroup.HeadOrdersNames(Holomorph(FG.DihedralWg(4)), setName: false);
        DisplayGroup.HeadOrdersNames(Holomorph(FG.QuaternionWg(8)), setName: false);
        DisplayGroup.HeadOrdersNames(Holomorph(FG.DihedralWg(5)), setName: false);
        DisplayGroup.HeadOrdersNames(Holomorph(FG.DiCyclic(3)), setName: false);
        DisplayGroup.HeadOrdersNames(Holomorph(FG.DihedralWg(6)), setName: false);
        DisplayGroup.HeadOrdersNames(Holomorph(FG.DihedralWg(7)), setName: false);
    }

    public static void Example4SimpleCyclicGroup()
    {
        foreach (var p in IntExt.Primes10000.Skip(1).Take(10))
        {
            var holCp = Holomorph(new Cn(p));
            var Fpq = FG.MetaCyclicPg(p, p - 1, IntExt.Solve_k_pow_m_equal_one_mod_n_strict(p, p - 1));
            var autD2p = Group.AutomorphismGroup(FG.Dihedral(p));
            DisplayGroup.Head(Fpq);
            DisplayGroup.AreIsomorphics(Fpq, holCp);
            DisplayGroup.AreIsomorphics(autD2p, holCp);
            Console.WriteLine();
        }
    }

    public static void Example5SimpleCyclicGroupPerm()
    {
        foreach (var p in IntExt.Primes10000.Skip(1).Take(10))
        {
            var (Cp, AutCp, holCp) = HolomorphPerm(new Cn(p));
            var Fpq = FG.MetaCyclicPg(p, p - 1, IntExt.Solve_k_pow_m_equal_one_mod_n_strict(p, p - 1));
            var (D2p, autD2p, _) = FG.RegPermAutGroup(FG.Dihedral(p));
            DisplayGroup.AreIsomorphics(FG.AbelianPerm(p - 1), AutCp);
            DisplayGroup.AreIsomorphics(Fpq, holCp);
            DisplayGroup.AreIsomorphics(autD2p, holCp);
            Console.WriteLine();
        }
    }

    public static void Example6Dihedral()
    {
        GlobalStopWatch.Restart();

        for (int n = 3; n <= 32; ++n)
        {
            var (cnD2n, c2D2n, cnAutD2n, phiAutD2n, c2nHolD2n) = HolD2nPerm(n);
            var holGens = IntExt.IsPrime(n)
                ? phiAutD2n.Prepend(c2nHolD2n).ToArray() // ord(g1) = 2n and gen of Un ord(g2) = n - 1
                : phiAutD2n.Append(c2D2n).Prepend(c2nHolD2n).ToArray(); // ord(g1) = 2n, ord(g2) = 2 and gens of Un 

            var d2n = Group.Generate($"D{2 * n}", cnD2n.Sn, cnD2n, c2D2n);
            var autD2n = Group.Generate($"Aut[{d2n}]", cnD2n.Sn, phiAutD2n.Append(cnAutD2n).ToArray());
            var holD2n = Group.Generate($"Hol[{d2n}]", cnD2n.Sn, holGens);
            DisplayGroup.HeadOrdersGenerators(d2n);
            DisplayGroup.HeadOrdersGenerators(autD2n);
            DisplayGroup.HeadOrdersGenerators(holD2n);
            Console.WriteLine($"{autD2n} = C{n} x: {phiAutD2n.Select(e => e.Order).ToAbString().WithParenthesis()}");
            Console.WriteLine($"{holD2n} = {d2n} x: {autD2n}");
            Console.WriteLine();

            var d2nNormal = Group.IsNormalSubgroup(holD2n, d2n);
            var autD2nNotNormal = Group.IsNormalSubgroup(holD2n, autD2n);
            var inter = d2n.Intersect(autD2n).ToXSet();
            Console.WriteLine($"{d2n,-10} is normal subgroup of {holD2n,-10} {d2nNormal}");
            Console.WriteLine($"{autD2n,-10} is normal subgroup of {holD2n,-10} {autD2nNotNormal}");
            Console.WriteLine($"{d2n,-10} ∩ {autD2n,10} = {inter}");
            Console.WriteLine();
            if (!d2nNormal || autD2nNotNormal || inter.Count != 1)
                throw new();
        }

        GlobalStopWatch.Show();
    }
}
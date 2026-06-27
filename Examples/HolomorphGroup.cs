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

    static (string name, Perm[] dicm, Perm[] aut, Perm[] hol) HolDicmPerm(int m)
    {
        if (int.IsPow2(m))
        {
            var (c2ma, _, c2mb, phi, c4m) = HolD2nPerm(2 * m);
            var Q4m = FG.QuaternionPg(4 * m);
            var c4 = Q4m.GetGenerators().First(e => !e.Equals(c2ma));
            var gensAut = phi.Append(c2mb).ToArray();
            var gensHol = phi.Append(c4).Prepend(c4m).ToArray();
            if (m == 2)
            {
                gensAut = new[] { c2mb, c2ma.Sn[(2, 6), (4, 8), (1, 3)] };
                gensHol = gensAut.Append(c4m).ToArray();
            }

            return (Q4m.Name, [c2ma, c4], gensAut, gensHol);
        }
        else if (m % 2 == 1)
        {
            var (c4a, c2a) = FG.MetaCyclicGens(4, 2, 3).gens.Deconstruct();
            var (cma, c2b, cmb, phi, _) = HolD2nPerm(m);
            var cm = FG.PaddingRight(cma, 4);
            var c4 = FG.ConcatPerm(c2b, c4a);
            var cmAut = FG.ConcatPerm(cmb, c2a);
            var phiDic = phi.Select(e => FG.ConcatPerm(e, c2a)).ToArray();
            var gensAut = phiDic.Append(cmAut).ToArray();
            var gensHol = phiDic.Append(cmAut).Append(c4).ToArray();
            return ($"Dic{m}", [cm, c4], gensAut, gensHol);
        }
        else
        {
            var pow2 = IntExt.PrimesDecomposition(m).Count(i => i == 2);
            var n = m / (1 << pow2);
            var q = 1 << pow2;
            var name = $"C{n} x: Q{4 * q}";

            var (cna, c2a, cnb, phi_n, c2n) = HolD2nPerm(n);
            var gensAutD2n = phi_n.Prepend(cnb).ToArray();
            var (c2qa, c2b, c2qb, phi_2q, c4q) = HolD2nPerm(2 * q);
            var gensAutQ4q = phi_2q.Prepend(c2qb).ToArray();
            var Q4q = FG.QuaternionPg(4 * q);
            var c4 = Q4q.GetGenerators().First(e => !e.Equals(c2qa));

            var cn = FG.PaddingRight(cna, c2qa.Sn.N);
            var hom = phi_2q.ToDictionary(e => e, _ => c2a.Sn.Neutral());
            hom[c2qb] = Group.GenerateElements(cna.Sn, phi_n).First(e => e.Order == 2);

            var gensAutD2nb = gensAutD2n.Select(c => FG.PaddingRight(c, c2qa.Sn.N)).ToArray();
            var gensAutQ4qb = gensAutQ4q.Select(c => FG.ConcatPerm(hom[c], c)).ToArray();
            var gensAut = gensAutD2nb.Concat(gensAutQ4qb).ToArray();
            var gensHol = gensAut.Append(FG.ConcatPerm(c2n, c4q)).ToArray();
            return (name, [cn, FG.PaddingLeft(c2qa, cna.Sn.N), FG.ConcatPerm(c2a, c4)], gensAut, gensHol);
        }
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

        for (int n = 3; n <= 48; ++n)
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
            var holOrd = d2n.Order * autD2n.Order;
            Console.WriteLine($"{($"|{holD2n}|"),-10} = {($"|{d2n}| * |{autD2n}|"),-30} {holD2n.Order == holOrd}");
            Console.WriteLine($"{d2n,-10} is normal subgroup of {holD2n,-10} {d2nNormal}");
            Console.WriteLine($"{autD2n,-10} is normal subgroup of {holD2n,-10} {autD2nNotNormal}");
            Console.WriteLine($"{d2n,-10} ∩ {autD2n,10} = {inter}");
            Console.WriteLine();
            if (!d2nNormal || autD2nNotNormal || inter.Count != 1 || holD2n.Order != holOrd)
                throw new();
            
            if (n <= 24)
            {
                var (d2nReg, autReg, holReg) = HolomorphPerm(FG.Dihedral(n));
                if (!d2nReg.ElementsOrdersList().SequenceEqual(d2n.ElementsOrdersList()) ||
                    !autReg.ElementsOrdersList().SequenceEqual(autD2n.ElementsOrdersList()) ||
                    !holReg.ElementsOrdersList().SequenceEqual(holD2n.ElementsOrdersList()))
                {
                    DisplayGroup.HeadOrdersGenerators(holReg);
                    throw new();
                }
            }
        }

        GlobalStopWatch.Show(); // Time:32.050s
    }
    // |Hol[D90]| = 97200
    // Type        NonAbelianGroup
    // BaseGroup   S28
    // 
    // Elements Orders : [1]:1, [2]:2311, [3]:26, [4]:5000, [5]:24, [6]:6302, [9]:216, [10]:2664, [12]:16600, [15]:624, [18]:8640, [30]:9648, [36]:27000, [45]:5184, [90]:12960
    // Generators of Hol[D90]
    // gen1 of order 2
    // [(7 14)(8 13)(9 12)(10 11)(21 28)(22 27)(23 26)(24 25)]
    // gen2 of order 12
    // [(2 3 5 4)(7 8 10 14 13 11)(9 12)(16 17 19 18)(21 22 24 28 27 25)(23 26)]
    // gen3 of order 90
    // [(1 15 2 16 3 17 4 18 5 19)(6 20 7 21 8 22 9 23 10 24 11 25 12 26 13 27 14 28)]
    // 
    // Aut[D90] = C45 x: (C12 x C2)
    // Hol[D90] = D90 x: Aut[D90]
    // 
    // |Hol[D90]| = |D90| * |Aut[D90]|             True
    // D90        is normal subgroup of Hol[D90]   True
    // Aut[D90]   is normal subgroup of Hol[D90]   False
    // D90        ∩   Aut[D90] = { [] }
    // 

    public static void Example7Dicyclic()
    {
        GlobalStopWatch.Restart();
        for (int m = 2; m <= 32; ++m)
        {
            var (name, gensDicm, gensAut, gensHol) = HolDicmPerm(m);
            var sn = gensDicm[0].Sn;
            var dicm = Group.Generate(name, sn, gensDicm);
            var autDicm = Group.Generate($"Aut[{name}]", sn, gensAut);
            var holDicm = Group.Generate($"Hol[{name}]", sn, gensHol);
            DisplayGroup.HeadOrdersGenerators(dicm);
            DisplayGroup.HeadOrdersGenerators(autDicm);
            DisplayGroup.HeadOrdersGenerators(holDicm);
            var q4mNormal = Group.IsNormalSubgroup(holDicm, dicm);
            var autQ4mNotNormal = Group.IsNormalSubgroup(holDicm, autDicm);
            var inter = dicm.Intersect(autDicm).ToXSet();
            var holOrd = dicm.Order * autDicm.Order;
            Console.WriteLine($"{($"|{holDicm}|"),-20} = {($"|{dicm}| * |{autDicm}|"),-40} {holDicm.Order == holOrd}");
            Console.WriteLine($"{dicm,-20} is normal subgroup of {holDicm,-20} {q4mNormal}");
            Console.WriteLine($"{autDicm,-20} is normal subgroup of {holDicm,-20} {autQ4mNotNormal}");
            Console.WriteLine($"{dicm,-20} ∩ {autDicm,20} = {inter}");
            Console.WriteLine();
            if (!q4mNormal || autQ4mNotNormal || inter.Count != 1 || holDicm.Order != holOrd)
                throw new();

            if (m <= 16)
            {
                var (dicmReg, autReg, holReg) = HolomorphPerm(FG.DicyclicPg(m));
                if (!dicmReg.ElementsOrdersList().SequenceEqual(dicm.ElementsOrdersList()) ||
                    !autReg.ElementsOrdersList().SequenceEqual(autDicm.ElementsOrdersList()) ||
                    !holReg.ElementsOrdersList().SequenceEqual(holDicm.ElementsOrdersList()))
                {
                    DisplayGroup.HeadOrdersGenerators(holReg);
                    throw new();
                }
            }
        }

        GlobalStopWatch.Show(); // Time:41.192s
    }
    // |Hol[C15 x: Q8]| = 115200
    // Type        NonAbelianGroup
    // BaseGroup   S24
    // 
    // Elements Orders : [1]:1, [2]:5679, [3]:8, [4]:24400, [5]:24, [6]:5112, [8]:5760, [10]:6696, [12]:25856, [15]:192, [20]:4800, [24]:11520, [30]:7488, [40]:3840, [60]:6144, [120]:7680
    // Generators of Hol[C15 x: Q8]
    // gen1 of order 2
    // [(18 20)(22 24)]
    // gen2 of order 2
    // [(2 3)(10 11)]
    // gen3 of order 4
    // [(5 6 8 7)(13 14 16 15)]
    // gen4 of order 4
    // [(2 3)(10 11)(17 18 19 20)]
    // gen5 of order 15
    // [(1 2 3)(4 5 6 7 8)]
    // gen6 of order 120
    // [(1 9 2 10 3 11)(4 12 5 13 6 14 7 15 8 16)(17 21 18 22 19 23 20 24)]
    // 
    // |Hol[C15 x: Q8]|     = |C15 x: Q8| * |Aut[C15 x: Q8]|           True
    // C15 x: Q8            is normal subgroup of Hol[C15 x: Q8]       True
    // Aut[C15 x: Q8]       is normal subgroup of Hol[C15 x: Q8]       False
    // C15 x: Q8            ∩       Aut[C15 x: Q8] = { [] }
}
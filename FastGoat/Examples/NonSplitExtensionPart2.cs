using System.Numerics;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Perms;
using FastGoat.UserGroup.Words;

namespace FastGoat.Examples;

public static class NonSplitExtensionPart2
{
    // Alejandro Adem, R. James Milgram
    // Cohomology of Finite Groups
    // Chapter I Group Extensions, Simple Algebras and Cohomology
    static bool CheckTwistedL<Tn, Tg>(ConcreteGroup<Tn> N, ConcreteGroup<Tg> G, MapGroups<Tg, Automorphism<Tn>> L,
        MapGroups<Ep2<Tg, Tg>, Tn> w)
        where Tg : struct, IElt<Tg>
        where Tn : struct, IElt<Tn>
    {
        // ω(s, e) = ω(e, s) = e.
        var e = G.Neutral();
        var wes = G.SelectMany(s => new[] { w[new(e, s)], w[new(s, e)] }).ToHashSet();
        if (!wes.SetEquals(new[] { N.Neutral() }))
            return false;

        foreach (var (r, s, t) in G.Grid3D(G, G))
        {
            var lti = L[t].Invert();
            var e0 = N.Op(lti[w[new(r, s)]], w[new(G.Op(r, s), t)]);
            var e1 = N.Op(w[new(s, t)], w[new(r, G.Op(s, t))]);
            if (!e0.Equals(e1))
                return false;
        }

        return true;
    }
        
    // Alejandro Adem, R. James Milgram
    // Cohomology of Finite Groups
    // Chapter I Group Extensions, Simple Algebras and Cohomology
    static HashSet<(MapGroups<Tg, Automorphism<Tn>>, MapGroups<Ep2<Tg, Tg>, Tn>)> TwistedL<Tn, Tg>(ConcreteGroup<Tn> N, ConcreteGroup<Tg> G)
        where Tg : struct, IElt<Tg> 
        where Tn : struct, IElt<Tn>
    {
        var og = G.Count();
        var nG = G.Neutral();
        var nN = N.Neutral();
        var GxG = Product.Generate(G, G);
        var autN = Group.AutomorphismGroup(N);
        var innN = Group.InnerAutomorphismGroup(N);
        var nAN = autN.Neutral();
        var arrG = G.Where(g0 => !g0.Equals(nG)).ToArray();

        var allTwisted = new HashSet<(MapGroups<Tg, Automorphism<Tn>>, MapGroups<Ep2<Tg, Tg>, Tn>)>();
        var we = GxG.Where(e => e.E1.Equals(nG) || e.E2.Equals(nG)).ToArray();
        var we_nN = we.Select(e => (e, nN)).ToArray();
        var rem = GxG.Except(we).ToArray();
        var nb = rem.Length;
    
        var max1 = BigInteger.Pow(autN.Count(), og - 1);
        var max2 = BigInteger.Pow(N.Count(), nb);
        if (max1 * max2 > 5000)
            throw new($"too much iterations, {max1}x{max2}={max1*max2} > 5000");

        var allL = autN.MultiLoop(og - 1).Select(l => arrG.Zip(l).Prepend((nG, nAN)).ToDictionary(e => e.Item1, e => e.Item2))
            .Where(L => innN.SuperSetOf(GxG.Select(e => autN.Op(L[e.E1].Invert(), autN.Op(L[e.E2].Invert(), L[G.Op(e.E1, e.E2)])))))
            .Select(m => new MapGroups<Tg, Automorphism<Tn>>(G, autN, m))
            .ToHashSet();
    
        foreach (var l in N.MultiLoop(nb))
        {
            var l0 = l.ToArray();
            var map = new MapGroups<Ep2<Tg, Tg>, Tn>(GxG, N, rem.Zip(l0).Concat(we_nN).ToDictionary(a => a.Item1, a => a.Item2));
            foreach (var L in allL)
            {
                if (CheckTwistedL(N, G, L, map))
                    allTwisted.Add((L, map));
            }
        }

        return allTwisted;
    }

    static IEnumerable<ExtensionGroup<Tn, Tg>> NonSplitExtensions<Tn, Tg>(ConcreteGroup<Tn> N, ConcreteGroup<Tg> G)
        where Tg : struct, IElt<Tg>
        where Tn : struct, IElt<Tn>
    {
        var twistedMaps = TwistedL(N, G);
        Console.WriteLine($"Twisted Actions : {twistedMaps.Count}");
        var count = 0;
        foreach (var (L, map) in twistedMaps)
        {
            var ext = new ExtensionGroup<Tn, Tg>(N, L, map, G);
            var isGroup = Group.IsGroup(ext);
            if (!isGroup)
                continue;

            Console.WriteLine($"{N} \u22b2 E{++count,-2} \u2192 {G}");
            yield return ext;
        }
    }

    public static void ExamplesOrder4()
    {
        var c2 = new Cn(2);
        foreach (var ext in NonSplitExtensions(c2, c2).ToArray())
        {
            DisplayGroup.HeadElements(ext);
            Console.WriteLine($"{ext} is isomorphic to {AbelianInvariantsFactors.Reduce(ext).Glue(" x ", "C{0}")}");

            Console.WriteLine("L : G -> Aut(N)");
            Console.WriteLine(ext.L);
            Console.WriteLine("f : G x G -> N");
            Console.WriteLine(ext.Map);
            NonSplitExtension.AllSplittingGroups(c2, ext, c2).First();
            Console.WriteLine();
        }
    }

    public static void ExamplesOrder8()
    {
        var (c2, c4, d8, q8) = (new Cn(2), new Cn(4), FG.Dihedral(4), FG.Quaternion(8));
        var sm8 = new dynamic[] { d8, q8 };
        var k = 1;
        foreach (var ext in NonSplitExtensions(c4, c2).ToHashSet(new IsomorphEquality<Ep2<ZnInt, ZnInt>>()))
        {
            ext.SetName($"({ext.Name})_{k++}");
            DisplayGroup.HeadElements(ext);
            if (ext.GroupType == GroupType.AbelianGroup)
            {
                Console.WriteLine($"{ext} is isomorphic to {AbelianInvariantsFactors.Reduce(ext).Glue(" x ", "C{0}")}");
            }
            else
            {
                var iso = sm8.First(e => ext.IsIsomorphicTo(e));
                Console.WriteLine($"{ext} is isomorphic to {iso}");
            }

            NonSplitExtension.AllSplittingGroups(c4, ext, c2).First();
            Console.WriteLine();
        }
    }

    public static void ExamplesOrder16()
    {
        var (q8, c2) = (FG.WordGroup("Q8", "a4, a2=b2, aba=b"), new Cn(2));
        var sm16 = new dynamic[]
        {
            FG.SemiDihedral(4), FG.Quaternion(16), Product.Generate(c2, q8),
            Group.SemiDirectProd(FG.Abelian(4, 2), c2)
        };
        var allExts = NonSplitExtensions(q8, c2).ToHashSet(new IsomorphEquality<Ep2<Word, ZnInt>>());
        var k = 1;
        foreach (var ext in allExts)
        {
            var q8ext = Group.IsomorphicsSubgroupsAll(ext, q8);
            if (q8ext.Count != 0)
            {
                ext.SetName($"({ext.Name})_{k++}");
                Console.WriteLine("#############################################");
                DisplayGroup.HeadElements(ext);
                var iso = sm16.FirstOrDefault(e => e.IsIsomorphicTo(ext), ext);
                DisplayGroup.AreIsomorphics(ext, iso);

                var e0 = q8ext.First();
                var conjs = Group.SubGroupsConjugates(ext, e0);
                if (conjs.Count == 1)
                    Console.WriteLine($"{q8} is normal subgroup of {iso}");
            
                Console.WriteLine();
                NonSplitExtension.AllSplittingGroups(q8, ext, c2).First();
                Console.WriteLine();
            }
        }
    }

    public static void ExampleOrder32()
    {
        var s32 = new Sn(32);
        var g1 = s32[(1, 2, 3, 4), (5, 6, 7, 8), (9, 10, 11, 12), (13, 14, 15, 16), (17, 18, 19, 20), (21, 22, 23, 24),
            (25, 26, 27, 28), (29, 30, 31, 32)];
        var g2 = s32[(1, 17, 15, 6), (2, 18, 16, 7), (3, 19, 13, 8), (4, 20, 14, 5), (9, 29, 24, 28), (10, 30, 21, 25),
            (11, 31, 22, 26), (12, 32, 23, 27)];
        var g3 = s32[(1, 26, 15, 31), (2, 32, 16, 27), (3, 28, 13, 29), (4, 30, 14, 25), (5, 12, 20, 23),
            (6, 24, 17, 9), (7, 10, 18, 21), (8, 22, 19, 11)];

        var sm3232 = Group.Generate("[(C4 x C4) . C2]pg", s32, g1, g2, g3);
        var ext = NonSplitExtensions(FG.Abelian(4, 4), new Cn(2)).First(g => g.IsIsomorphicTo(sm3232));
        DisplayGroup.HeadElements(ext);
        DisplayGroup.AreIsomorphics(ext, sm3232);
        NonSplitExtension.AllSplittingGroups(FG.Abelian(4, 4), ext, new Cn(2)).First();
        Console.WriteLine();
    }
    /*
       C4 x C4 ⊲ E48 → C2 (step 48)
       |(C4 x C4) . C2| = 32
       Type        NonAbelianGroup
       BaseGroup   (C4 x C4) . C2
       
       Elements
       ( 1)[1] = ((0, 0), 0)
       ( 2)[2] = ((0, 2), 0)
       ( 3)[2] = ((2, 0), 0)
       ( 4)[2] = ((2, 2), 0)
       ( 5)[4] = ((0, 0), 1)
       ( 6)[4] = ((0, 1), 0)
       ( 7)[4] = ((0, 1), 1)
       ( 8)[4] = ((0, 2), 1)
       ( 9)[4] = ((0, 3), 0)
       (10)[4] = ((0, 3), 1)
       (11)[4] = ((1, 0), 0)
       (12)[4] = ((1, 0), 1)
       (13)[4] = ((1, 1), 0)
       (14)[4] = ((1, 1), 1)
       (15)[4] = ((1, 2), 0)
       (16)[4] = ((1, 2), 1)
       (17)[4] = ((1, 3), 0)
       (18)[4] = ((1, 3), 1)
       (19)[4] = ((2, 0), 1)
       (20)[4] = ((2, 1), 0)
       (21)[4] = ((2, 1), 1)
       (22)[4] = ((2, 2), 1)
       (23)[4] = ((2, 3), 0)
       (24)[4] = ((2, 3), 1)
       (25)[4] = ((3, 0), 0)
       (26)[4] = ((3, 0), 1)
       (27)[4] = ((3, 1), 0)
       (28)[4] = ((3, 1), 1)
       (29)[4] = ((3, 2), 0)
       (30)[4] = ((3, 2), 1)
       (31)[4] = ((3, 3), 0)
       (32)[4] = ((3, 3), 1)
       
       (C4 x C4) . C2 IsIsomorphicTo [(C4 x C4) . C2]pg : True
       1 -----> G:C4 x C4 -----> GH:(C4 x C4) . C2 -----> H:C2 -----> 1
       G  --i--> GH
       GH --p--> H
       H  --s--> GH
       Searching 1 -----> G --i--> GH --p--> H -----> 1 with Im(i)=Ker(p)
       Morphism 'i' from G to GH
       [(0, 0)->[((0, 0), 0)]; (0, 1)->[((0, 1), 0)]; (0, 2)->[((0, 2), 0)]; (0, 3)->[((0, 3), 0)]; (1, 0)->[((1, 0), 0)]; (1, 1)->[((1, 1), 0)]; (1, 2)->[((1, 2), 0)]; (1, 3)->[((1, 3), 0)]; (2, 0)->[((2, 0), 0)]; (2, 1)->[((2, 1), 0)]; (2, 2)->[((2, 2), 0)]; (2, 3)->[((2, 3), 0)]; (3, 0)->[((3, 0), 0)]; (3, 1)->[((3, 1), 0)]; (3, 2)->[((3, 2), 0)]; (3, 3)->[((3, 3), 0)]]
       Morphism 'p' from GH to H
       [((0, 0), 0)->[0]; ((0, 0), 1)->[1]; ((0, 1), 0)->[0]; ((0, 1), 1)->[1]; ((0, 2), 0)->[0]; ((0, 2), 1)->[1]; ((0, 3), 0)->[0]; ((0, 3), 1)->[1]; ((1, 0), 0)->[0]; ((1, 0), 1)->[1]; ((1, 1), 0)->[0]; ((1, 1), 1)->[1]; ((1, 2), 0)->[0]; ((1, 2), 1)->[1]; ((1, 3), 0)->[0]; ((1, 3), 1)->[1]; ((2, 0), 0)->[0]; ((2, 0), 1)->[1]; ((2, 1), 0)->[0]; ((2, 1), 1)->[1]; ((2, 2), 0)->[0]; ((2, 2), 1)->[1]; ((2, 3), 0)->[0]; ((2, 3), 1)->[1]; ((3, 0), 0)->[0]; ((3, 0), 1)->[1]; ((3, 1), 0)->[0]; ((3, 1), 1)->[1]; ((3, 2), 0)->[0]; ((3, 2), 1)->[1]; ((3, 3), 0)->[0]; ((3, 3), 1)->[1]]
       Isomorphism 's' from H to GH
       Not Found
     */
}
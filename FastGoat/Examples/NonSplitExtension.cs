using FastGoat.Structures.CartesianProduct;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Perms;

namespace FastGoat.Examples;

public static class NonSplitExtension
{
    static void SearchQuartenionInSymmetric(int n)
    {
        var sn = Group.Create(new Sn(n));
        var allC4 = sn.Where(e => sn.ElementsOrders[e] == 4).Ascending().ToArray();

        var query = from a in allC4
            from b in allC4
            select (a, b);

        var tuple = query.FirstOrDefault(e =>
                (e.a ^ 2) == (e.b ^ 2) && ((e.b * e.a * (e.b ^ -1)) == (e.a ^ -1))
            , (a: sn.Neutral(), b: sn.Neutral()));

        if (tuple == (sn.Neutral(), sn.Neutral()))
        {
            Console.WriteLine("Quartenion Not found in {0}", sn);
        }
        else
        {
            var q8 = Group.Generate("Q8", sn, tuple.a, tuple.b);
            DisplayGroup.HeadElementsTable(q8);
        }
    }

    static void SplittingGroups<T1, T2, T3>(ConcreteGroup<T1> G, ConcreteGroup<T2> GH, ConcreteGroup<T3> H)
        where T1 : struct, IElt<T1>
        where T2 : struct, IElt<T2>
        where T3 : struct, IElt<T3>
    {
        Console.WriteLine("1 -----> G:{0} -----> GH:{1} -----> H:{2} -----> 1", G, GH, H);
        Console.WriteLine("G  --i--> GH");
        Console.WriteLine("GH --p--> H");
        Console.WriteLine("H  --s--> GH");
        var homI = Group.AllHomomorphisms(G, GH);
        var homP = Group.AllHomomorphisms(GH, H);
        var isoS = Group.AllIsomorphisms(H, GH);

        // Im(i) = Ker(p)
        var allImI = homI.Select(hi => (i: hi, im: hi.Image().ToHashSet())).Where(e => e.im.Count > 1).ToArray();
        var allKerP = homP.Select(hp => (p: hp, ker: hp.Kernel().ToHashSet())).Where(e => e.ker.Count > 1)
            .ToArray();
        var ImIeqKerP = allImI
            .SelectMany(e1 => allKerP.Where(e2 => e1.im.SetEquals(e2.ker)).Select(e2 => (e1.i, e2.p)))
            .ToArray();

        Console.WriteLine("Searching 1 -----> G --i--> GH --p--> H -----> 1 with Im(i)=Ker(p)");
        foreach (var (i, p) in ImIeqKerP)
        {
            Console.WriteLine("Morphism 'i' from G to GH");
            Console.WriteLine("    [{0}]", i);
            Console.WriteLine("Morphism 'p' from GH to H");
            Console.WriteLine("    [{0}]", p);

            Console.WriteLine("Isomorphism 's' from H to GH");
            var allIsoS = isoS
                .Where(s0 => s0.Count == H.Count() && H.All(h => p.Domain.Contains(s0[h]) && p[s0[h]].Equals(h)))
                .ToArray();
            foreach (var s in allIsoS)
            {
                Console.WriteLine("    [{0}]", s);
            }

            if (allIsoS.Length == 0)
                Console.WriteLine("    Not Found");

            Console.WriteLine();
        }

        Console.WriteLine("End");
        Console.WriteLine();
    }

    public static void SearchingQuartenion()
    {
        SearchQuartenionInSymmetric(4);
        SearchQuartenionInSymmetric(5);
        SearchQuartenionInSymmetric(6);
        SearchQuartenionInSymmetric(7);
        SearchQuartenionInSymmetric(8);
    }

    public static void SplittingSemiDirectProduct()
    {
        var s7 = new Symm(7);
        var a = s7[(1, 2, 3, 4, 5, 6, 7)];
        var b = s7[(2, 3, 5), (4, 7, 6)];
        var sdp = Group.Generate("C7 x: C3", s7, a, b);
        // var c7 = Group.Generate("C7", sdp, a);
        // var c3 = Group.Generate("C3", sdp, b);
        var c7 = new Cn(7);
        var c3 = new Cn(3);

        SplittingGroups(c3, sdp, c7); // c3 isnot a normal subgroup
        SplittingGroups(c7, sdp, c3);
    }

    public static void SplittingS3andC6()
    {
        var c2 = new Cn(2);
        var c3 = new Cn(3);
        var c6 = new Cn(6);
        var s3 = new Symm(3);

        SplittingGroups(c2, c6, c3);
        SplittingGroups(c3, c6, c2);

        SplittingGroups(c2, s3, c3); // c2 isnot a normal subgroup
        SplittingGroups(c3, s3, c2);
    }

    public static void SplittinDirectProduct()
    {
        var c12 = new Cn(12);
        var c4 = new Cn(4);
        var c3 = new Cn(3);

        SplittingGroups(c3, c12, c4);
        SplittingGroups(c4, c12, c3);
    }

    public static void SplittingDihedral()
    {
        var s4 = new Symm(4);
        var d8 = Group.Generate("D8", s4, s4[(1, 4), (2, 3)], s4[(1, 2, 3, 4)]);
        var c4 = Group.Generate("C4", d8, s4[(1, 2, 3, 4)]);
        var c2 = Group.Generate("C2", d8, s4[(1, 3), (2, 4)]);
        var v = Group.Generate("V", d8, s4[(1, 3), (2, 4)], s4[(1, 2), (3, 4)]);

        SplittingGroups(c2, d8, v);
        SplittingGroups(c2, d8, c4);
        SplittingGroups(v, d8, c2);
        SplittingGroups(c4, d8, c2);
    }

    public static void SplittingQuartenion()
    {
        var s8 = new Sn(8);
        var a = s8[(1, 2, 3, 4), (5, 6, 7, 8)];
        var b = s8[(1, 4, 3, 2), (5, 8, 7, 6)];
        var q8 = Group.Generate("Q8", s8, a, b);

        var c2 = Group.Generate("C2", s8, s8[(1, 2)]);
        var c4 = Group.Generate("C4", s8, s8[(1, 2, 3, 4)]);
        var v = Group.Generate("V", s8, s8[(1, 3), (2, 4)], s8[(1, 4), (2, 3)]);

        SplittingGroups(c2, q8, c4);
        SplittingGroups(c4, q8, c2);
        SplittingGroups(v, q8, c2);
        SplittingGroups(c2, q8, v);
    }

    public static void SplittingSmallGroup_32_32()
    {
        // https://people.maths.bris.ac.uk/~matyd/GroupNames/1/C4%5E2.C2.html
        // Generators and relations for C42.C2
        // G = < a,b,c | a4=b4=1, c2=b2, ab=ba, cac-1=ab2, cbc-1=a2b >
        // G:=Group(    (1,2,3,4)(5,6,7,8)(9,10,11,12)(13,14,15,16)(17,18,19,20)(21,22,23,24)(25,26,27,28)(29,30,31,32),
        //              (1,17,15,6)(2,18,16,7)(3,19,13,8)(4,20,14,5)(9,29,24,28)(10,30,21,25)(11,31,22,26)(12,32,23,27),
        //              (1,26,15,31)(2,32,16,27)(3,28,13,29)(4,30,14,25)(5,12,20,23)(6,24,17,9)(7,10,18,21)(8,22,19,11) );

        var s32 = new Sn(32);
        var g1 = s32[(1, 2, 3, 4), (5, 6, 7, 8), (9, 10, 11, 12), (13, 14, 15, 16), (17, 18, 19, 20), (21, 22, 23, 24),
            (25, 26, 27, 28), (29, 30, 31, 32)];
        var g2 = s32[(1, 17, 15, 6), (2, 18, 16, 7), (3, 19, 13, 8), (4, 20, 14, 5), (9, 29, 24, 28), (10, 30, 21, 25),
            (11, 31, 22, 26), (12, 32, 23, 27)];
        var g3 = s32[(1, 26, 15, 31), (2, 32, 16, 27), (3, 28, 13, 29), (4, 30, 14, 25), (5, 12, 20, 23),
            (6, 24, 17, 9), (7, 10, 18, 21), (8, 22, 19, 11)];

        var sm3232 = Group.Generate("(C4 x C4) . C2", s32, g1, g2, g3);
        DisplayGroup.HeadElements(sm3232);

        var c2 = new Cn(2);
        var c4 = new Cn(4);
        var c4c4 = Product.Generate(c4, c4);
        SplittingGroups(c4c4, sm3232, c2); // Founding 96 non unique extensions but they never split
    }
}
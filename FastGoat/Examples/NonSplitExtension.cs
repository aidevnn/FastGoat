using FastGoat.Gp;
using FastGoat.UserGroup;

namespace FastGoat.Examples;

public static class NonSplitExtension
{
    static void SearchQuartenionInSymmetric(int n)
    {
        var sn = new Symm(n);
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
        var homS = Group.AllHomomorphisms(H, GH);

        // Im(i)= Ker(p)
        var allImI = homI.Where(hi => Group.Image(hi).Count > 1).Select(hi => (i: hi, im: Group.Image(hi))).ToArray();
        var allKerP = homP.Where(hp => Group.Kernel(H.Neutral(), hp).Count > 1)
            .Select(hp => (p: hp, ker: Group.Kernel(H.Neutral(), hp))).ToArray();
        var sols = allImI
            .SelectMany(e1 => allKerP.Where(e2 => e1.im.SetEquals(e2.ker)).Select(e2 => (e1.i, e2.p)))
            .ToArray();

        Console.WriteLine("Searching 1 -----> G --i--> GH --p--> H -----> 1 with Im(i)=Ker(p)");
        foreach (var (i, p) in sols)
        {
            Console.WriteLine("Morphism 'i' from G to GH");
            Console.WriteLine("    [{0}]", i.GlueMap());
            Console.WriteLine("Morphism 'p' from GH to H");
            Console.WriteLine("    [{0}]", p.GlueMap());


            Console.WriteLine("Isomorphism 's' from H to GH");
            var allS = homS.Where(s0 => s0.Count == H.Count() && H.All(h => p.ContainsKey(s0[h]) && p[s0[h]].Equals(h)))
                .ToArray();
            foreach (var s in allS)
            {
                Console.WriteLine("    [{0}]", s.GlueMap());
            }

            if (allS.Length == 0)
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
        var s8 = new Symm(8);
        var a = s8[(1, 2, 3, 4), (5, 6, 7, 8)];
        var b = s8[(1, 4, 3, 2), (5, 8, 7, 6)];
        var q8 = Group.Generate("Q8", s8, a, b);

        var c2 = Group.Generate("C2", s8[(1, 2)]);
        var c4 = Group.Generate("C4", s8[(1, 2, 3, 4)]);
        var v = Group.Generate("V", s8[(1, 3), (2, 4)], s8[(1, 4), (2, 3)]);

        SplittingGroups(c2, q8, c4);
        SplittingGroups(c4, q8, c2);
        SplittingGroups(v, q8, c2);
        SplittingGroups(c2, q8, v);
    }

    public static void SplittingSmallGroup_32_32()
    {
        // TODO
    }
}
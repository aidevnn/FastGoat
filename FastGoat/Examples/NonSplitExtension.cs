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
        var homI = Group.AllHomomorphisms(G, GH);
        var homP = Group.AllHomomorphisms(GH, H);

        // Im(i)= Ker(p)
        var allImI = homI.Select(hi => (hi, hi.Values.Distinct().Ascending().ToArray())).ToArray();
        var allKerP = homP.Select(hp =>
            (hp, hp.Where(p1 => p1.Value.Equals(H.Neutral())).Select(p2 => p2.Key).Ascending().ToArray())
        ).ToArray();

        var homI0 = allImI.Where(e1 => allKerP.Any(e2 => e1.Item2.SequenceEqual(e2.Item2))).Select(e1 => e1.hi)
            .Where(e1 => e1.Values.Distinct().Count() != 1)
            .ToArray();

        var homP0 = allKerP.Where(e1 => allImI.Any(e2 => e1.Item2.SequenceEqual(e2.Item2))).Select(e1 => e1.hp)
            .Where(e1 => e1.Values.Distinct().Count() != 1)
            .ToArray();

        Console.WriteLine("Searching 1 -----> G --i--> GH --p--> H -----> 1 with Im(i)=Ker(p)");
        Console.WriteLine("Searching a morphism 'i' from G to GH");
        foreach (var hi in homI0)
        {
            Console.WriteLine(hi.GlueMap());
        }

        Console.WriteLine();
        Console.WriteLine("Searching a morphism 'p' from GH to H");
        foreach (var hi in homP0)
        {
            Console.WriteLine(hi.GlueMap());
        }

        Console.WriteLine();
        Console.WriteLine("Searching if extension can split H --s--> GH with s(H)=H");
        var allS = Group.AllHomomorphisms(H, GH);
        foreach (var s in allS.Where(s0 =>
                     s0.Values.Distinct().Count() != 1 && H.Select(h => s0[h]).Distinct().Count() == H.Count()))
        {
            Console.WriteLine("[{0}]", s.GlueMap());
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
        var c4 = Group.Generate("C4", s4, s4[(1, 2, 3, 4)]);
        var c2 = Group.Generate("C2", s4, s4[(1, 3), (2, 4)]);
        var v = Group.Generate("V", s4, s4[(1, 3), (2, 4)], s4[(1, 2), (3, 4)]);
        var d8 = Group.Generate("D8", s4, s4[(1, 3), (2, 4)], s4[(1, 2, 3, 4)]);

        SplittingGroups(v, d8, c2);
        SplittingGroups(c2, d8, v);
        SplittingGroups(c2, d8, c4);
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
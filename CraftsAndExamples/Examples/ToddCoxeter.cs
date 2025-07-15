using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Perms;
using FastGoat.UserGroup.Words;
using FastGoat.UserGroup.Words.Tools;

namespace CraftsAndExamples.Examples;

public static class ToddCoxeter
{
    static ToddCoxeter()
    {
        Logger.Level = LogLevel.Level2;
    }
    public static void CyclicGroup()
    {
        Graph.RunToddCoxeterAlgo("a5");
    }

    public static void KleinGroup()
    {
        Graph.RunToddCoxeterAlgo("a2, b2, ab = ba");
    }

    public static void Symm3Group()
    {
        Graph.RunToddCoxeterAlgo("a2, b3, abab");
    }

    public static void MoreExamples()
    {
        // Algebre Tome 1, Daniel Guin – Thomas Hausberger

        Graph.RunToddCoxeterAlgo("a", "a4, b3, abab"); // p101 step by step algorithm with subgroup H=<a>
        Graph.RunToddCoxeterAlgo("a", "a3, b3, abab"); // p106 step by step algorithm with subgroup H=<a>
        Graph.RunToddCoxeterAlgo("a", "a3, b3, aba2b"); // p107 step by step algorithm with subgroup H=<a>

        Graph.RunToddCoxeterAlgo("b", "a7, b3, a2=bab-1"); // C7 : C3 step by step algorithm with subgroup H=<b>
        Graph.RunToddCoxeterAlgo("a7, b3, a2=bab-1"); // C7 : C3 step by step algorithm with subgroup H=<Id>

        Graph.RunToddCoxeterAlgo("a2, b4, ab=ba"); // Dihedral 8 step by step algorithm
        Graph.Run("a4, a2=b2, b-1aba").DisplayTableOps(); // Quartenion Table
    }

    public static void SmallGroup32_32()
    {
        Logger.Level = LogLevel.Off;
        Graph.Run("a4, b4, c2=b2, ab=ba, cac-1=ab2, cbc-1=a2b").DisplayTableOps(); // (C4 x C4) . C2 Table
        Graph.Run("a, b", "a4, b4, c2=b2, ab=ba, cac-1=ab2, cbc-1=a2b").DisplayTableOps(); // C2 Table with subgroup H = C4 x C4=<a,b>
    }

    public static void Quaternion()
    {
        Logger.Level = LogLevel.Off;
        var gname = "Q8";
        var relators = "a4, a2=b2, b-1aba";
        var wg = new WordGroup(gname, relators);
        DisplayGroup.HeadElementsTable(wg);

        var s8 = new Sn(8);
        var a = s8[(1, 2, 3, 4), (5, 6, 7, 8)];
        var b = s8[(1, 5, 3, 7), (2, 8, 4, 6)];
        var q8 = Group.Generate("Q8", s8, a, b);
        Console.WriteLine("({0}) IsIsomorphicTo ({1}) : {2}", wg, q8, wg.IsIsomorphicTo(q8));

        GlobalStopWatch.Restart();
        var wg0 = new WordGroup(gname, relators);
        GlobalStopWatch.Show($"{wg0.Name}");

        GlobalStopWatch.Restart();
        var g0 = Group.Generate("Q8", s8, a, b);
        GlobalStopWatch.Show($"{g0.Name}");
    }

    public static void Frobenius20()
    {
        Logger.Level = LogLevel.Off;
        var gname = "F20";
        var relators = "a5, b4, b-1ab=a2";
        var wg = new WordGroup(gname, relators);
        DisplayGroup.HeadElementsTable(wg);

        var gr = Group.SemiDirectProd(new Cn(5), new Cn(4));
        Console.WriteLine("({0}) IsIsomorphicTo ({1}) : {2}", wg, gr, wg.IsIsomorphicTo(gr));

        GlobalStopWatch.Restart();
        var wg0 = new WordGroup(gname, relators);
        GlobalStopWatch.Show($"{wg0.Name}");

        GlobalStopWatch.Restart();
        var g0 = Group.SemiDirectProd(new Cn(5), new Cn(4));
        GlobalStopWatch.Show($"{g0.Name}");
    }

    public static void DiCyclic3()
    {
        Logger.Level = LogLevel.Off;
        var wg = new WordGroup("Dic3", "a6, b2=a3, bab-1=a-1");
        var g = Group.SemiDirectProd(new Cn(3), new Cn(4));
        DisplayGroup.HeadElements(wg);
        DisplayGroup.HeadElementsSdp(g);
        Console.WriteLine("({0}) IsIsomorphicTo ({1}) : {2}", wg, g, wg.IsIsomorphicTo(g));
    }

    public static void DiCyclic12()
    {
        Logger.Level = LogLevel.Off;
        var gname = "wgDic12";
        var relators = "a12, b2=a6, bab-1=a-1";
        var wg = new WordGroup(gname, relators);
        DisplayGroup.HeadElements(wg);

        var s8 = new Sn(8);
        var a0 = s8[(1, 2, 3, 4), (5, 6, 7, 8)];
        var b0 = s8[(1, 5, 3, 7), (2, 8, 4, 6)];
        var q8 = Group.Generate("Q8", s8, a0, b0);
        var sp = Group.SemiDirectProd(new Cn(3), q8);
        DisplayGroup.HeadElements(sp);
        Console.WriteLine("({0}) IsIsomorphicTo ({1}) : {2}", wg, sp, wg.IsIsomorphicTo(sp));
        
        var s24 = new Sn(24);
        var a = s24[(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), (13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24)];
        var b = s24[(1, 17, 7, 23), (2, 16, 8, 22), (3, 15, 9, 21), (4, 14, 10, 20), (5, 13, 11, 19), (6, 24, 12, 18)];
        var dic12 = Group.Generate("pgDic12", s24, a, b);
        DisplayGroup.HeadElements(dic12);
        Console.WriteLine("({0}) IsIsomorphicTo ({1}) : {2}", wg, dic12, wg.IsIsomorphicTo(dic12));
    }

    public static void Quaternion32()
    {
        Logger.Level = LogLevel.Off;
        var gname = "wgQ32";
        var relators = "a16, b2=a8, bab-1=a-1";
        var wg = new WordGroup(gname, relators);
        DisplayGroup.HeadElements(wg);

        var s32 = new Sn(32);
        var a = s32[(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16),
            (17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32)];
        var b = s32[(1, 30, 9, 22), (2, 29, 10, 21), (3, 28, 11, 20), (4, 27, 12, 19), (5, 26, 13, 18), (6, 25, 14, 17),
            (7, 24, 15, 32), (8, 23, 16, 31)];
        var q32 = Group.Generate("pgQ32", s32, a, b);
        DisplayGroup.HeadElements(q32);
        Console.WriteLine("({0}) IsIsomorphicTo ({1}) : {2}", wg, q32, wg.IsIsomorphicTo(q32));
    }

    public static void Coincidences1()
    {
        // Algebre Tome 1, Daniel Guin – Thomas Hausberger
        // p107 step by step algorithm
        Graph.RunToddCoxeterAlgo("a", "a3, b3, aba2b");
        Graph.RunToddCoxeterAlgo("a3, b3, aba2b");
    }

    public static void Coincidences2()
    {
        // Ken Brown paper toddcox.pdf
        Graph.RunToddCoxeterAlgo("aba-1 = b2, bab-1 = a2");
    }

    public static void Group_576_8282()
    {
        // Ken Brown paper toddcox.pdf
        Graph.RunToddCoxeterAlgo("a,b", "a3,b2,c2,abababab,acac,bcbcbc");

        Logger.Level = LogLevel.Off;
        GlobalStopWatch.Restart();
        GlobalStopWatch.AddLap();
        var g = FG.WordGroup("a3,b2,c2,abababab,acac,bcbcbc");
        DisplayGroup.HeadOrders(g);
        GlobalStopWatch.Show("WG");
        GlobalStopWatch.AddLap();
        DisplayGroup.HeadOrdersNames(g);
        GlobalStopWatch.Show("Names");
        GlobalStopWatch.Show("End");
    }
    /*
       |SL(2,3) x: S4| = 576
       Type        NonAbelianGroup
       BaseGroup   WG[a,b,c]

       Elements Orders : [1]:1, [2]:91, [3]:80, [4]:84, [6]:80, [8]:144, [12]:96
       SubGroupsInfos { AllSubGr = 1731, AllConjsCl = 127, AllNorms = 13 }
       Group names
           SL(2,3) x: S4
           (SL(2,3) x: A4) x: C2
           (SL(2,3) x: (C2 x C2)) x: S3
           (Q8 x: A4) x: S3
           (D8 x: (C2 x C2)) x: ((C3 x C3) x: C2)
           C2 . (A4 x: S4)
           Q8 . (C3 x: S4)

       #  End Time:39.947s
     */

    public static void Symm6a()
    {
        Logger.Level = LogLevel.Level1;
        Graph.RunToddCoxeterAlgo("a", "a2, b6, ababababab, ab2ab-2ab2ab-2, abab-1abab-1abab-1");
    }

    public static void Symm6b()
    {
        Logger.Level = LogLevel.Level1;
        Graph.RunToddCoxeterAlgo("b", "a2, b6, ababababab, ab2ab-2ab2ab-2, abab-1abab-1abab-1");
    }

    public static void Symm6()
    {
        Logger.Level = LogLevel.Level1;
        Graph.RunToddCoxeterAlgo("a2, b6, ababababab, ab2ab-2ab2ab-2, abab-1abab-1abab-1");
    }

    public static void Symm6Orders()
    {
        Logger.Level = LogLevel.Off;
        GlobalStopWatch.Restart();
        var g = FG.WordGroup("a2, b6, ababababab, ab2ab-2ab2ab-2, abab-1abab-1abab-1");
        DisplayGroup.HeadOrders(g);
        GlobalStopWatch.Show();
    }

    public static void PermGroup7()
    {
        Logger.Level = LogLevel.Level1;
        Graph.DefiningRelatorsOfGroup(FG.Alternate(7));
        Graph.RunToddCoxeterAlgo("a3, b3, c3, d3, e3, abab, acac, adad, aeae, bcbc, bdbd, bebe, cdcd, cece, dede");
        // Step:2521 NbClasses:2520
        // Time:20.648s

        Graph.DefiningRelatorsOfGroup(FG.Symmetric(7));
        Graph.RunToddCoxeterAlgo("a7, b2, abababababab, baba-1baba-1baba-1, a2ba-2ba2ba-2b"); // Symm7
        // Step:5787 NbClasses:5040
        // Time:44.875s
    }
    
    public static void L2p()
    {
        Logger.Level = LogLevel.Level1;
        Graph.DefiningRelatorsOfGroup(FG.L2p(11));
        /*
           |L2(11)| = 660
        
           #  Time:518ms
           All Relators
               a3
               b2
               ababababababababababab
               baba-1baba-1baba-1baba-1baba-1
               abababa-1ba-1ba-1babababa-1ba-1ba-1b
         */
        
        Graph.RunToddCoxeterAlgo("a3,b2,ababababababababababab,baba-1baba-1baba-1baba-1baba-1,abababa-1ba-1ba-1babababa-1ba-1ba-1b");
        // Step:683 NbClasses:660
        // #  Time:2.355s
        
        Graph.DefiningRelatorsOfGroup(FG.L2p(13));
        /* |L2(13)| = 1092
        
           #  Time:210ms
           All Relators
               a3
               b2
               ababababababababababababab
               abababa-1baba-1babababa-1baba-1b
               ababababa-1bababababa-1ba-1bababa-1ba-1b
         */
        
        Graph.RunToddCoxeterAlgo(
            "a3,b2,ababababababababababababab,abababa-1baba-1babababa-1baba-1b,ababababa-1bababababa-1ba-1bababa-1ba-1b");
        // Step:1099 NbClasses:1092
        // #  Time:3.939s
        
        Graph.DefiningRelatorsOfGroup(FG.L2p(17));
        /* |L2(17)| = 2448
        
           #  Time:547ms
           All Relators
               a3
               b2
               ababababa-1ba-1bababa-1baba-1bababa-1ba-1b
               abababababa-1ba-1ba-1babababababa-1ba-1ba-1b
         */
        
        Graph.RunToddCoxeterAlgo("a3,b2,ababababa-1ba-1bababa-1baba-1bababa-1ba-1b,abababababa-1ba-1ba-1babababababa-1ba-1ba-1b");
        // Step:2499 NbClasses:2448
        // #  Time:14.651s
    }
    
    public static void u33()
    {
        Logger.Level = LogLevel.Level1;
        var s28 = new Sn(28);
        var a2 = s28[(1, 5, 7, 3, 12, 24, 11), (2, 23, 4, 27, 13, 14, 26), (6, 20, 18, 8, 25, 21, 28), (9, 10, 17, 15, 22, 16, 19)];
        var b2 = s28[(3, 4), (5, 17, 7, 16, 8, 20, 6, 13), (9, 19, 11, 14, 12, 18, 10, 15), (21, 23, 26, 28, 24, 22, 27, 25)];
        var u3_3pg = Group.Generate("U3(3)pg", s28, a2, b2);
        Graph.DefiningRelatorsOfGroup(u3_3pg);
        /* |U3(3)pg| = 6048
        
           #  Time:2.945s
           All Relators
               a8
               b7
               ba-1ba-1ba-1
               a4ba3b-1a-1b2
               a2b-1ab2abab-1a-1b
         */
        
        Graph.RunToddCoxeterAlgo("a8,b7,ba-1ba-1ba-1,a4ba3b-1a-1b2,a2b-1ab2abab-1a-1b");
        // Step:8541 NbClasses:6048
        // #  Time:1m46s
    }
}
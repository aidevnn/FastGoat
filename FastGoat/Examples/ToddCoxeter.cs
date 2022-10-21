using FastGoat.ToddCoxeter;
using FastGoat.UserGroup;

namespace FastGoat.Examples;

public static class ToddCoxeter
{
    public static void CyclicGroup()
    {
        ToddCoxeterAlgo.Run("a5", details: true);
    }
    
    public static void KleinGroup()
    {
        ToddCoxeterAlgo.Run("a2, b2, ab = ba", details: true);
    }
    
    public static void Symm3Group()
    {
        ToddCoxeterAlgo.Run("a2, b3, abab", details: true);
    }

    public static void MoreExamples()
    {
        // Algebre Tome 1, Daniel Guin – Thomas Hausberger
        
        ToddCoxeterAlgo.Run("a", "a4, b3, abab", details: true); // p101 step by step algorithm with subgroup H=<a>
        ToddCoxeterAlgo.Run("a", "a3, b3, abab", details: true); // p106 step by step algorithm with subgroup H=<a>
        ToddCoxeterAlgo.Run("a", "a3, b3, aba2b", details: true); // p107 step by step algorithm with subgroup H=<a>
    
        ToddCoxeterAlgo.Run("b", "a7, b3, a2=bab-1", details: true); // C7 : C3 step by step algorithm with subgroup H=<b>
        ToddCoxeterAlgo.Run("a7, b3, a2=bab-1", details: true); // C7 : C3 step by step algorithm with subgroup H=<Id>
    
        ToddCoxeterAlgo.Run("a2, b4, ab=ba", details: true); // Dihedral 8 step by step algorithm
        ToddCoxeterAlgo.Run("a4, a2=b2, b-1aba").DisplayOps(); // Quartenion Table
    }

    public static void SmallGroup32_32()
    {
        ToddCoxeterAlgo.Run("a4, b4, c2=b2, ab=ba, cac-1=ab2, cbc-1=a2b").DisplayOps(); // (C4 x C4) . C2 Table
        ToddCoxeterAlgo.Run("a, b", "a4, b4, c2=b2, ab=ba, cac-1=ab2, cbc-1=a2b").DisplayOps(); // C2 Table with subgroup H = C4 x C4=<a,b>
    }

    public static void Quaternion()
    {
        var gname = "Q8";
        var relators = "a4, a2=b2, b-1aba";
        var wg = Group.Words(gname, relators);
        DisplayGroup.HeadElementsTable(wg);
        
        var s8 = new Symm(8);
        var a = s8[(1, 2, 3, 4), (5, 6, 7, 8)];
        var b = s8[(1, 5, 3, 7), (2, 8, 4, 6)];
        var q8 = Group.Generate("Q8", s8, a, b);
        Console.WriteLine("({0}) IsIsomorphicTo ({1}) : {2}", wg, q8, wg.IsIsomorphicTo(q8));
    
        GlobalStopWatch.Restart();
        var wg0 = Group.Words(gname, relators);
        GlobalStopWatch.Show($"{wg0.Name}");
    
        GlobalStopWatch.Restart();
        var g0 = Group.Generate("Q8", s8, a, b);
        GlobalStopWatch.Show($"{g0.Name}");
    }

    public static void Frobenius20()
    {
        var gname = "F20";
        var relators = "a5, b4, abababab, a2ba-1b-1";
        var wg = Group.Words(gname, relators);
        DisplayGroup.HeadElementsTable(wg);
        
        var gr = Group.SemiDirectProd(new Cn(5), new Cn(4));
        Console.WriteLine("({0}) IsIsomorphicTo ({1}) : {2}", wg, gr, wg.IsIsomorphicTo(gr));
    
        GlobalStopWatch.Restart();
        var wg0 = Group.Words(gname, relators);
        GlobalStopWatch.Show($"{wg0.Name}");
    
        GlobalStopWatch.Restart();
        var g0 = Group.SemiDirectProd(new Cn(5), new Cn(4));
        GlobalStopWatch.Show($"{g0.Name}");
    }

    public static void DiCyclic3()
    {
        var wg = Group.Words("Dic3", "a6, b2=a3, bab-1=a-1");
        var g = Group.SemiDirectProd(new Cn(3), new Cn(4));
        DisplayGroup.HeadElements(wg);
        DisplayGroup.HeadElementsSdp(g);
        Console.WriteLine("({0}) IsIsomorphicTo ({1}) : {2}", wg, g, wg.IsIsomorphicTo(g));
    }

    public static void DiCyclic12()
    {
        var gname = "wgDic12";
        var relators = "a12, b2=a6, bab-1=a-1";
        var wg = Group.Words(gname, relators);
        DisplayGroup.HeadElements(wg);
        
        var s8 = new Sn(8);
        var a0 = s8[(1, 2, 3, 4), (5, 6, 7, 8)];
        var b0 = s8[(1, 5, 3, 7), (2, 8, 4, 6)];
        var q8 = Group.Generate("Q8", s8, a0, b0);
        var sp = Group.SemiDirectProd(new Cn(3), q8);
        DisplayGroup.HeadElements(sp);
        Console.WriteLine("({0}) IsIsomorphicTo ({1}) : {2}", wg, sp, wg.IsIsomorphicTo(sp));

        // G:=Group(    (1,2,3,4,5,6,7,8,9,10,11,12)(13,14,15,16,17,18,19,20,21,22,23,24),
        //              (1,17,7,23)(2,16,8,22)(3,15,9,21)(4,14,10,20)(5,13,11,19)(6,24,12,18) );
        var s24 = new Sn(24);
        var a = s24.Cycle(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12) *
                s24.Cycle(13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24);
        var b = s24[(1, 17, 7, 23), (2, 16, 8, 22), (3, 15, 9, 21), (4, 14, 10, 20), (5, 13, 11, 19), (6, 24, 12, 18)];
        var dic12 = Group.Generate("pgDic12", s24, a, b);
        DisplayGroup.HeadElements(dic12);
        Console.WriteLine("({0}) IsIsomorphicTo ({1}) : {2}", wg, dic12, wg.IsIsomorphicTo(dic12));
    }

    public static void Quaternion32()
    {
        var gname = "wgQ32";
        var relators = "a16, b2=a8, bab-1=a-1";
        var wg = Group.Words(gname, relators);
        DisplayGroup.HeadElements(wg);

        var s32 = new Sn(32);
        var a = s32.Cycle(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16) *
                s32.Cycle(17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32);
        var b = s32[(1, 30, 9, 22), (2, 29, 10, 21), (3, 28, 11, 20), (4, 27, 12, 19), (5, 26, 13, 18), (6, 25, 14, 17),
            (7, 24, 15, 32), (8, 23, 16, 31)];
        var q32 = Group.Generate("pgQ32", s32, a, b);
        DisplayGroup.HeadElements(q32);
        Console.WriteLine("({0}) IsIsomorphicTo ({1}) : {2}", wg, q32, wg.IsIsomorphicTo(q32));
    }
}
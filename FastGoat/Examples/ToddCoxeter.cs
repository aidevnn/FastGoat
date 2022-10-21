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
        Console.WriteLine(wg.IsIsomorphicTo(q8));
    
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
        Console.WriteLine(wg.IsIsomorphicTo(gr));
    
        GlobalStopWatch.Restart();
        var wg0 = Group.Words(gname, relators);
        GlobalStopWatch.Show($"{wg0.Name}");
    
        GlobalStopWatch.Restart();
        var g0 = Group.SemiDirectProd(new Cn(5), new Cn(4));
        GlobalStopWatch.Show($"{g0.Name}");
    }
}
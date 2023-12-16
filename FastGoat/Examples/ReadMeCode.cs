using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Perms;
using FastGoat.UserGroup.Polynoms;
using FastGoat.UserGroup.Words;

namespace FastGoat.Examples;

public static class ReadMeCode
{
    public static void ExamplePermGroup21()
    {
        GlobalStopWatch.Restart();
        var s7 = new Sn(7);
        var a = s7[(1, 2, 3, 4, 5, 6, 7)];
        var allC3 = s7.Where(p => (p ^ 3) == s7.Neutral()).ToArray();
        var b = allC3.First(p => Group.GenerateElements(s7, a, p).Count() == 21);
        GlobalStopWatch.Stop();

        Console.WriteLine("|S7|={0}, |{{b in S7 with b^3 = 1}}| = {1}", s7.Count(), allC3.Count());
        Console.WriteLine("First Solution |HK| = 21 : h = {0} and k = {1}", a, b);
        Console.WriteLine();

        var h = Group.Generate("H", s7, a);
        var g21 = Group.Generate("G21", s7, a, b);
        DisplayGroup.Head(g21);
        DisplayGroup.Head(g21.Over(h));
        GlobalStopWatch.Show("Group21");
    }
    
    public static void ExampleWordGroup21()
    {
        GlobalStopWatch.Restart();
        var wg = new WordGroup("a7, b3, a2 = bab-1");
        GlobalStopWatch.Stop();

        DisplayGroup.Head(wg);
        var n = Group.Generate("<a>", wg, wg["a"]);
        DisplayGroup.Head(wg.Over(n));
        GlobalStopWatch.Show($"{wg}");
        Console.WriteLine();
    }
    
    public static void ExampleSdpGroup21()
    {
        GlobalStopWatch.Restart();
        var c7 = new Cn(7);
        var c3 = new Cn(3);
        var g21 = Group.SemiDirectProd(c7, c3);
        GlobalStopWatch.Stop();

        var n = Group.Generate("N", g21, g21[1, 0]);
        DisplayGroup.HeadSdp(g21);
        DisplayGroup.Head(g21.Over(n));
        GlobalStopWatch.Show("Group21");
    }
    
    public static void ExampleCharacterTableGroup21()
    {
        var c7 = new Cn(7);
        var c3 = new Cn(3);
        var g21 = Group.SemiDirectProd(c7, c3);
        FG.CharacterTable(g21).DisplayCells();
    }

    public static void ExampleGaloisGroup21()
    {
        Ring.DisplayPolynomial = MonomDisplay.StarCaret;
        var x = FG.QPoly('X');
        var P = x.Pow(7) - 8 * x.Pow(5) - 2 * x.Pow(4) + 16 * x.Pow(3) + 6 * x.Pow(2) - 6 * x - 2; // GroupNames website
        GaloisApplicationsPart2.GaloisGroupChebotarev(P, details: true);
    }
    
    public static void ExampleGaloisCyclicGroup5()
    {
        Ring.DisplayPolynomial = MonomDisplay.StarCaret;
        var x = FG.QPoly('X');
        var P = x.Pow(5) + x.Pow(4) - 4 * x.Pow(3) - 3 * x.Pow(2) + 3 * x + 1;
        var roots = IntFactorisation.AlgebraicRoots(P, details: true);
        var gal = GaloisTheory.GaloisGroup(roots, details: true);
        DisplayGroup.AreIsomorphics(gal, FG.Abelian(5));
    }

    public static void ExampleGaloisKleinGroup()
    {
        Ring.DisplayPolynomial = MonomDisplay.Caret;
        var x = FG.QPoly('X');
        var (X, _) = FG.EPolyXc(x.Pow(2) - 2, 'a');
        var (minPoly, a0, b0) = IntFactorisation.PrimitiveElt(X.Pow(2) - 3);
        var roots = IntFactorisation.AlgebraicRoots(minPoly);
        Console.WriteLine("Q(√2, √3) = Q(α)");
        var gal = GaloisTheory.GaloisGroup(roots, details: true);
        DisplayGroup.AreIsomorphics(gal, FG.Abelian(2, 2));
    }

    public static void ExampleGaloisDihedralGroup8()
    {
        Ring.DisplayPolynomial = MonomDisplay.Caret;
        var x = FG.QPoly('X');
        var (X, i) = FG.EPolyXc(x.Pow(2) + 1, 'i');
        var (minPoly, _, _) = IntFactorisation.PrimitiveElt(X.Pow(4) - 2);
        var roots = IntFactorisation.AlgebraicRoots(minPoly);
        Console.WriteLine("With α^4-2 = 0, Q(α, i) = Q(β)");
        var gal = GaloisTheory.GaloisGroup(roots, primEltChar: 'β', details: true);
        DisplayGroup.AreIsomorphics(gal, FG.Dihedral(4));
    }
}
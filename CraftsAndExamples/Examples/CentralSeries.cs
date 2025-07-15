using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Matrix;
using FastGoat.UserGroup.Perms;
using FastGoat.UserGroup.Words;

namespace CraftsAndExamples.Examples;

public static class CentralSeries
{
    static void Chains<T>(ConcreteGroup<T> g, bool details = false) where T : struct, IElt<T>
    {
        var comChain = Group.CommutatorsChain(g);
        var zentrumsChain = Group.ZentrumsChain(g);
        var derivedChain = Group.DerivedChain(g);

        if (details)
        {
            Console.WriteLine("############### Upper Central Series");
            foreach (var z in zentrumsChain)
                DisplayGroup.Head(z);

            Console.WriteLine("############### Lower Central Series");
            foreach (var d in comChain)
                DisplayGroup.Head(d);

            Console.WriteLine("############### Derived Series");
            foreach (var d in derivedChain)
                DisplayGroup.Head(d);
        }

        if (comChain.Last().Count() == 1)
            Console.WriteLine($"{g} is Nilpotent");
        else if (derivedChain.Last().Count() == 1)
            Console.WriteLine($"{g} is soluble");

        var hc = zentrumsChain.Last();
        if (hc.Count() == 1)
            Console.WriteLine($"{g} is centerless");
        else if (hc.Count() != g.Count() && hc.Count() != 1)
            Console.WriteLine($"{g} has hypercenter {hc} of order {hc.Count()}");
    }

    public static void NonAbelian2Groups()
    {
        var d8 = new WordGroup("D8", "a4, b2, abab");
        var q8 = new WordGroup("Q8", "a4, a2 = b2, abab-1");

        Chains(d8, details: true);
        Chains(q8, details: true);

        Chains(new WordGroup("D16", "a8, b2, abab"));
        Chains(new WordGroup("D32", "a16, b2, abab"));
        Chains(new WordGroup("Q16", "a8, a4 = b2, abab-1"));
        Chains(new WordGroup("Q32", "a16, a8 = b2, abab-1"));
    }

    public static void NonAbelian3Groups()
    {
        var c3 = new Cn(3);
        var c9 = new Cn(9);
        var g1 = new WordGroup("C9 x: C3", "a9, b3, a4 = bab-1");
        var g2 = Group.SemiDirectProd(Product.Generate(c3, c3), c3);
        var g3 = new WordGroup("C9 x: C9", "a9, b9, a4 = bab-1");
        var g4 = Group.SemiDirectProd(Product.Generate(c9, c3), c3);
        DisplayGroup.HeadOrders(g1);
        DisplayGroup.HeadOrders(g2);
        DisplayGroup.HeadOrders(g3);
        DisplayGroup.HeadOrders(g4);

        Chains(g1, details: true);
        Chains(g2, details: true);
        Chains(g3);
        Chains(g4);
    }

    public static void NonAbelianOrder24()
    {
        var gl = new GL(2, 3);

        var s4 = new Symm(4);
        var e24 = new WordGroup("E24", "a4, b2, c3, bab = a3, bcb = c, aca3 = c2");
        var sl23 = Group.Generate("SL23", gl, gl[1, 1, 0, 1], gl[0, 1, 2, 0]);
        var c3q8 = new WordGroup("C3 x Q8", "a4, b2=a2, bab-1=a-1, c3, ac = ca, bc = cb");
        DisplayGroup.HeadOrders(s4);
        DisplayGroup.HeadOrders(e24);
        DisplayGroup.HeadOrders(sl23);
        DisplayGroup.HeadOrders(c3q8);

        Chains(s4);
        Chains(e24);
        Chains(sl23);
        Chains(c3q8);
    }

    static void UpperSeriesFast<T>(ConcreteGroup<T> gr) where T : struct, IElt<T>
    {
        Console.WriteLine($"######## Start {gr}");
        GlobalStopWatch.Restart();
        foreach (var g in Group.ZentrumsChain(gr))
            DisplayGroup.Head(g);

        GlobalStopWatch.Show("Definition");
        Console.WriteLine("######## Fast");
        GlobalStopWatch.Restart();
        foreach (var g in Group.ZentrumsChainFast(gr))
            DisplayGroup.Head(g);

        GlobalStopWatch.Show("Fast");
        Console.WriteLine("######## End");
        Console.WriteLine();
    }

    public static void UpperSeries()
    {
        var c3 = new Cn(3);
        var c9 = new Cn(9);

        var g1 = new WordGroup("C9 x: C3", "a9, b3, a4 = bab-1");
        var g2 = Group.SemiDirectProd(Product.Generate(c3, c3), c3);
        var g3 = Group.SemiDirectProd(Product.Generate(c9, c3), c3);
        var d8 = new WordGroup("D8", "a4, b2, abab");
        var d12 = new WordGroup("D12", "a6, b2, abab");
        var e24 = new WordGroup("E24", "a4, b2, c3, bab = a3, bcb = c, aca3 = c2");

        UpperSeriesFast(g1);
        UpperSeriesFast(g2);
        UpperSeriesFast(g3);
        UpperSeriesFast(d8);
        UpperSeriesFast(d12);
        UpperSeriesFast(e24);
    }
}
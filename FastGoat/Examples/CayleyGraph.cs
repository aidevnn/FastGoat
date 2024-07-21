using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Words.Tools;

namespace FastGoat.Examples;

public static class CayleyGraph
{
    public static void ExamplesCayleyGraph()
    {
        DisplayGroup.HeadElementsCayleyGraph(FG.AbelianPerm(5));
        DisplayGroup.HeadElementsCayleyGraph(FG.Symmetric(3));
        DisplayGroup.HeadElementsCayleyGraph(FG.Dihedral(4));
        DisplayGroup.HeadElementsCayleyGraph(FG.QuaternionWg(8));
        DisplayGroup.HeadElementsCayleyGraph(FG.MetaCyclicSdpWg(7, 3, 2)); 
    }

    public static void GroupsUptoOrder12()
    {
        foreach (var g in 16.Range(1).SelectMany(o => FG.AllGroupsOfOrder(o)))
            DisplayGroup.HeadElementsCayleyGraph(g);
    }

    public static void Q8ThreeGenerators()
    {
        var q8a = FG.DiCyclicThreeGens(2);

        var (a, b, c) = q8a.GetGenerators().Deconstruct();
        Console.WriteLine($"{q8a}:{q8a.Definition}");
        Console.WriteLine($"b2 = {q8a.Times(b, 2)} and c2 = {q8a.Times(c, 2)}");
        Console.WriteLine($"ab = {q8a.Op(a, b)}  and (ba)-1 = {q8a.Invert(q8a.Op(b, a))}");
        Console.WriteLine($"bc = {q8a.Op(b, c)}  and (cb)-1 = {q8a.Invert(q8a.Op(c, b))}");
        Console.WriteLine($"ca = {q8a.Op(c, a)}  and (ac)-1 = {q8a.Invert(q8a.Op(a, c))}");
        Console.WriteLine();
        
        DisplayGroup.HeadElementsCayleyGraph(q8a);
    }

    public static void A4ThreeGenerators()
    {
        var a4 = FG.Alternate(4);
        DisplayGroup.HeadElementsCayleyGraph(a4);

        var gens = a4.Where(e => a4.ElementsOrders[e] == 2).Take(2)
            .Append(a4.First(e => a4.ElementsOrders[e] == 3))
            .ToArray();

        DisplayGroup.HeadElementsCayleyGraph(a4, gens: gens);

        var a4wg = FG.WordGroup("A4", Graph.DefiningRelatorsOfGroup(a4, gens));
        Console.WriteLine($"{a4wg}:{a4wg.Definition}");
        DisplayGroup.HeadElementsCayleyGraph(a4wg);
    }

    public static void S4CayleyGraph1()
    {
        var s4 = FG.Symmetric(4);
        DisplayGroup.HeadElementsCayleyGraph(s4);

        var (c2, c4) = (s4.First(e => s4.ElementsOrders[e] == 2), s4.First(e => s4.ElementsOrders[e] == 4));
        var s4wg = FG.WordGroup("S4", Graph.DefiningRelatorsOfGroup(s4, [c2, c4]));
        Console.WriteLine($"{s4wg}:{s4wg.Definition}");
        DisplayGroup.HeadElementsCayleyGraph(s4wg);
    }

    public static void S4CayleyGraph2()
    {
        var s4 = FG.Symmetric(4);
        DisplayGroup.HeadElementsCayleyGraph(s4);

        var (c3, c4) = (s4.First(e => s4.ElementsOrders[e] == 3), s4.First(e => s4.ElementsOrders[e] == 4));
        var s4wg = FG.WordGroup("S4", Graph.DefiningRelatorsOfGroup(s4, [c3, c4]));
        Console.WriteLine($"{s4wg}:{s4wg.Definition}");
        DisplayGroup.HeadElementsCayleyGraph(s4wg);
    }
}
using System.Security.Cryptography;
using FastGoat.Gp;
using FastGoat.UserGroup;

namespace FastGoat.Examples;

public static class AbelianInvariantsFactors
{
    public static Stack<int> Reduce<T>(ConcreteGroup<T> g) where T : struct, IElt<T>
    {
        if (g.GroupType == GroupType.NonAbelianGroup)
            throw new Exception("Only Abelian Group");

        var g0 = g;
        Stack<int> facts = new Stack<int>();
        while (g0.Count()!=1)
        {
            var p = g0.ElementsOrders.OrderByDescending(e => e.Value).ThenBy(e => e.Key).First();
            var h = Group.Generate(g0, p.Key);
            g0 = g0.Over(h);
            facts.Push(p.Value);
        }

        return facts;
    }
    
    public static void InvariantFactors294()
    {
        var c14 = new Cn(14);
        var c21 = new Cn(21);
        
        var bg = Product.Group(c14, c21);
        var g = Group.Create(bg.Name, bg);
        var decomposition = Reduce(g);
        Console.WriteLine("{0} ~ {1}", g, decomposition.Glue(" x ", "C{0}"));
    }

    public static void InvariantFactors600()
    {
        var c20 = new Cn(20);
        var c30 = new Cn(30);
        
        var bg = Product.Group(c20, c30);
        var g = Group.Create(bg.Name, bg);
        var decomposition = Reduce(g);
        Console.WriteLine("{0} ~ {1}", g, decomposition.Glue(" x ", "C{0}"));
    }

    public static void InvariantFactors4320()
    {
        var c8 = new Cn(8);
        var c18 = new Cn(18);
        var c30 = new Cn(30);
        
        var bg = Product.Group(c8, c18, c30);
        var g = Group.Create(bg.Name, bg);
        var decomposition = Reduce(g);
        Console.WriteLine("{0} ~ {1}", g, decomposition.Glue(" x ", "C{0}"));
    }
}
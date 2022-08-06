using FastGoat.Structures.SetTheory;

namespace FastGoat.Structures.GroupTheory;

public static partial class GroupExt
{
    public static List<U> InvariantsFactors<U>(SubGroup<U> gr) where U : struct, IElt<U>
    {
        if (!gr.IsCommutative())
        {
            Console.WriteLine("Abelian Group Only");
            return new List<U>();
        }

        var G = new GroupElement<U>(gr, gr.AllElements);
        var H = gr.Singleton();
        QuotientGroup<U> k0 = G.Over(H);
        k0.SetName("G");
        Console.WriteLine($"Invariants factors of G = {gr.Infos.Name}");
        k0.Infos.SetDetails(name: "G");
        k0.DisplayHead();
        List<string> chainFactors = new();
        List<U> allFactors = new();

        while (true)
        {
            k0.ComputeOrders();
            var Fact0 = k0.AllElements.OrderByDescending(e => k0.GetOrder(e)).First();
            allFactors.Add(Fact0);
            var C0 = k0.Monogenic(Fact0);
            chainFactors.Add($"C{C0.Count}");
            C0.SetName($"C{C0.Count}");
            Console.WriteLine($"C{k0.GetOrder(Fact0)} = {Fact0}; |<C{C0.Count}>|={C0.Count}");
            Console.WriteLine("{0} is SubGroup of {1} : {2}", C0.Infos.Name, k0.Infos.Name, C0.AllElements.All(k0.Contains));

            k0 = k0.Over(C0);
            k0.DisplayHead();

            if (k0.Count == 1) break;
        }

        chainFactors.Reverse();
        Console.WriteLine("{0} = G[{1}] ~ {2}", gr.Infos.Name, G.Count, chainFactors.Glue(sep: " x "));
        Console.WriteLine("-----------------------------");

        return allFactors;
    }

}

using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Words;

namespace Examples;

public static class HolomorphC7
{
    public static void Sample()
    {
        // Saunders MacLane, Garrett Birkhoff. Algebra (3rd ed.) page 416
        var c7 = new Cn(7);
        var c6 = new Cn(6);
        var thetas = Group.AllOpsByAutomorphisms(c6, c7);
        List<SemiDirectProduct<ZnInt, ZnInt>> nonIsomorphics = new();
        char k = 'a';
        var nm = "(C7 x: C6)";
        foreach (var theta in thetas)
        {
            var sdp = Group.SemiDirectProd($"{nm}{k++}", c7, theta, c6);
            if (nonIsomorphics.All(g => !g.IsIsomorphicTo(sdp)))
                nonIsomorphics.Add(sdp);
        }

        nonIsomorphics.Sort((a, b) => -a.ElementsOrdersList().Max().CompareTo(b.ElementsOrdersList().Max()));
        foreach (var sdp in nonIsomorphics)
        {
            DisplayGroup.HeadSdp(sdp);
            Console.WriteLine(sdp.ElementsOrdersList().Glue(", "));
            Console.WriteLine();
        }

        var c7c6 = Product.Generate(c6, c7);
        DisplayGroup.Head(c7c6);
        Console.WriteLine(c7c6.ElementsOrdersList().Glue(", "));
        Console.WriteLine();

        var wgC3D14 = new WordGroup("C3 x D14", "a7, b2, c3, abab, ac = ca, bc = cb");
        DisplayGroup.Head(wgC3D14);
        Console.WriteLine(wgC3D14.ElementsOrdersList().Glue(", "));
        Console.WriteLine();

        var wgC2H21 = new WordGroup("C2 x H21", "a7, b3, c2, a2 = bab-1, ac = ca, bc = cb");
        DisplayGroup.Head(wgC2H21);
        Console.WriteLine(wgC2H21.ElementsOrdersList().Glue(", "));
        Console.WriteLine();

        var u7 = new Un(7);
        var a = u7[(c7[1], c7[3])];
        var theta1 = Group.HomomorphismMap(u7, u7, new() { [a] = a });
        var homTheta1 = Group.Hom(u7, theta1);
        var hol7 = Group.SemiDirectProd("Hol7", c7, homTheta1, u7);
        DisplayGroup.HeadSdp(hol7);
        Console.WriteLine(hol7.ElementsOrdersList().Glue(", "));
        Console.WriteLine();

        var g0 = nonIsomorphics[0];
        Console.WriteLine("{0} IsIsomorphicTo {1} : {2}", c7c6, g0, c7c6.IsIsomorphicTo(g0));
        var g1 = nonIsomorphics[1];
        Console.WriteLine("{0} IsIsomorphicTo {1} : {2}", wgC3D14, g1, wgC3D14.IsIsomorphicTo(g1));
        var g2 = nonIsomorphics[2];
        Console.WriteLine("{0} IsIsomorphicTo {1} : {2}", wgC2H21, g2, wgC2H21.IsIsomorphicTo(g2));
        var g3 = nonIsomorphics[3];
        Console.WriteLine("{0} IsIsomorphicTo {1} : {2}", hol7, g3, hol7.IsIsomorphicTo(g3));
    }
}
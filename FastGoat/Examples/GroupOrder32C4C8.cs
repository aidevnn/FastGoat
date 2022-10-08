using FastGoat;
using FastGoat.Gp;
using FastGoat.UserGroup;

namespace FastGoat.Examples;

public static class GroupOrder32C4C8
{
    private static Cn C4 { get; } = new Cn(4);
    private static Un U4 { get; } = new Un(4);
    private static Cn C8 { get; } = new Cn(8);
    private static Un U8 { get; } = new Un(8);

    public static void NonIsomorphic()
    {
        var yC8C4 = U8[(C8[1], C8[5])];
        var mapC8C4 = Group.PartialMap((C4[1], yC8C4));
        var thetaC8C4 = Group.HomomorphismMap(C4, U8, mapC8C4);
        var gC8C4 = Group.SemiDirectProd(C8, thetaC8C4, C4);
        DisplayGroup.HeadSdp(gC8C4);
        Console.WriteLine("Sorted Orders {0} : ({1})", gC8C4, gC8C4.ElementsOrders.Values.Ascending().Glue(", "));

        var yC4C8 = U4[(C4[1], C4[3])];
        var mapC4C8 = Group.PartialMap((C8[1], yC4C8));
        var thetaC4C8 = Group.HomomorphismMap(C8, U4, mapC4C8);
        var gC4C8 = Group.SemiDirectProd(C4, thetaC4C8, C8);
        DisplayGroup.HeadSdp(gC4C8);
        Console.WriteLine("Sorted Orders {0} : ({1})", gC4C8, gC4C8.ElementsOrders.Values.Ascending().Glue(", "));

        Console.WriteLine();
        Console.WriteLine("Isomorphic {0} {1} {2}", gC4C8, gC8C4, gC4C8.IsIsomorphicTo(gC8C4));

        var homs = Group.AllHomomorphisms(gC8C4, gC4C8, true);
        Console.WriteLine("Homomorphisms Count = {0}", homs.Count);
        Console.WriteLine("Isomorphisms  Count = {0}", homs.Count(h => h.Values.Distinct().Count() == h.Count));
    }
}
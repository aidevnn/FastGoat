using FastGoat;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.UserGroup.Integers;

namespace FastGoat.Examples;

public static class GroupOrder32C4C8
{
    private static Cn C4 { get; } = new Cn(4);
    private static Un U4 { get; } = new Un(4);
    private static Cn C8 { get; } = new Cn(8);
    private static Un U8 { get; } = new Un(8);
    private static Cn C12 { get; } = new Cn(12);

    public static void NonIsomorphic()
    {
        var yC8C4 = U8[(C8[1], C8[5])];
        var mapC8C4 = Group.PartialMap((C4[1], yC8C4));
        var thetaC8C4 = Group.HomomorphismMap(C4, U8, mapC8C4);
        var homThetaC8C4 = Group.Hom(C4, thetaC8C4);
        var gC8C4 = Group.SemiDirectProd(C8, homThetaC8C4, C4);
        DisplayGroup.HeadSdp(gC8C4);
        Console.WriteLine("Sorted Orders {0} : ({1})", gC8C4, gC8C4.ElementsOrders.Values.Ascending().Glue(", "));

        var yC4C8 = U4[(C4[1], C4[3])];
        var mapC4C8 = Group.PartialMap((C8[1], yC4C8));
        var thetaC4C8 = Group.HomomorphismMap(C8, U4, mapC4C8);
        var homThetaC4C8 = Group.Hom(C8, thetaC4C8);
        var gC3C4 = Group.SemiDirectProd(C4, homThetaC4C8, C8);
        DisplayGroup.HeadSdp(gC3C4);
        Console.WriteLine("Sorted Orders {0} : ({1})", gC3C4, gC3C4.ElementsOrders.Values.Ascending().Glue(", "));

        Console.WriteLine();
        Console.WriteLine("Isomorphic {0} ~ {1} {2}", gC3C4, gC8C4, gC3C4.IsIsomorphicTo(gC8C4));

        var homs = Group.AllHomomorphisms(gC8C4, gC3C4);
        Console.WriteLine("Homomorphisms Count = {0}", homs.Count);
        Console.WriteLine("Isomorphisms  Count = {0}", homs.Count(h => h.Image().Count() == h.Count));
    }

    public static void NonIsomorphicAnotherExample()
    {
        var C2 = new Cn(2);
        var yC8C2 = U8[(C8[1], C8[5])];
        var mapC8C2 = Group.PartialMap((C2[1], yC8C2));
        var thetaC8C2 = Group.HomomorphismMap(C2, U8, mapC8C2);
        var homThetaC8C2 = Group.Hom(C2, thetaC8C2);
        var gC8C2 = Group.SemiDirectProd(C8, homThetaC8C2, C2);
        DisplayGroup.HeadSdp(gC8C2);
        var oC8C2 = gC8C2.ElementsOrders.Values.Ascending().ToArray();
        Console.WriteLine("Sorted Orders {0} : ({1})", gC8C2, oC8C2.Glue(", "));

        var gdC8C2 = Product.Generate(C8, C2);
        DisplayGroup.Head(gdC8C2);
        var odC8C2 = gdC8C2.ElementsOrders.Values.Ascending().ToArray();
        Console.WriteLine("Sorted Orders {0} : ({1})", gdC8C2, odC8C2.Glue(", "));
        Console.WriteLine("Sorted Order Equal : {0}", oC8C2.SequenceEqual(odC8C2));

        Console.WriteLine();
        Console.WriteLine("Isomorphic {0} ~ {1} {2}", gdC8C2, gC8C2, gdC8C2.IsIsomorphicTo(gC8C2));

        var homs = Group.AllHomomorphisms(gC8C2, gdC8C2);
        Console.WriteLine("Homomorphisms Count = {0}", homs.Count);
        Console.WriteLine("Isomorphisms  Count = {0}", homs.Count(h => h.Image().Count() == h.Count));
    }

    public static void Isomorphic()
    {
        var gC3C4 = Product.Generate(new Cn(3), C4);
        DisplayGroup.Head(C12);
        Console.WriteLine("Sorted Orders {0} : ({1})", C12, C12.ElementsOrders.Values.Ascending().Glue(", "));
        DisplayGroup.Head(gC3C4);
        Console.WriteLine("Sorted Orders {0} : ({1})", gC3C4, gC3C4.ElementsOrders.Values.Ascending().Glue(", "));

        Console.WriteLine();
        Console.WriteLine("Isomorphic {0} ~ {1} {2}", gC3C4, C12, gC3C4.IsIsomorphicTo(C12));

        var homs = Group.AllHomomorphisms(C12, gC3C4);
        Console.WriteLine("Homomorphisms Count = {0}", homs.Count);
        Console.WriteLine("Isomorphisms  Count = {0}", homs.Count(h => h.Image().Count() == h.Count));
    }

    public static void IsomorphicAnotherExample()
    {
        var gC8C12 = Product.Generate(C8, C12);
        var gC4C24 = Product.Generate(C4, new Cn(24));
        DisplayGroup.Head(gC4C24);
        Console.WriteLine("Sorted Orders {0} : ({1})", gC4C24, gC4C24.ElementsOrders.Values.Ascending().Glue(", "));
        DisplayGroup.Head(gC8C12);
        Console.WriteLine("Sorted Orders {0} : ({1})", gC8C12, gC8C12.ElementsOrders.Values.Ascending().Glue(", "));

        Console.WriteLine();
        Console.WriteLine("Isomorphic {0} ~ {1} {2}", gC8C12, gC4C24, gC8C12.IsIsomorphicTo(gC4C24));

        var homs = Group.AllHomomorphisms(gC4C24, gC8C12);
        Console.WriteLine("Homomorphisms Count = {0}", homs.Count);
        Console.WriteLine("Isomorphisms  Count = {0}", homs.Count(h => h.Image().Count() == h.Count));
    }
}
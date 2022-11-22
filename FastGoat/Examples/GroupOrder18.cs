using FastGoat.Structures.CartesianProduct;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Perms;
using FastGoat.UserGroup.Words;

namespace FastGoat.Examples;

public static class GroupOrder18
{
    public static void Abelians18()
    {
        var c18 = new Cn(18);
        DisplayGroup.Head(c18);
        Console.WriteLine("Elements Orders : {0}", c18.ElementsOrdersList().Glue(", "));
        var c3c6 = Product.Generate(new Cn(3), new Cn(6));
        DisplayGroup.Head(c3c6);
        Console.WriteLine("Elements Orders : {0}", c3c6.ElementsOrdersList().Glue(", "));
    }

    public static void C3xS3()
    {
        var c3s6 = Product.Generate(new Cn(3), new Symm(3));
        DisplayGroup.Head(c3s6);
        Console.WriteLine("Elements Orders : {0}", c3s6.ElementsOrdersList().Glue(", "));

        var s6 = new Sn(6);
        var a = s6[(1, 2, 3)];
        var b = s6[(1, 2)];
        var c = s6[(4, 5, 6)];
        var pg = Group.Generate("pg(C3 x S3)", s6, a, b, c);
        DisplayGroup.Head(pg);
        Console.WriteLine("Elements Orders : {0}", pg.ElementsOrdersList().Glue(", "));
        Console.WriteLine("({0}) IsIsomorphicTo ({1}) : {2}", c3s6, pg, c3s6.IsIsomorphicTo(pg));

        var wg = new WordGroup("wg(C3 x S3)", "a3, b3, c2, ab=ba, ac=ca, bcbc");
        DisplayGroup.Head(wg);
        Console.WriteLine("Elements Orders : {0}", wg.ElementsOrdersList().Glue(", "));
        Console.WriteLine("({0}) IsIsomorphicTo ({1}) : {2}", wg, pg, wg.IsIsomorphicTo(pg));
    }

    public static void D18()
    {
        var n = new Cn(9);
        var g = new Cn(2);
        var autN = Group.AutomorphismGroup(n);
        var inv = autN[(n[1], n[8])];
        var theta = Group.HomomorphismMap(g, autN, new() { [g[1]] = inv });
        var homTheta = Group.Hom(g, theta);
        var sdpD18 = Group.SemiDirectProd("sdp(D18)", n, homTheta, g);
        DisplayGroup.HeadSdp(sdpD18);
        Console.WriteLine("Elements Orders : {0}", sdpD18.ElementsOrdersList().Glue(", "));
        Console.WriteLine();

        var s9 = new Sn(9);
        var a = s9[(1, 2, 3, 4, 5, 6, 7, 8, 9)];
        var b = s9[(2, 9), (3, 8), (4, 7), (5, 6)];
        var pgD18 = Group.Generate("pg(D18)", s9, a, b);
        DisplayGroup.Head(pgD18);
        Console.WriteLine("Elements Orders : {0}", pgD18.ElementsOrdersList().Glue(", "));
        Console.WriteLine("({0}) IsIsomorphicTo ({1}) : {2}", sdpD18, pgD18, sdpD18.IsIsomorphicTo(pgD18));
        Console.WriteLine();

        var wg = new WordGroup("wg(D18)", "a9, b2, abab");
        DisplayGroup.Head(wg);
        Console.WriteLine("Elements Orders : {0}", wg.ElementsOrdersList().Glue(", "));
        Console.WriteLine("({0}) IsIsomorphicTo ({1}) : {2}", wg, pgD18, wg.IsIsomorphicTo(pgD18));
    }

    public static void C3sdpS3()
    {
        var n = Product.Generate(new Cn(3), new Cn(3));
        var g = new Cn(2);
        var autN = Group.AutomorphismGroup(n);
        var inv = autN[(n[1, 0], n[2, 0]), (n[0, 1], n[0, 2])];
        var theta = Group.HomomorphismMap(g, autN, new() { [g[1]] = inv });
        var homTheta = Group.Hom(g, theta);
        var sdp = Group.SemiDirectProd("sdp((C3 x C3) x: C2)", n, homTheta, g);
        DisplayGroup.HeadSdp(sdp);
        Console.WriteLine("Elements Orders : {0}", sdp.ElementsOrdersList().Glue(", "));
        Console.WriteLine();

        var wg = new WordGroup("wg(C3 x: S3)", "a3, b3, c2, ab=ba, cac=a-1, cbc=b-1");
        DisplayGroup.Head(wg);
        Console.WriteLine("Elements Orders : {0}", wg.ElementsOrdersList().Glue(", "));
        Console.WriteLine("({0}) IsIsomorphicTo ({1}) : {2}", wg, sdp, wg.IsIsomorphicTo(sdp));
        Console.WriteLine();

        var sdp2 = Group.SemiDirectProd("sdp(C3 x: S3)", new Cn(3), new Symm(3));
        DisplayGroup.HeadSdp(sdp2);
        Console.WriteLine("Elements Orders : {0}", sdp2.ElementsOrdersList().Glue(", "));
        Console.WriteLine("({0}) IsIsomorphicTo ({1}) : {2}", wg, sdp2, wg.IsIsomorphicTo(sdp2));
        Console.WriteLine();

        var s9 = new Sn(9);
        var a = s9[(1, 2, 3), (4, 5, 6), (7, 8, 9)];
        var b = s9[(1, 5, 8), (2, 6, 9), (3, 4, 7)];
        var c = s9[(2, 3), (4, 9), (5, 8), (6, 7)];
        var pg = Group.Generate("pg(C3 x: S3)", s9, a, b, c);
        DisplayGroup.Head(pg);
        Console.WriteLine("Elements Orders : {0}", pg.ElementsOrdersList().Glue(", "));
        Console.WriteLine("({0}) IsIsomorphicTo ({1}) : {2}", wg, pg, wg.IsIsomorphicTo(pg));
    }
}
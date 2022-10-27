using FastGoat.UserGroup;

namespace FastGoat.Examples;

public static class PSLnp
{
    public static void L23()
    {
        var gl = new GL(2, 3);
        var a = gl[1, 1, 0, 1];
        var b = gl[0, 1, 2, 0];

        var sl23 = Group.Generate("SL2(3)", gl, a, b);
        DisplayGroup.Head(sl23);
        var zg = Group.Zentrum(sl23);
        DisplayGroup.Head(zg);
        var l23 = sl23.Over(zg, "L2(3)");
        DisplayGroup.Head(l23);

        var sn = new Sn(4);
        var a4 = Group.Generate("A4", sn, sn[(1, 2, 3)], sn[(1, 2), (3, 4)]);
        DisplayGroup.Head(a4);
        Console.WriteLine("{0} IsIsomorphicTo {1} : {2}",l23, a4, l23.IsIsomorphicTo(a4));
    }
    
    public static void L27_L32()
    {
        var gl27 = new GL(2, 7);
        var a = gl27[3, 0, 0, 5];
        var b = gl27[6, 1, 6, 0];

        var sl27 = Group.Generate("SL2(7)", gl27, a, b);
        DisplayGroup.Head(sl27);
        var zg27 = Group.Zentrum(sl27);
        DisplayGroup.Head(zg27);
        var l27 = sl27.Over(zg27, "L2(7)");
        DisplayGroup.Head(l27);

        var gl32 = new GL(3, 2);
        var a3 = gl32.At((0, 4, 8, 1), 1);
        var b3 = gl32.At((2, 3, 7), 1);

        var sl32 = Group.Generate("SL3(2)", gl32, a3, b3);
        DisplayGroup.Head(sl32);
        var zg32 = Group.Zentrum(sl32);
        DisplayGroup.Head(zg32);
        var l32 = sl32.Over(zg32, "L3(2)");
        DisplayGroup.Head(l32);
        
        var sn = new Sn(7);
        var a2 = sn[(1, 2), (4, 5)];
        var b2 = sn[(2, 3, 4), (5, 6, 7)];
        var pg = Group.Generate("pg(GL3(2))", a2, b2);
        DisplayGroup.Head(pg);
        
        DisplayGroup.AreIsomorphics(pg, l32);
        DisplayGroup.AreIsomorphics(l27, l32);
    }
}
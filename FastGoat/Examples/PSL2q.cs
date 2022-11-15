using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Matrix;
using FastGoat.UserGroup.Perms;

namespace FastGoat.Examples;

public static class PSL2q
{
    public static void L2_4()
    {
        var gl = new GLnq(2, 4);
        var a = gl['x', 0, 0, 1];
        var b = gl[1, 1, 1, 0];

        var gl2_4 = Group.Generate(gl, a, b);
        DisplayGroup.HeadOrders(gl2_4);

        var one = gl.Fq.One;
        var detPositive = gl2_4.Where(e => gl.Determinant(e).Equals(one)).ToArray();
        Console.WriteLine("SL count :{0}", detPositive.Length);
        var a1 = gl['x', 0, 0, ('x', 2)];
        var sl2_4 = Group.Generate("SL2(4)", gl, a1, b);
        DisplayGroup.HeadOrders(sl2_4);

        var zg2_4 = Group.Zentrum(sl2_4);
        DisplayGroup.Head(zg2_4);
        var l2_4 = sl2_4.Over(zg2_4, "L2(4)");
        DisplayGroup.Head(l2_4);
    }

    public static void L2_8()
    {
        var gl = new GLnq(2, 8);
        var a = gl['x', 0, 0, 1];
        var b = gl[1, 1, 1, 0];

        var gl2_8 = Group.Generate(gl, a, b);
        DisplayGroup.HeadOrders(gl2_8); // |Gl(2, F8)| = 3528

        var one = gl.Fq.One;
        var detPositive = gl2_8.Where(e => gl.Determinant(e).Equals(one)).ToArray();
        Console.WriteLine("SL count :{0}", detPositive.Length);
        var a1 = gl['x', 0, 0, ('x', 6)];
        var sl2_8 = Group.Generate("SL2(8)", gl, a1, b);
        DisplayGroup.HeadOrders(sl2_8);

        var zg2_8 = Group.Zentrum(sl2_8);
        DisplayGroup.Head(zg2_8);
        var l2_8 = sl2_8.Over(zg2_8, "L2(8)");
        DisplayGroup.Head(l2_8);
    }

    public static void L2_9()
    {
        var gl = new GLnq(2, 9);
        var a = gl['x', 0, 0, 1];
        var b = gl[2, 1, 2, 0];

        var gl2_9 = Group.Generate(gl, a, b);
        DisplayGroup.HeadOrders(gl2_9); // |Gl(2, F9)| = 5760
        
        var a1 = gl['x', 0, 0, ('x', 7)];
        var sl2_9 = Group.Generate("SL2(9)", gl, a1, b);
        DisplayGroup.HeadOrders(sl2_9); // |SL2(9)| = 720

        var s6 = new Symm(6);
        DisplayGroup.HeadOrders(s6);
        DisplayGroup.AreIsomorphics(s6, sl2_9);

        var zg2_9 = Group.Zentrum(sl2_9);
        var l2_9 = sl2_9.Over(zg2_9, "L2(9)");
        DisplayGroup.HeadOrders(l2_9);
        
        var a6 = Group.Generate("A6", s6, s6[(4, 5, 6)], s6[(1, 2, 3, 4, 5)]);
        DisplayGroup.HeadOrders(a6);

        DisplayGroup.AreIsomorphics(l2_9, a6);
    }

    public static void L2_16()
    {
        var gl = new GLnq(2, 16);
        var a = gl['x', 0, 0, 1];
        var b = gl[1, 1, 1, 0];

        // var gl216 = Group.Generate(gl, a, b);
        // DisplayGroup.HeadOrders(gl216); // |Gl(2, F16)| = 61200

        var a1 = gl['x', 0, 0, ('x', 14)];
        var sl2_16 = Group.Generate("SL2(16)", gl, a1, b);
        DisplayGroup.HeadOrders(sl2_16); // |SL2(16)| = 4080

        var zg2_16 = Group.Zentrum(sl2_16);
        DisplayGroup.Head(zg2_16);
        var l2_16 = sl2_16.Over(zg2_16, "L2(16)");
        DisplayGroup.Head(l2_16); // |L2(16)| = 4080
    }

    public static void L2_25()
    {
        var gl = new GLnq(2, 25);

        var a = gl['x', 0, 0, ('x', 23)];
        var b = gl[4, 1, 4, 0];
        var sl2_25 = Group.Generate("SL2(25)", gl, a, b);
        DisplayGroup.HeadOrders(sl2_25); // |SL2(25)| = 15600


        var zg2_25 = Group.Zentrum(sl2_25);
        DisplayGroup.Head(zg2_25); // |Z(SL2(25))| = 2
        var l2_25 = sl2_25.Over(zg2_25, "L2(25)");
        DisplayGroup.Head(l2_25); // |L2(25)| = 7800
    }

    public static void L2_27()
    {
        var gl = new GLnq(2, 27);

        var a = gl['x', 0, 0, ('x', 25)];
        var b = gl[2, 1, 2, 0];
        var sl2_27 = Group.Generate("SL2(27)", gl, a, b);
        DisplayGroup.HeadOrders(sl2_27); // |SL2(27)| = 19656
        
        var zg2_27 = Group.Zentrum(sl2_27);
        DisplayGroup.Head(zg2_27); // |Z(SL2(27))| = 2
        var l2_27 = sl2_27.Over(zg2_27, "L2(27)");
        DisplayGroup.Head(l2_27); // |L2(27)| = 9828
    }

}
using FastGoat.Structures;
using FastGoat.UserGroup.Matrix;

namespace FastGoat.Examples;

public static class UnitaryGroup
{
    public static void U3_2()
    {
        var gl34 = new GLnq(3, 4);
        var x = gl34.Fq['x'];
        var a = gl34[
            1, x, x,
            0, 1, x + 1,
            0, 0, 1
        ];
        var b = gl34[
            x, 1, 1,
            1, 1, 0,
            1, 0, 0
        ];

        var u3_2 = Group.Generate("U3(2)", gl34, a, b);
        DisplayGroup.Head(u3_2);
    }

    public static void U3_3()
    {
        var gl39 = new GLnq(3, 9);
        var x = gl39.Fq['x'];
        var a = gl39[
            x + 2, 1, 1,
            1, 2, 0,
            1, 0, 0
        ];
        var b = gl39[
            2 * x + 1, 2 * x + 1, 1,
            x, 2, 0,
            1, 0, 0
        ];

        var u3_3 = Group.Generate("U3(3)", gl39, a, b);
        DisplayGroup.Head(u3_3); // |U3(3)| = 6048
    }

    public static void U3_4()
    {
        var gl316 = new GLnq(3, 16);
        var x = gl316.Fq['x'];
        var x3 = x.Pow(3);
        var x11 = x.Pow(3) + x.Pow(2) + x;
        var a = gl316[
            x, 0, 0,
            0, x3, 0,
            0, 0, x11
        ];
        var b = gl316[
            x, 1, 1,
            1, 1, 0,
            1, 0, 0
        ];

        var u3_4 = Group.Generate("U3(4)", gl316, a, b);
        DisplayGroup.Head(u3_4); // |U3(4)| = 62400
    }

    public static void U4_2()
    {
        var gl44 = new GLnq(4, 4);
        var x = gl44.Fq['x'];

        var a = gl44[
            x, 0, 0, 0,
            0, x * x, 0, 0,
            0, 0, x * x, 0,
            0, 0, 0, x];
        var b = gl44[
            1, 0, 1, 0,
            1, 0, 0, 0,
            0, 1, 0, 1,
            0, 1, 0, 0];

        var u4_2 = Group.Generate("U4(2)", gl44, a, b);
        DisplayGroup.HeadOrders(u4_2);
        
        var gl43 = new GL(4, 3);
        var a1 = gl43[
            2, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 1, 0,
            0, 0, 0, 2
        ];
        var b1 = gl43[
            1, 0, 1, 0,
            1, 0, 0, 0,
            0, 1, 0, 1,
            0, 2, 0, 0
        ];

        var sp43 = Group.Generate("Sp4(3)", gl43, a1, b1);
        var z = Group.Zentrum(sp43);
        var PSp43 = sp43.Over(z, "PSp4(3)");
        DisplayGroup.HeadOrders(PSp43);
    }
}
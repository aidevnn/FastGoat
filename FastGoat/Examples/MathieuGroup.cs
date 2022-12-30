using FastGoat.Structures;
using FastGoat.UserGroup.Matrix;
using FastGoat.UserGroup.Perms;

namespace FastGoat.Examples;

public static class MathieuGroup
{
    public static void M11()
    {
        var s11 = new Sn(11);
        var a = s11[(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11)];
        var b = s11[(3, 7, 11, 8), (4, 10, 5, 6)];
        var m11 = Group.Generate("M11", s11, a, b);
        DisplayGroup.Head(m11);
    }

    public static void M12()
    {
        var s12 = new Sn(12);
        var a = s12[(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11)];
        var b = s12[(3, 7, 11, 8), (4, 10, 5, 6)];
        var c = s12[(1, 12), (2, 11), (3, 6), (4, 8), (5, 9), (7, 10)];
        var m12 = Group.Generate("M12", s12, a, b, c);
        DisplayGroup.Head(m12);
    }

    public static void M21()
    {
        var s21 = new Sn(21);
        // (1,4,5,9,3)(2,8,10,7,6)(12,15,16,20,14)(13,19,21,18,17), (1,21,5,12,20)(2,16,3,4,17)(6,18,7,19,15)(8,13,9,14,11)
        var a = s21[(1, 4, 5, 9, 3), (2, 8, 10, 7, 6), (12, 15, 16, 20, 14), (13, 19, 21, 18, 17)];
        var b = s21[(1, 21, 5, 12, 20), (2, 16, 3, 4, 17), (6, 18, 7, 19, 15), (8, 13, 9, 14, 11)];
        var m21 = Group.Generate("M21", s21, a, b);
        DisplayGroup.Head(m21);

        var gl = new GLnq(3, 4);

        var x = gl.Fq['x'];
        var a1 = gl[
            x, 0, 0,
            0, x * x, 0,
            0, 0, 1
        ];
        var b1 = gl[
            1, 0, 1,
            1, 0, 0,
            0, 1, 0
        ];

        var SL34 = Group.Generate($"SL(3,4)", gl, a1, b1);
        var z = Group.Zentrum(SL34);
        var L34 = SL34.Over(z, "L3(4)");
        DisplayGroup.Head(L34);

        DisplayGroup.AreIsomorphics(m21, L34);

        // |M21| = 20160
        // Type        NonAbelianGroup
        // BaseGroup   S21
        //
        // |L3(4)| = 20160
        // Type        NonAbelianGroup
        // BaseGroup   SL(3,4)/Z(SL(3,4))
        // Group           |SL(3,4)| = 60480
        // NormalSubGroup  |Z(SL(3,4))| = 3
        //
        // M21 IsIsomorphicTo L3(4) : True
    }
}
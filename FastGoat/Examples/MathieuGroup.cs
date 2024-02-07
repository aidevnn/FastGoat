using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.UserGroup.Matrix;
using FastGoat.UserGroup.Perms;
using FastGoat.UserGroup.Words.Tools;

namespace FastGoat.Examples;

public static class MathieuGroup
{
    public static void M11()
    {
        GlobalStopWatch.Restart();
        var s11 = new Sn(11);
        var a = s11[(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11)];
        var b = s11[(3, 7, 11, 8), (4, 10, 5, 6)];
        var m11 = Group.Generate("M11", s11, a, b);
        DisplayGroup.Head(m11);
        Graph.DefiningRelatorsOfGroup(m11);
    }
    
    /* |M11| = 7920
       All Relators
           b4
           a11
           ab2ab2ab-2
           a2ba-2b2ab-1a2b-1
           a2ba2b2a-1ba-1b-1a-1b-1
           aba-1ba-1b-1aba-1ba-1b-1
       #  Time:5.116s
     */

    public static void M12()
    {
        GlobalStopWatch.Restart();
        var s12 = new Sn(12);
        var a = s12[(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11)];
        var b = s12[(3, 7, 11, 8), (4, 10, 5, 6)];
        var c = s12[(1, 12), (2, 11), (3, 6), (4, 8), (5, 9), (7, 10)];
        var m12 = Group.Generate("M12", s12, a, b, c);
        DisplayGroup.Head(m12);
        Graph.DefiningRelatorsOfGroup(m12);
    }
    
    /* |M12| = 95040
       All Relators
           b4
           c2
           a11
           acacac
           bcbcbcbc
           ab2ab2ab-2
           b2cb2cbcb-1c
           abcb-1ca2cbcb-1
           cab-1cab-1cab-1
           a3cb-1cba-1b-1cbc
           ab-1a-1caba-1bcb-1
       #  Time:2m25s
     */
    
    public static void M21()
    {
        GlobalStopWatch.Restart();
        var s21 = new Sn(21);
        var a = s21[(1, 4, 5, 9, 3), (2, 8, 10, 7, 6), (12, 15, 16, 20, 14), (13, 19, 21, 18, 17)];
        var b = s21[(1, 21, 5, 12, 20), (2, 16, 3, 4, 17), (6, 18, 7, 19, 15), (8, 13, 9, 14, 11)];
        var m21 = Group.Generate("M21", s21, a, b);
        DisplayGroup.Head(m21);
        Graph.DefiningRelatorsOfGroup(m21);
    }
    
    /* |M21| = 20160
       All Relators
           a5
           b5
           baba-1baba-1baba-1
           abab-1a-1b2a-1ba-1b
           a2b-1a2b-1a2b-1a2b-1
           ababa-1b-1ababa-1b-1
           ba-1ba-1ba-1ba-1ba-1ba-1ba-1
       #  Time:13.611s
     */

    public static void M21Iso()
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
    
    public static void M22()
    {
        var s22 = new Sn(22);
        
        // var a1 = s22[(4, 15),(6, 18),(7, 8),(9, 13),(10, 16),(11, 21),(12, 19),(20, 22)];
        // var a2 = s22[(4, 10),(6, 11),(7, 12),(8, 19),(9, 22),(13, 20),(15, 16),(18, 21)];
        // var a3 = s22[(4, 20),(6, 8),(7, 18),(9, 16),(10, 13),(11, 19),(12, 21),(15, 22)];
        // var a4 = s22[(4, 8),(6, 20),(7, 15),(9, 21),(10, 19),(11, 13),(12, 16),(18, 22)];
        // var b = s22[(3, 4),(5, 22),(6, 18),(7, 12),(8, 21),(11, 19),(13, 14),(16, 17)];
        // var c = s22[(2, 3),(6, 18),(9, 21),(10, 12),(11, 13),(14, 17),(16, 19),(20, 22)];
        // var d = s22[(1, 2),(6, 15),(7, 10),(8, 11),(9, 19),(14, 17),(16, 22),(20, 21)];
        
        var x = s22[(1, 4, 5, 9, 3), (2, 8, 10, 7, 6), (12, 15, 16, 20, 14), (13, 19, 21, 18, 17)];
        var y = s22[(1, 21), (2, 10, 8, 6), (3, 13, 4, 17), (5, 19, 9, 18), (11, 22), (12, 14, 16, 20)];
        var z = s22[(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11), (12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22)];

        GlobalStopWatch.Restart();
        var m22 = Group.Generate("M22", s22, x, y, z);
        DisplayGroup.Head(m22);
        GlobalStopWatch.Show("M22"); // # M22 Time:2.884s
        Graph.DefiningRelatorsOfGroup(m22);
    }
    /* |M22| = 443520
        ...
        ...
       All Relators
           b5
           c4
           b2cbc-1
           a3ba-1b-1
           a4b-1a-1b
           acacacacacac
           a2ca2ca2ca2ca2c
           ac2a-1c2ac-2a-1c-2
           acac-1acac-1acac-1
           c2a2c-1bcac-1b-1c-1
           abc2acacb-1a-1bc-1ac-1
           acb-1cac-1bc-1ac-1bc-1
       #  Time:24m4s
     */
}
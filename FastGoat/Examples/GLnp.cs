using FastGoat.Structures;
using FastGoat.UserGroup.Matrix;
using FastGoat.UserGroup.Perms;
using FastGoat.UserGroup.Words;

namespace FastGoat.Examples;

public static class GLnp
{
    public static void Minimal()
    {
        {
            var gl = new GL(2, 3);
            var a = gl[2, 0, 0, 1];
            var b = gl[2, 1, 2, 0];

            var g1 = Group.Generate(gl, a, b);
            DisplayGroup.HeadElements(g1); // 48 elements
        }

        {
            var gl = new GL(2, 3);
            var a = gl[1, 1, 0, 1];
            var b = gl[0, 1, 2, 0];

            var g1 = Group.Generate(gl, a, b);
            DisplayGroup.Head(g1); // 24 elements
        }

        {
            var gl = new GL(2, 5);
            var a = gl[2, 0, 0, 1];
            var b = gl[4, 1, 4, 0];
            var g1 = Group.Generate(gl, a, b);
            DisplayGroup.Head(g1); // 480 elements
        }

        {
            var gl = new GL(2, 7);
            var a = gl[3, 0, 0, 1];
            var b = gl[6, 1, 6, 0];
            var g1 = Group.Generate(gl, a, b);
            DisplayGroup.Head(g1); // 2016 elements
        }

        {
            var gl = new GL(2, 11);
            var a = gl[2, 0, 0, 1];
            var b = gl[10, 1, 10, 0];
            var g1 = Group.Generate(gl, a, b);
            DisplayGroup.Head(g1); // 13200 elements
        }

        {
            var gl = new GL(3, 2);
            var a = gl.At((0, 1, 4, 8), (1, 1, 1, 1));
            var b = gl.At((2, 3, 7), (1, 1, 1));
            var g1 = Group.Generate(gl, a, b);
            DisplayGroup.Head(g1); // 168 elements
        }

        {
            var gl = new GL(3, 3);
            var a = gl.At((0, 4, 8), (2, 1, 1));
            var b = gl.At((0, 2, 3, 7), (2, 1, 2, 2));
            var g1 = Group.Generate(gl, a, b);
            DisplayGroup.Head(g1); // 11232 elements
        }

        {
            var gl = new GL(4, 2);
            var a = gl.At((0, 1, 5, 10, 15), (1, 1, 1, 1, 1));
            var b = gl.At((3, 4, 9, 14), (1, 1, 1, 1));
            var g1 = Group.Generate(gl, a, b);
            DisplayGroup.Head(g1); // 20160 elements
        }
    }

    public static void SmallGroup3244()
    {
        // GroupName https://people.maths.bris.ac.uk/~matyd/GroupNames/1/C8.C2%5E2.html
        // G:=sub<GL(4,GF(3))| [0,1,0,0,2,0,2,0,0,1,0,1,0,0,2,0],[0,0,0,2,2,0,1,0,0,1,0,2,2,0,0,0],[2,0,0,0,0,1,0,0,0,0,2,0,0,0,0,1] >;
        // G:=Group( (1,2,3,4,5,6,7,8)(9,10,11,12,13,14,15,16), (2,4)(3,7)(6,8)(10,12)(11,15)(14,16), (1,11)(2,16)(3,13)(4,10)(5,15)(6,12)(7,9)(8,14) );
        var gl = new GL(4, 3);
        var m1 = gl[0, 1, 0, 0, 2, 0, 2, 0, 0, 1, 0, 1, 0, 0, 2, 0];
        var m2 = gl[0, 0, 0, 2, 2, 0, 1, 0, 0, 1, 0, 2, 2, 0, 0, 0];
        var m3 = gl[2, 0, 0, 0, 0, 1, 0, 0, 0, 0, 2, 0, 0, 0, 0, 1];
        var g1 = Group.Generate("gl(C8 . (C2 x C2))", gl, m1, m2, m3);
        DisplayGroup.Head(g1);

        var s16 = new Sn(16);
        var p1 = s16[(1, 2, 3, 4, 5, 6, 7, 8), (9, 10, 11, 12, 13, 14, 15, 16)];
        var p2 = s16[(2, 4), (3, 7), (6, 8), (10, 12), (11, 15), (14, 16)];
        var p3 = s16[(1, 11), (2, 16), (3, 13), (4, 10), (5, 15), (6, 12), (7, 9), (8, 14)];

        var g2 = Group.Generate("pg(C8 . (C2 x C2))", s16, p1, p2, p3);
        DisplayGroup.Head(g2);

        var g3 = new WordGroup("wg(C8 . (C2 x C2))", "a8, b2, c2, bab=a3, cac=a5, cbc=a4b");
        DisplayGroup.Head(g3);

        Console.WriteLine("({0}) IsIsomorphicTo ({1}) : {2}", g1, g2, g1.IsIsomorphicTo(g2));
        Console.WriteLine("({0}) IsIsomorphicTo ({1}) : {2}", g1, g3, g1.IsIsomorphicTo(g3));
    }
}
using FastGoat.UserGroup;

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
}
using System;
using System.Linq;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.UserGroup.Matrix;
using Xunit;

namespace Tests;

public class GLnpUnitTest
{
    [Fact]
    public void Test1Ops()
    {
        Random rnd = new Random();
        for (int n = 2; n < 6; n++) // n = 2, 3, 4, 5
        {
            for (int p0 = 0; p0 < 5; p0++) // p = 2, 3, 5, 7, 11
            {
                int p = IntExt.Primes10000[p0];

                var gl = new GL(n, p);
                var rg = Enumerable.Range(0, n * n).ToArray();
                for (int j = 0; j < 5; j++)
                {
                    var mat0 = rg.Select(i => rnd.Next(p)).ToArray();
                    var mat1 = rg.Select(i => rnd.Next(p)).ToArray();
                    var det0 = IntExt.AmodP(MatrixExt.ComputeDeterminant(mat0), p);
                    var det1 = IntExt.AmodP(MatrixExt.ComputeDeterminant(mat1), p);
                    if (det0 == 0 || det1 == 0)
                    {
                        --j;
                        continue;
                    }

                    var inv0 = MatrixExt.InvertModP(p, mat0);
                    var inv1 = MatrixExt.InvertModP(p, mat1);
                    var dot = MatrixExt.DotModP(p, mat0, mat1);

                    var e0 = gl.Create(mat0);
                    var ei0 = gl.Invert(e0);
                    var f0 = gl.Create(inv0);
                    var e1 = gl.Create(mat1);
                    var ei1 = gl.Invert(e1);
                    var f1 = gl.Create(inv1);
                    var e2 = gl.Op(e0, e1);
                    var f2 = gl.Create(dot);
                    Assert.Equal(det0, gl.Det(e0));
                    Assert.Equal(det1, gl.Det(e1));
                    Assert.True(inv0.SequenceEqual(ei0.Table));
                    Assert.True(inv1.SequenceEqual(ei1.Table));
                    Assert.True(dot.SequenceEqual(e2.Table));
                    Assert.Equal(f0, ei0);
                    Assert.Equal(f1, ei1);
                    Assert.Equal(f2, e2);
                }
            }
        }
    }

    [Fact]
    public void Test2GL23()
    {
        var gl = new GL(2, 3);
        var a1 = gl[2, 0, 0, 1];
        var b1 = gl[2, 1, 2, 0];
        var g1 = Group.Generate(gl, a1, b1);
        Assert.Equal(48, g1.Count());
        Assert.Equal(24, g1.Count(e => IntExt.AmodP(gl.Det(e), 3) == 1));
        Assert.Equal(24, g1.Count(e => IntExt.AmodP(gl.Det(e), 3) == 2));

        var a2 = gl[1, 1, 0, 1];
        var b2 = gl[0, 1, 2, 0];
        var g2 = Group.Generate(gl, a2, b2);
        Assert.Equal(24, g2.Count());
        Assert.Equal(24, g2.Count(e => gl.Det(e) == 1));
        Assert.True(g1.ToHashSet().IsSupersetOf(g2));
    }

    [Fact]
    public void Test3GL3p()
    {
        var gl32 = new GL(3, 2);
        var a1 = gl32.At((0, 1, 4, 8), (1, 1, 1, 1));
        var b1 = gl32.At((2, 3, 7), (1, 1, 1));
        var g1 = Group.Generate(gl32, a1, b1);
        Assert.Equal(168, g1.Count());

        var gl33 = new GL(3, 3);
        var a2 = gl33.At((0, 4, 8), (2, 1, 1));
        var b2 = gl33.At((0, 2, 3, 7), (2, 1, 2, 2));
        var g2 = Group.Generate(gl33, a2, b2);
        Assert.Equal(11232, g2.Count());
    }
}
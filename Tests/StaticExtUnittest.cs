using System;
using System.Linq;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Polynoms;
using Xunit;

namespace Tests;

public class StaticExtUnittest
{
    [Fact]
    public void Test1Matrix()
    {
        Random rnd = new Random();
        for (int n = 2; n < 6; n++)
        {
            var p = rnd.Next(2, 10);
            var diag0 = MatrixExt.Diagonal(n, p);
            var diag1 = Enumerable.Range(0, n * n).Select(i => i % (n + 1) == 0 ? p : 0).ToArray();
            Assert.True(diag0.SequenceEqual(diag1));

            var rg = Enumerable.Range(0, n * n).ToArray();
            for (int k = 0; k < 5; k++)
            {
                var mat = rg.Select(i => rnd.Next(p)).ToArray();
                var det = MatrixExt.ComputeDeterminant(mat);
                var com = MatrixExt.Comatrix(mat);
                var tcom = MatrixExt.Transpose(com);
                var dot = MatrixExt.Dot(tcom, mat);
                Assert.True(dot.SequenceEqual(MatrixExt.Diagonal(n, det)));
            }
        }
    }
    
    [Fact]
    public void Test2BezoutGcd()
    {
        var zn = new Zn(3 * 5 * 17);
        var seq = Group.Generate(zn, zn[3]).Grid2D(Group.Generate(zn, zn[5]));
        bool TestBz(int a, int b)
        {
            var d = IntExt.Gcd(a, b);
            var (x, y) = IntExt.Bezout(a, b);
            return d >= 0 && d == a * x + b * y;
        }

        Assert.True(seq.All(e => TestBz(e.t1.K, e.t2.K)));
    }

    [Fact]
    public void Test3FpPolynom()
    {
        var x = new FpPolynom(5);
        var a0 = (x + 1) * (x + 2);
        var a1 = (x + 1) * (x + 3);
        var a2 = (x + 2) * (x + 3);
        var a3 = (x + 1) * (x + 2) * (x + 3);

        var q0 = a3 / (x + 1);
        var q1 = a3 / (x + 2);
        var q2 = a3 / (x + 3);
        Assert.Equal(q0, a2);
        Assert.Equal(q1, a1);
        Assert.Equal(q2, a0);

        var x4 = x + 4;
        var qr = a3.DivideBy(x4);
        Assert.Equal(a3, qr.quo * x4 + qr.rem);
    }

    [Fact]
    public void Test4FpBezout()
    {
        var x = new FpPolynom(7);
        var a = x * x + 3 * x;
        var b = x.Pow(4) + x * x;
        var gcd = FpPolynom.Gcd(a, b);
        var (x0, y0) = FpPolynom.Bezout(a, b);
        Assert.Equal(gcd, a * x0 + b * y0);
    }
}
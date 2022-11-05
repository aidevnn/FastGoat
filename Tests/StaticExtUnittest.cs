using System;
using System.Linq;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.UserGroup.Integers;
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

}
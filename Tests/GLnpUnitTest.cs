using System;
using System.Linq;
using FastGoat;
using FastGoat.UserGroup;
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
                    var e1 = gl.Create(mat1);
                    var ei1 = gl.Invert(e1);
                    var e2 = gl.Op(e0, e1);
                    Assert.Equal(det0, gl.Det(e0));
                    Assert.Equal(det1, gl.Det(e1));
                    Assert.True(inv0.SequenceEqual(ei0.Table));
                    Assert.True(inv1.SequenceEqual(ei1.Table));
                    Assert.True(dot.SequenceEqual(e2.Table));
                }
            }
        }
    }
}
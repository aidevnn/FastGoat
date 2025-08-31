using System;
using System.Linq;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.UserGroup.Matrix;
using Xunit;

namespace Tests;

public class DPGLUnitTest
{
    [Fact]
    public void Test1Constructors()
    {
        var (n, p) = (4, 17);
        var dpgl = new DPGL(n, p);

        Assert.True(dpgl.Neutral().SequenceEqual([1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]));
        Assert.True(dpgl[1, 3, 2, 5].SequenceEqual([1, 0, 0, 0, 0, 3, 0, 0, 0, 0, 2, 0, 0, 0, 0, 5]));

        int[] d1 = [7, 3, 11, 15];
        var p1 = dpgl.Sn[(4, 3, 2, 1)];
        var e1 = dpgl[d1, p1];
        var e2 = dpgl[[7, 3, 11, 15], (4, 3, 2, 1)];
        Assert.Equal(e1, e2);
        Assert.True(e1.SequenceEqual([0, 0, 0, 7, 3, 0, 0, 0, 0, 11, 0, 0, 0, 0, 15, 0]));
    }

    [Fact]
    public void Test2Invert()
    {
        var (n, p) = (4, 17);
        var dpgl = new DPGL(n, p);
        var gl = dpgl.GL;
        
        foreach (var permSn in dpgl.Sn)
        {
            var diag = n.SeqLazy().Select(_ => IntExt.Rng.Next(1, p)).ToArray();
            
            var e = dpgl[diag, permSn];
            var ei = dpgl.Invert(e);

            var matPerm = gl.Create(MatrixExt.Permutation(permSn.Table));
            var matDiag = gl.Create(MatrixExt.Diagonal(diag));
            var m = gl.Op(matDiag, matPerm);
            var mi = gl.Invert(m);
            Assert.Equal(e.ToGL(), m);
            Assert.Equal(ei.ToGL(), mi);
        }
    }

    [Fact]
    public void Test3Multiplication()
    {
        var (n, p) = (4, 17);
        var dpgl = new DPGL(n, p);
        var gl = dpgl.GL;
        var allPerms = dpgl.Sn.ToArray();

        for (int i = 0; i < 10; i++)
        {
            var d1 = n.SeqLazy().Select(_ => IntExt.Rng.Next(1, p)).ToArray();
            var p1 = DistributionExt.Dice(allPerms);
            var e1 = dpgl[d1, p1];
            var m1 = gl.Op(gl.Create(MatrixExt.Diagonal(d1)), gl.Create(MatrixExt.Permutation(p1.Table)));

            var d2 = n.SeqLazy().Select(_ => IntExt.Rng.Next(1, p)).ToArray();
            var p2 = DistributionExt.Dice(allPerms);
            var e2 = dpgl[d2, p2];
            var m2 = gl.Op(gl.Create(MatrixExt.Diagonal(d2)), gl.Create(MatrixExt.Permutation(p2.Table)));
        
            Assert.Equal(gl.Op(m1, m2), dpgl.Op(e1, e2).ToGL());
        }
    }
}
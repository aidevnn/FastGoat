using System;
using FastGoat;
using FastGoat.Gp;
using FastGoat.UserGroup;
using Xunit;

namespace Tests;

public class SemiDirectProdUnitTest
{
    [Fact]
    public void Test1ThrowException()
    {
        var z3 = new Zn(3);
        var z2Xz2 = Product.Group(new Zn(2), new Zn(2));
        var v = Group.Generate(z2Xz2[1, 0], z2Xz2[0, 1]);
        var g3 = Group.Generate(z3[1]);

        Assert.Throws<GroupException>(() => Group.SemiDirectProd(v, g3));
        Assert.Throws<GroupException>(() => Group.SemiDirectProd(g3, g3));
    }

    [Fact]
    public void Test2SemiDirectProduct()
    {
        var g2 = Group.Generate(new Zn(2)[1]);
        var g3 = Group.Generate(new Zn(3)[1]);
        var g4 = Group.Generate(new Zn(4)[1]);
        var d3 = Group.SemiDirectProd(g3, g2);
        var d4 = Group.SemiDirectProd("D4", g4, g2);

        var s3 = new Sn(3);
        var g6 = Group.Generate(s3[(1, 2)], s3[(1, 2, 3)]);
        Assert.True(d3.IsIsomorphicTo(g6));
        Assert.Equal("D4", d4.Name);
    }
}
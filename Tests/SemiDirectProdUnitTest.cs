using System;
using System.Linq;
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
        var g3 = Group.Generate(z3[1]);
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

    [Fact]
    public void Test3ComplexSemiDirectProduct()
    {
        var g1 = Group.SemiDirectProd(Product.Generate(new Cn(3), new Cn(3)), new Cn(2));
        Assert.Equal(GroupType.NonAbelianGroup, g1.GroupType);
        Assert.Equal(18, g1.Count());
        Assert.True(Group.IsGroup(g1));

        var g2 = Group.SemiDirectProd(new Cn(8), Product.Generate(new Cn(2), new Cn(2)));
        Assert.Equal(GroupType.NonAbelianGroup, g2.GroupType);
        Assert.Equal(32, g2.Count());
        Assert.True(Group.IsGroup(g2));

        var g3 = Group.SemiDirectProd(Product.Generate(new Cn(2), new Cn(2), new Cn(2)),
            Product.Generate(new Cn(2), new Cn(2)));
        Assert.Equal(GroupType.NonAbelianGroup, g3.GroupType);
        Assert.Equal(32, g3.Count());
        Assert.True(Group.IsGroup(g3));

        var g4 = Group.SemiDirectProd(new Cn(3), new Cn(8));
        var g5 = Group.SemiDirectProd(g4, new Cn(2));
        Assert.Equal(GroupType.NonAbelianGroup, g5.GroupType);
        Assert.Equal(48, g5.Count());
        Assert.True(Group.IsGroup(g5));
    }
}
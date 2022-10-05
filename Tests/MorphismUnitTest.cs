using System.Collections.Generic;
using System.Linq;
using FastGoat;
using FastGoat.Gp;
using FastGoat.UserGroup;
using Xunit;

namespace Tests;

public class MorphismUnitTest
{
    [Fact]
    public void Test1NotMorphism()
    {
        var sn = new Sn(3);
        var a = sn[(1, 2)];
        var b = sn[(1, 2, 3)];
        var g = Group.Generate(a, b);
        var pMap = Group.PartialMap((a, sn.Neutral()), (b, sn[(1, 2)]));
        var hMap = Group.EndomorphismMap(g, pMap);
        Assert.Empty(hMap);
    }

    [Fact]
    public void Test2Alternate3()
    {
        var sn = new Sn(3);
        var a = sn[(1, 2, 3)];
        var b = sn[(1, 2)];
        var g = Group.Generate(a, b);
        var pMap = Group.PartialMap((a, sn.Neutral()), (b, sn[(1, 2)]));
        var hMap = Group.EndomorphismMap(g, pMap);
        Assert.Equal(6, hMap.Count);
        Assert.Equal(2, hMap.Values.Distinct().Count());
    }

    [Fact]
    public void Test3Alternate4()
    {
        var sn = new Sn(4);
        var a = sn[(1, 2, 3)];
        var b = sn[(1, 2, 3, 4)];
        var g = Group.Generate(a, b);
        var pMap = Group.PartialMap((a, sn.Neutral()), (b, sn[(1, 2)]));
        var hMap = Group.EndomorphismMap(g, pMap);
        Assert.Equal(24, hMap.Count);
        Assert.Equal(2, hMap.Values.Distinct().Count());
    }

    [Fact]
    public void Test4Alternate5()
    {
        var sn = new Sn(5);
        var a = sn[(1, 2, 3, 4, 5)];
        var b = sn[(1, 2, 3, 4)];
        var g = Group.Generate(a, b);
        var pMap = Group.PartialMap((a, sn.Neutral()), (b, sn[(1, 2)]));
        var hMap = Group.EndomorphismMap(g, pMap);
        Assert.Equal(120, hMap.Count);
        Assert.Equal(2, hMap.Values.Distinct().Count());
    }

    [Fact]
    public void Test5NotMorphism()
    {
        var c20 = new Cn(20);
        var pMap = Group.PartialMap((c20[2], c20[3]));
        var hMap = Group.EndomorphismMap(c20, pMap);
        Assert.Empty(hMap);
    }

    [Fact]
    public void Test6C20()
    {
        var c20 = new Cn(20);
        var pMap = Group.PartialMap((c20[1], c20[2]));
        var hMap = Group.EndomorphismMap(c20, pMap);
        Assert.Equal(20, hMap.Count);
        Assert.Equal(10, hMap.Values.Distinct().Count());
    }

    [Fact]
    public void Test7G600()
    {
        var g = Group.Create(Product.Group(new Cn(20), new Cn(30)));
        var pMap = Group.PartialMap((g[1, 1], g[0, 0]), (g[3, 4], g[2, 3]));
        var hMap = Group.EndomorphismMap(g, pMap);
        Assert.Equal(600, hMap.Count);
        Assert.Equal(10, hMap.Values.Distinct().Count());
    }

    [Fact]
    public void Test7G4320()
    {
        var g = Group.Create(Product.Group(new Cn(8), new Cn(18), new Cn(30)));
        var h = Group.Create(Product.Group(new Cn(2), new Cn(6)));
        var pMap = Group.PartialMap((g[1, 1, 1], h[0, 0]), (g[5, 1, 6], h[0, 1]), (g[1, 10, 1], h[1, 0]));
        var hMap = Group.HomomorphismMap(g, h, pMap);
        Assert.Equal(4320, hMap.Count);
        Assert.Equal(12, hMap.Values.Distinct().Count());
    }
}
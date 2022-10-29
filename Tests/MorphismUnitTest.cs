using System.Collections.Generic;
using System.Linq;
using FastGoat;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Perms;
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
        var g = Group.Generate(sn, a, b);
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
        var g = Group.Generate(sn, a, b);
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
        var g = Group.Generate(sn, a, b);
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
        var g = Group.Generate(sn, a, b);
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

    [Fact]
    public void Test8Isomorphism()
    {
        var s4 = new Sn(4);
        var g1 = Group.SemiDirectProd(new Cn(4), new Cn(2));
        var h1 = Group.Generate(s4, s4[(1, 2), (3, 4)], s4[(1, 2, 3, 4)]);
        var pMap1 = Group.PartialMap((g1[1, 0], h1[(1, 2, 3, 4)]), (g1[0, 1], h1[(1, 2), (3, 4)]));
        var hMap1 = Group.IsomorphismMap(g1, h1, pMap1);
        Assert.Equal(8, hMap1.Count);
        Assert.Equal(8, hMap1.Values.Distinct().Count());

        var s5 = new Sn(5);
        var c = s5[(1, 2), (3, 4, 5)];
        var g2 = Group.Generate(s5, c);
        var h2 = Group.Create(Product.Group(new Cn(2), new Cn(3)));
        var pMap2 = Group.PartialMap((g2[(1, 2)], h2[1, 0]), (g2[(3, 4, 5)], h2[0, 1]));
        var hMap2 = Group.IsomorphismMap(g2, h2, pMap2);
        Assert.Equal(6, hMap2.Count);
        Assert.Equal(6, hMap2.Values.Distinct().Count());
    }

    [Fact]
    public void Test9Automorphism()
    {
        var s7 = new Sn(7);
        var a = s7[(1, 2, 3, 4, 5, 6, 7)];
        var b = s7[(1, 3, 5, 7, 2, 4, 6)];
        var g = Group.Generate(s7, a);
        var pMap1 = Group.PartialMap((a, b));
        var hMap1 = Group.AutomorphismMap(g, pMap1);
        Assert.Equal(7, hMap1.Count);
        Assert.Equal(7, hMap1.Values.Distinct().Count());

        var h = Group.Generate(s7, s7[(1, 2, 3)], s7[(4, 5, 6)]);
        var pMap2 = Group.PartialMap((s7[(1, 2, 3)], s7[(4, 5, 6)]), (s7[(4, 6, 5)], s7[(1, 2, 3)]));
        var hMap2 = Group.AutomorphismMap(h, pMap2);
        Assert.Equal(9, hMap2.Count);
        Assert.Equal(9, hMap2.Values.Distinct().Count());
    }
}
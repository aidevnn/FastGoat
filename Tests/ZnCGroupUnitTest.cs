using System.Linq;
using Xunit;
using FastGoat;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Integers;

namespace Tests;

public class ZnCGroupUnitTest
{
    [Fact]
    public void Test1GenerateElement()
    {
        var zn = new Zn(30);
        var g1 = Group.GenerateElements(zn, zn[6]);
        var g2 = Group.GenerateElements(zn, zn[10]);
        var g3 = Group.GenerateElements(zn, zn[6], zn[10]);
        Assert.Equal(5, g1.Count());
        Assert.Equal(3, g2.Count());
        Assert.Equal(15, g3.Count());
    }

    [Fact]
    public void Test2Cycle()
    {
        var zn = new Zn(8);
        var cycle = Group.Cycle(zn, zn[2]);
        Assert.Equal(4, cycle.Count);
        Assert.Equal(4, cycle[zn[0]]);
        Assert.Equal(1, cycle[zn[2]]);
        Assert.Equal(2, cycle[zn[4]]);
        Assert.Equal(3, cycle[zn[6]]);
    }

    [Fact]
    public void Test3LongestCycle()
    {
        var g = Product.Group(new Zn(2), new Zn(4));
        var elements = Group.GenerateElements(g, g[1, 0], g[0, 1]);
        var longCycles = Group.LongestCycles(g, elements);
        Assert.Equal(4, longCycles.Count);

        Assert.Equal(1, longCycles[g[1, 0]][g[1, 0]]);
        Assert.Equal(2, longCycles[g[1, 0]][g[0, 0]]);

        Assert.Equal(1, longCycles[g[1, 2]][g[1, 2]]);
        Assert.Equal(2, longCycles[g[1, 2]][g[0, 0]]);

        Assert.Equal(1, longCycles[g[0, 1]][g[0, 1]]);
        Assert.Equal(2, longCycles[g[0, 1]][g[0, 2]]);
        Assert.Equal(3, longCycles[g[0, 1]][g[0, 3]]);
        Assert.Equal(4, longCycles[g[0, 1]][g[0, 0]]);

        Assert.Equal(1, longCycles[g[1, 1]][g[1, 1]]);
        Assert.Equal(2, longCycles[g[1, 1]][g[0, 2]]);
        Assert.Equal(3, longCycles[g[1, 1]][g[1, 3]]);
        Assert.Equal(4, longCycles[g[1, 1]][g[0, 0]]);
    }

    [Fact]
    public void Test4ConcreteGroup()
    {
        var zn = new Zn(40);
        var g0 = Group.Generate(zn, zn[1]);
        var g1 = Group.Generate(g0, zn[10]);
        var g2 = Group.Generate(g0, zn[8]);
        var g3 = Group.Create(Product.Group(g1, g2));
        var g4 = Group.Generate(g0, zn[2]);
        var g5 = Group.DirectProduct(g1, g2);
        Assert.Equal(40, g0.Count());
        Assert.Equal(4, g1.Count());
        Assert.Equal(5, g2.Count());
        Assert.Equal(20, g3.Count());
        Assert.Equal(20, g4.Count());
        Assert.Equal(20, g5.Count());
        Assert.Equal(GroupType.AbelianGroup, g0.GroupType);
        Assert.Equal(GroupType.AbelianGroup, g1.GroupType);
        Assert.True(g3.IsIsomorphicTo(g4));
        Assert.True(g3.IsIsomorphicTo(g5));
    }
}
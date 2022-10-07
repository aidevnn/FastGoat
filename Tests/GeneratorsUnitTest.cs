using System.Linq;
using FastGoat;
using FastGoat.Gp;
using FastGoat.UserGroup;
using Xunit;

namespace Tests;

public class GeneratorsUnitTest
{
    [Fact]
    public void Test1ZnGenerators()
    {
        var g0 = new Zn(5);
        var g0G = g0.GetGenerators().ToArray();
        var gr0 = Group.Generate(g0, g0G);
        Assert.Contains(g0[1], g0G);
        Assert.Single(g0G);
        Assert.Equal(5, gr0.Count());

        var g1 = Product.Group(g0, g0);
        var g1G = g1.GetGenerators().ToArray();
        var gr1 = Group.Generate(g1, g1G);
        Assert.Contains(g1[1, 0], g1G);
        Assert.Contains(g1[0, 1], g1G);
        Assert.Equal(2, g1G.Length);
        Assert.Equal(25, gr1.Count());

        var g2 = Product.Group(g0, g0, g0);
        var g2G = g2.GetGenerators().ToArray();
        var gr2 = Group.Generate(g2, g2G);
        Assert.Contains(g2[1, 0, 0], g2G);
        Assert.Contains(g2[0, 1, 0], g2G);
        Assert.Contains(g2[0, 0, 1], g2G);
        Assert.Equal(3, g2G.Length);
        Assert.Equal(125, gr2.Count());

        var g3 = Product.Group(g1, g0);
        var g3G = g3.GetGenerators().ToArray();
        var gr3 = Group.Generate(g3, g3G);
        Assert.Contains(g3[(1, 0), 0], g3G);
        Assert.Contains(g3[(0, 1), 0], g3G);
        Assert.Contains(g3[(0, 0), 1], g3G);
        Assert.Equal(3, g3G.Length);
        Assert.Equal(125, gr3.Count());
    }

    [Fact]
    public void Test2SnGenerators()
    {
        var s3 = new Sn(3);
        var s3G = s3.GetGenerators().ToArray();
        var gr3 = Group.Generate(s3, s3G);
        Assert.Contains(s3[(1, 2)], s3G);
        Assert.Contains(s3[(2, 3)], s3G);
        Assert.Equal(2, s3G.Length);
        Assert.Equal(6, gr3.Count());

        var s4 = new Sn(4);
        var s4G = s4.GetGenerators().ToArray();
        var gr4 = Group.Generate(s4, s4G);
        Assert.Contains(s4[(1, 2)], s4G);
        Assert.Contains(s4[(2, 3)], s4G);
        Assert.Contains(s4[(3, 4)], s4G);
        Assert.Equal(3, s4G.Length);
        Assert.Equal(24, gr4.Count());

        var s5 = new Sn(5);
        var s5G = s5.GetGenerators().ToArray();
        var gr5 = Group.Generate(s5, s5G);
        Assert.Contains(s5[(1, 2)], s5G);
        Assert.Contains(s5[(2, 3)], s5G);
        Assert.Contains(s5[(3, 4)], s5G);
        Assert.Contains(s5[(5, 4)], s5G);
        Assert.Equal(4, s5G.Length);
        Assert.Equal(120, gr5.Count());
    }

    [Fact]
    public void Test3Monogenics()
    {
        var c0 = new Cn(40);
        var c1 = Group.Generate(c0, c0[4]);
        var c1G = c1.GetGenerators().ToArray();
        Assert.Single(c1G);
        Assert.Contains(c0[4], c1G);

        var c2 = Group.Generate(c0, c0[8]);
        var g1 = Product.Group(c1, c2);
        var g1G = g1.GetGenerators().ToArray();
        var gr1 = Group.Generate(g1, g1G);
        Assert.Equal(2, g1G.Length);
        Assert.Equal(50, gr1.Count());

        var g2 = Product.Generate(c1, c2);
        var act = () => g2.GetGenerators().ToArray();
        Assert.Throws<GroupException>(act);

        var gr2 = Group.Generate(g2, g2[4, 8]);
        var gr2G = gr2.GetGenerators().ToArray();
        Assert.Single(gr2G);
        Assert.Contains(g2[4, 8], gr2G);
        Assert.Equal(10, gr2.Count());
    }

    [Fact]
    public void Test4GapAllAutomorphisms()
    {
        // gap> Size(AllAutomorphisms(AbelianGroup([2,3]))); = 2
        var n1 = Product.Group(new Zn(2), new Zn(3));
        var autN1 = Group.Aut(n1, n1.GetGenerators().ToArray());
        Assert.Equal(2, autN1.Count());

        // gap> Size(AllAutomorphisms(AbelianGroup([3,3]))); = 48
        var n2 = Product.Group(new Zn(3), new Zn(3));
        var autN2 = Group.Aut(n2, n2.GetGenerators().ToArray());
        Assert.Equal(48, autN2.Count());

        // gap> Size(AllAutomorphisms(AbelianGroup([2,2,2]))); = 168
        var n3 = Product.Group(new Zn(2), new Zn(2), new Zn(2));
        var autN3 = Group.Aut(n3, n3.GetGenerators().ToArray());
        Assert.Equal(168, autN3.Count());

        // gap> Size(AllAutomorphisms(AbelianGroup([4,4]))); = 96
        var n4 = Product.Group(new Zn(4), new Zn(4));
        var autN4 = Group.Aut(n4, n4.GetGenerators().ToArray());
        Assert.Equal(96, autN4.Count());
    }
}
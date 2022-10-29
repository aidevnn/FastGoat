using System.Linq;
using FastGoat;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Perms;
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

        var g4 = Product.Gp(g0, g0);
        var g4G = g4.GetGenerators().ToArray();
        var gr4 = Group.Generate(g4, g4G);
        Assert.Contains(g4[1, 0], g4G);
        Assert.Contains(g4[0, 1], g4G);
        Assert.Equal(2, g4G.Length);
        Assert.Equal(25, gr4.Count());
    }

    [Fact]
    public void Test2SnGenerators()
    {
        var s3 = new Sn(3);
        var s3G = s3.GetGenerators().ToArray();
        var gr3 = Group.Generate(s3, s3G);
        Assert.Contains(s3[(1, 2)], s3G);
        Assert.Contains(s3[(1, 2, 3)], s3G);
        Assert.Equal(2, s3G.Length);
        Assert.Equal(6, gr3.Count());

        var s4 = new Sn(4);
        var s4G = s4.GetGenerators().ToArray();
        var gr4 = Group.Generate(s4, s4G);
        Assert.Contains(s4[(1, 2)], s4G);
        Assert.Contains(s4[(1, 2, 3, 4)], s4G);
        Assert.Equal(2, s4G.Length);
        Assert.Equal(24, gr4.Count());

        var s5 = new Sn(5);
        var s5G = s5.GetGenerators().ToArray();
        var gr5 = Group.Generate(s5, s5G);
        Assert.Contains(s5[(1, 2)], s5G);
        Assert.Contains(s5[(1, 2, 3, 4, 5)], s5G);
        Assert.Equal(2, s5G.Length);
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
        var n1 = Product.Generate(new Cn(2), new Cn(3));
        var autN1 = Group.AllAutomorphisms(n1);
        Assert.Equal(2, autN1.Count());

        // gap> Size(AllAutomorphisms(AbelianGroup([3,3]))); = 48
        var n2 = Product.Generate(new Cn(3), new Cn(3));
        var autN2 = Group.AllAutomorphisms(n2);
        Assert.Equal(48, autN2.Count());

        // gap> Size(AllAutomorphisms(AbelianGroup([2,2,2]))); = 168
        var n3 = Product.Generate(new Cn(2), new Cn(2), new Cn(2));
        var autN3 = Group.AllAutomorphisms(n3);
        Assert.Equal(168, autN3.Count());

        // gap> Size(AllAutomorphisms(AbelianGroup([4,4]))); = 96
        var n4 = Product.Generate(new Cn(4), new Cn(4));
        var autN4 = Group.AllAutomorphisms(n4);
        Assert.Equal(96, autN4.Count());

        // gap> Size(AllAutomorphisms(AbelianGroup([3,5]))); = 8
        var n5 = Product.Generate(new Cn(3), new Cn(5));
        var autN5 = Group.AllAutomorphisms(n5);
        Assert.Equal(8, autN5.Count());

        // gap> Size(AllAutomorphisms(AbelianGroup([5,5]))); = 480
        var n6 = Product.Generate(new Cn(5), new Cn(5));
        var autN6 = Group.AllAutomorphisms(n6);
        Assert.Equal(480, autN6.Count());

        // // gap> Size(AllAutomorphisms(AbelianGroup([3,3,3]))); = 11232
        // var n7 = Product.Generate(new Cn(3), new Cn(3), new Cn(3));
        // var autN7 = Group.AllAutomorphisms(n7);
        // Assert.Equal(11232, autN7.Count()); // 2sec
    }

    [Fact]
    public void Test5GapAllOpsByAutomorphisms()
    {
        var c2 = new Cn(2);
        var c3 = new Cn(3);
        var c4 = new Cn(4);
        // gap> Size(AllHomomorphisms(AbelianGroup([2,2]),AutomorphismGroup(AbelianGroup([2,2,2])))); = 148
        var g1 = Product.Generate(c2, c2);
        var n1 = Product.Generate(c2, c2, c2);
        var allHom1 = Group.AllOpsByAutomorphisms(g1, n1);
        Assert.Equal(148, allHom1.Count);

        // gap> Size(AllHomomorphisms(AbelianGroup([2,2,2]),AutomorphismGroup(AbelianGroup([2,2])))); = 22
        var allHom2 = Group.AllOpsByAutomorphisms(n1, g1);
        Assert.Equal(22, allHom2.Count);

        // gap> Size(AllHomomorphisms(AbelianGroup([2,4]),AutomorphismGroup(AbelianGroup([3,3])))); = 88
        var g3 = Product.Generate(c2, c4);
        var n3 = Product.Generate(c3, c3);
        var allHom3 = Group.AllOpsByAutomorphisms(g3, n3);
        Assert.Equal(88, allHom3.Count);

        // gap> Size(AllHomomorphisms(AbelianGroup([2,2]),AutomorphismGroup(AbelianGroup([3,3])))); = 76
        var g4 = Product.Generate(c2, c2);
        var n4 = Product.Generate(c3, c3);
        var allHom4 = Group.AllOpsByAutomorphisms(g4, n4);
        Assert.Equal(76, allHom4.Count);
    }
}
using System.Linq;
using FastGoat;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Integers;
using Xunit;

namespace Tests;

public class ZnUnitTest
{
    [Fact]
    public void Test1ThrowException()
    {
        var z4 = new Zn(3);
        var z5 = new Zn(5);
        var e1 = z4[2];
        var e2 = z5[3];
        var g = Product.Group(z4, z5);

        Assert.Throws<GroupException>(() => new Zn(-5));
        Assert.Throws<GroupException>(() => z5.Invert(e1));
        Assert.Throws<GroupException>(() => z5.Op(e1, e2));
        Assert.Throws<GroupException>(() => e1.CompareTo(e2));
        Assert.Throws<GroupException>(() => z4.Times(e2, 3));
        Assert.Throws<GroupException>(() => g[(1, 2), 2]);
    }

    [Fact]
    public void Test2CreateElementsAndGroups()
    {
        var z7 = new Zn(7);
        var e1 = z7[3];
        var e2 = z7[6];
        Assert.Equal(e1.BaseGroup, z7);
        Assert.Equal(e1, z7[10]);
        Assert.Equal(z7.Invert(e1), z7[4]);
        Assert.Equal(z7.Op(e1, e2), z7[2]);
        Assert.Equal(z7.Times(e1, 4), z7[5]);
        Assert.Equal(z7.Times(e1, -2), z7[1]);
    }

    [Fact]
    public void Test3SetOfElements()
    {
        var z7 = new Zn(7);
        var set0 = (new[] { 1, 2, 6, 11, 13, 23, 32 }).Select(i => z7[i]).ToHashSet();
        var set1 = (new[] { 1, 2, 6, 4 }).Select(i => z7[i]).ToHashSet();
        Assert.True(set0.SetEquals(set1));

        var arr1 = new[] { 1, 11, 13, 2, 6, 23, 32 };
        var arr2 = new[] { 1, 2, 2, 4, 4, 6, 6 };
        var es1 = arr1.Select(i => z7[i]).Ascending();
        var es2 = arr2.Select(i => z7[i]);
        Assert.True(es1.SequenceEqual(es2));
    }

    [Fact]
    public void Test4Tuples()
    {
        var z5 = new Zn(5);
        var z8 = new Zn(8);
        var g1 = Product.Group(z5, z8);

        var e1 = z5[2];
        var e2 = z8[5];
        var ep1 = g1[2, 5];
        var ep2 = Product.Elt(e1, e2);

        Assert.Equal(g1.G1, z5);
        Assert.Equal(g1.G2, z8);
        Assert.Equal(ep1.E1, e1);
        Assert.Equal(ep1.E2, e2);

        var ep3 = g1[1, 7];
        Assert.Equal(g1.Invert(ep1), g1[3, 3]);
        Assert.Equal(g1.Op(ep1, ep3), g1[3, 4]);
        Assert.Equal(g1.Times(ep1, 3), ep3);

        var z12 = new Zn(12);
        var g2 = Product.Group(z5, z8, z12);
        Assert.Equal(g2.Neutral(), g2[0, 0, 0]);
        Assert.Equal(g2.Invert(g2[2, 3, 10]), g2[3, 5, 2]);
        Assert.Equal(g2.Op(g2[2, 3, 10], g2[4, 7, 5]), g2[1, 2, 3]);
        Assert.Equal(g2.Times(g2[2, 3, 10], 4), g2[3, 4, 4]);

        var g3 = Product.Group(g1, z12);
        Assert.Equal(g3.Neutral(), g3[(0, 0), 0]);
        Assert.Equal(g3.Invert(g3[(2, 3), 10]), g3[(3, 5), 2]);
        Assert.Equal(g3.Op(g3[(2, 3), 10], g3[(4, 7), 5]), g3[(1, 2), 3]);
        Assert.Equal(g3.Times(g3[(2, 3), 10], 4), g3[(3, 4), 4]);

        var g4 = Product.Group(z12, z12, z12, z12);
        var ep4 = Product.Elt(z12[1], z12[3], z12[5], z12[6]);
        Assert.Equal(g4.Neutral(), g4[0, 0, 0, 0]);
        Assert.Equal(g4[11, 9, 7, 6], g4.Invert(ep4));
        Assert.Equal(g4[2, 6, 10, 0], g4.Op(ep4, ep4));
        Assert.Equal(ep4, g4[1, 3, 5, 6]);
        Assert.True(ep4.Equals(g4[1, 3, 5, 6]));
        Assert.Equal(1, ep4.CompareTo(g4[1, 3, 4, 6]));
        Assert.Single(g4);
        Assert.Equal(4, g4.GetGenerators().Count());

        var g5 = Product.Group(z12, z12, z12, z12, z12);
        var ep5 = Product.Elt(z12[1], z12[3], z12[5], z12[6], z12[2]);
        Assert.Equal(g5.Neutral(), g5[0, 0, 0, 0, 0]);
        Assert.Equal(g5[11, 9, 7, 6, 10], g5.Invert(ep5));
        Assert.Equal(g5[2, 6, 10, 0, 4], g5.Op(ep5, ep5));
        Assert.Equal(ep5, g5[1, 3, 5, 6, 2]);
        Assert.True(ep5.Equals(g5[1, 3, 5, 6, 2]));
        Assert.Equal(1, ep5.CompareTo(g5[1, 3, 4, 6, 1]));
        Assert.Single(g5);
        Assert.Equal(5, g5.GetGenerators().Count());
    }

    [Fact]
    public void Test5InfixOperators()
    {
        var z6 = new Zn(6);
        Assert.True(z6[4] == z6[10]);
        Assert.True(z6[4] != z6[5]);
        Assert.False(z6[4] != z6[10]);
        Assert.False(z6[4] == z6[5]);
        Assert.Equal(z6[2], z6[5] + z6[3]);
        Assert.Equal(z6[2], z6[5] * 4);
    }

    [Fact]
    public void Test6Tuples()
    {
        var z5 = new Zn(5);
        var z8 = new Zn(8);
        var g1 = Product.Gp(z5, z8);

        var e1 = z5[2];
        var e2 = z8[5];
        var ep1 = g1[2, 5];
        var ep2 = Product.Ep(e1, e2);

        Assert.Equal(g1.Gi[0], z5);
        Assert.Equal(g1.Gi[1], z8);
        Assert.Equal(ep1.Ei[0], e1);
        Assert.Equal(ep1.Ei[1], e2);

        var ep3 = g1[1, 7];
        Assert.Equal(g1.Invert(ep1), g1[3, 3]);
        Assert.Equal(g1.Op(ep1, ep3), g1[3, 4]);
        Assert.Equal(g1.Times(ep1, 3), ep3);
    }
}
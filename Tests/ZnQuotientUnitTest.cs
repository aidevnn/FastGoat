using System.Linq;
using Xunit;
using FastGoat;
using FastGoat.Gp;
using FastGoat.UserGroup;

namespace Tests;

public class ZnQuotientUnitTest
{
    [Fact]
    public void Test1Cosets()
    {
        var g0 = Group.Generate(new Zn(2)[1]);
        var h0 = Product.Group(g0, g0, g0);
        var h1 = Group.Create("H", h0);
        var h2 = Group.Generate("K", h1, h1[1, 0, 0]);
        var cosets = Group.Cosets(h1, h2);
        Assert.Equal(4, cosets.Count);

        var c000 = cosets[h1[0, 0, 0]];
        Assert.Contains(h1[0, 0, 0], c000);
        Assert.Contains(h1[1, 0, 0], c000);

        var c010 = cosets[h1[0, 1, 0]];
        Assert.Contains(h1[0, 1, 0], c010);
        Assert.Contains(h1[1, 1, 0], c010);

        var c001 = cosets[h1[0, 0, 1]];
        Assert.Contains(h1[0, 0, 1], c001);
        Assert.Contains(h1[1, 0, 1], c001);

        var c011 = cosets[h1[0, 1, 1]];
        Assert.Contains(h1[0, 1, 1], c011);
        Assert.Contains(h1[1, 1, 1], c011);
    }

    [Fact]
    public void Test2QuotientGroup()
    {
        var g0 = Group.Generate(new Zn(2)[1]);
        var klein = Group.Create(Product.Group(g0, g0));
        var h0 = Product.Group(new Zn(4), new Zn(10));
        var h1 = Group.Generate("H", h0, h0[1, 0], h0[0, 1]);
        var h2 = Group.Generate("K", h1, h1[2, 2]);
        var hQuo = h1.Over(h2);
        Assert.True(hQuo.IsIsomorphicTo(klein));
    }
}
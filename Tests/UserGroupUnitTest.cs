using System.Linq;
using FastGoat.Structures;
using FastGoat.UserGroup;
using Xunit;

namespace Tests;

public class UserGroupUnitTest
{
    [Fact]
    public void Test1Dihedral()
    {
        for (int n = 3; n <= 8; n++)
        {
            var dn0 = FG.Dihedral(n);
            var dn1 = FG.DihedralSdp(n);
            var dn2 = FG.DihedralWg(n);
            Assert.True(dn1.IsIsomorphicTo(dn0));
            Assert.True(dn1.IsIsomorphicTo(dn2));
        }
    }

    [Fact]
    public void Test2Alternate()
    {
        Assert.Equal(3, FG.Alternate(3).Count());
        Assert.Equal(6, FG.Symmetric(3).Count());
        
        Assert.Equal(12, FG.Alternate(4).Count());
        Assert.Equal(24, FG.Symmetric(4).Count());
        
        Assert.Equal(60, FG.Alternate(5).Count());
        Assert.Equal(120, FG.Symmetric(5).Count());
        
        Assert.Equal(360, FG.Alternate(6).Count());
        Assert.Equal(720, FG.Symmetric(6).Count());
    }

    [Fact]
    public void Test3Abelian()
    {
        Assert.Equal(9, FG.Abelian(3, 3).Count());
        Assert.Equal(150, FG.Abelian(10, 15).Count());
        Assert.Equal(180, FG.Abelian(3, 6, 10).Count());
        Assert.Equal("H1", FG.Abelian("H1", 3, 5).Name);
        Assert.Equal("K2", FG.Abelian("K2", 7, 11, 15).Name);
    }

    [Fact]
    public void Test4PermGroup()
    {
        Assert.Equal(4, FG.PermGroup(4, ((1, 2), (3, 4)), ((1, 3), (2, 4))).Count());
        Assert.Equal(8, FG.PermGroup(4, (1, 2), ((1, 3), (2, 4))).Count());
        Assert.Equal(12, FG.PermGroup(4, ((1, 2), (3, 4)), (1, 2, 3)).Count());
        Assert.Equal("S3", FG.PermGroup("S3", 3, (1, 2), (1, 2, 3)).Name);
    }

    [Fact]
    public void Test5Frobenius()
    {
        var f20 = FG.Frobenius(20);
        Assert.Equal(2, f20.Count);
        Assert.Equal(20, f20[0].Count());
        Assert.Equal(20, f20[1].Count());
        
        var f20sdp = FG.FrobeniusSdp(20);
        Assert.True(f20[0].IsIsomorphicTo(f20sdp[0]));
        Assert.True(f20[1].IsIsomorphicTo(f20sdp[1]));
        
        var f21 = FG.Frobenius(21);
        Assert.Single(f21);
        Assert.Equal(21, f21[0].Count());

        var f21sdp = FG.FrobeniusSdp(21);
        Assert.True(f21[0].IsIsomorphicTo(f21sdp[0]));
    }

    [Fact]
    public void Test6DiCyclic()
    {
        var d2 = FG.DiCyclic(2);
        var d2sdp = FG.DiCyclicSdp(2);
        Assert.Equal(8, d2.Count());
        Assert.True(d2.IsIsomorphicTo(d2sdp));

        var d3 = FG.DiCyclic(3);
        var d3sdp = FG.DiCyclicSdp(3);
        Assert.Equal(12, d3.Count());
        Assert.True(d3.IsIsomorphicTo(d3sdp));

        var d5 = FG.DiCyclic(5);
        var d5sdp = FG.DiCyclicSdp(5);
        Assert.Equal(20, d5.Count());
        Assert.True(d5.IsIsomorphicTo(d5sdp));
    }
}
using System.Linq;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.UserGroup.Perms;
using Xunit;

namespace Tests;

public class SnUnitTest
{
    [Fact]
    public void Test1ThrowException()
    {
        var s3 = new Sn(3);
        var e1 = s3.Neutral();
        var s4 = new Sn(4);
        var e2 = s4.Neutral();
        Assert.Throws<GroupException>(() => new Sn(-5));
        Assert.Throws<GroupException>(() => s4.Invert(e1));
        Assert.Throws<GroupException>(() => s3.Op(e1, e2));
        Assert.Throws<GroupException>(() => e1.CompareTo(e2));
        Assert.Throws<GroupException>(() => s3.Times(e2, 3));
        Assert.Throws<GroupException>(() => s4[(4, 3, 4)]);
        Assert.Throws<GroupException>(() => s4[(2, 6)]);
        Assert.Throws<GroupException>(() => s4[(2, 1), 3]);
        Assert.Throws<GroupException>(() => s4[((2, 1), (3, 4))]);
    }

    [Fact]
    public void Test2Create()
    {
        var s4 = new Sn(4);
        var e1 = s4[(1, 2)];
        var e2 = s4[(1, 3), (2, 4)];
        var e3 = s4[(1, 2, 4)];
        var e4 = s4[(4, 2, 1)];
        Assert.True(e1.Table.SequenceEqual(new[] { 1, 0, 2, 3 }));
        Assert.True(e2.Table.SequenceEqual(new[] { 2, 3, 0, 1 }));
        Assert.Equal(e1, s4.Invert(e1));
        Assert.Equal(e3, s4.Invert(e4));
        Assert.Equal(s4[(1, 4)], s4.Op(e1, e3));
        Assert.Equal(e4, s4.Times(e3, 2));
        Assert.Equal(s4.Neutral(), s4.Times(e3, -6));
        Assert.Equal(s4, e3.BaseGroup);
    }

    [Fact]
    public void Test3AllElements()
    {
        var s3 = new Sn(3);
        var s4 = new Sn(4);
        var s5 = new Sn(5);
        var s6 = new Sn(6);
        var s7 = new Sn(7);

        var S3 = Group.Generate(s3, s3[(1, 2)], s3[(1, 2, 3)]).ToHashSet();
        var S4 = Group.Generate(s4, s4[(1, 2)], s4[(1, 2, 3, 4)]).ToHashSet();
        var S5 = Group.Generate(s5, s5[(1, 2)], s5[(1, 2, 3, 4, 5)]).ToHashSet();
        var S6 = Group.Generate(s6, s6[(1, 2)], s6[(1, 2, 3, 4, 5, 6)]).ToHashSet();
        var S7 = Group.Generate(s7, s7[(1, 2)], s7[(1, 2, 3, 4, 5, 6, 7)]).ToHashSet();

        Assert.True(S3.SetEquals(IntExt.GetPermutations(3).Select(s3.CreateElement)));
        Assert.True(S4.SetEquals(IntExt.GetPermutations(4).Select(s4.CreateElement)));
        Assert.True(S5.SetEquals(IntExt.GetPermutations(5).Select(s5.CreateElement)));
        Assert.True(S6.SetEquals(IntExt.GetPermutations(6).Select(s6.CreateElement)));
        Assert.True(S7.SetEquals(IntExt.GetPermutations(7).Select(s7.CreateElement)));
    }

    [Fact]
    public void Test4InfixOperators()
    {
        var s3 = new Sn(3);
        Assert.True(s3[(1, 2)] == s3[(2, 1)]);
        Assert.True(s3[(1, 2)] != s3[(1, 3)]);
        Assert.False(s3[(1, 2)] != s3[(2, 1)]);
        Assert.False(s3[(1, 2)] == s3[(1, 3)]);
        Assert.Equal(s3[(1, 2, 3)], s3[(1, 2)] * s3[(1, 3)]);
        Assert.Equal(s3.Neutral(), s3[(1, 2, 3)] ^ 3);
    }

    [Fact]
    public void Test5Symm()
    {
        var s3 = new Symm(3);
        Assert.Contains(new[]
            {
                s3.Neutral(),
                s3[(1, 2)], s3[(1, 3)], s3[(2, 3)],
                s3[(1, 2, 3)], s3[(1, 3, 2)]
            }
            , s3.Contains
        );
        Assert.Equal(3, s3.SubGroupsGenerators().Count());
    }
}
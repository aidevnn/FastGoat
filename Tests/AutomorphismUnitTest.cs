using System.Linq;
using Xunit;
using FastGoat;
using FastGoat.UserGroup;

namespace Tests;

public class AutomorphismUnitTest
{
    [Fact]
    public void Test1AutC2()
    {
        var c2 = new Cn(2);
        var gbAut = Group.Aut(c2);
        var id2 = gbAut.Neutral();
        Assert.Equal(c2[0], id2[c2[0]]);
        Assert.Equal(c2[1], id2[c2[1]]);
        Assert.Throws<GroupException>(() => gbAut[(c2[0], c2[1])]);
    }

    [Fact]
    public void Test2AutC3()
    {
        var c3 = new Cn(3);
        var gbAut = Group.Aut(c3);
        var id3 = gbAut.Neutral();
        var y = gbAut[(c3[1], c3[2])];
        var yi = gbAut.Invert(y);
        Assert.Equal(id3, gbAut.Op(y, yi));

        var gAut = Group.Generate(y);
        Assert.NotEqual(3, gAut.Count());
        Assert.Contains(id3, gAut);
        Assert.Contains(y, gAut);
        Assert.Contains(yi, gAut);
        Assert.Equal(yi, y);
        Assert.Equal(2, gAut.Count());
    }

    [Fact]
    public void Test3AutC5()
    {
        var c5 = new Cn(5);
        var gbAut = Group.Aut(c5);
        var y = gbAut[(c5[1], c5[2])];
        Assert.Equal(c5[0], y[c5[0]]);
        Assert.Equal(c5[2], y[c5[1]]);
        Assert.Equal(c5[4], y[c5[2]]);
        Assert.Equal(c5[1], y[c5[3]]);
        Assert.Equal(c5[3], y[c5[4]]);

        var gAut = Group.Generate(y);
        Assert.Equal(4, gAut.Count());
        var y2 = gbAut[(c5[1], c5[4])];
        var y3 = gbAut[(c5[1], c5[3])];
        Assert.Contains(y2, gAut);
        Assert.Contains(y3, gAut);
        Assert.Equal(y2, gAut.Times(y, 2));
        Assert.Equal(y3, gAut.Times(y, 3));
        Assert.Equal(y3, gAut.Invert(y));
        Assert.Equal(y2, gAut.Invert(y2));

        Assert.Throws<GroupException>(() => gbAut[(c5[1], c5[2]), (c5[2], c5[3])]);
    }

    [Fact]
    public void Test4AutS3()
    {
        var s3 = new Sn(3);
        var S3 = Group.Generate(s3[(1, 2)], s3[(1, 3)]);
        var gbAut = Group.Aut(S3);
        var y1 = gbAut[(S3[(1, 2)], S3[(2, 3)]), (S3[(1, 3)], S3[(1, 2)])];
        var y2 = gbAut[(S3[(1, 2)], S3[(1, 3)]), (S3[(1, 3)], S3[(1, 2)])];
        var AutS3 = Group.Generate(y1, y2);
        Assert.True(AutS3.IsIsomorphicTo(S3));
    }
}
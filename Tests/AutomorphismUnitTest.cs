using System.Linq;
using Xunit;
using FastGoat;
using FastGoat.Gp;
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

    [Fact]
    public void Test5Un()
    {
        Assert.Single(new Un(2));
        Assert.Equal(2, new Un(3).Count());
        Assert.True(new Un(4).IsIsomorphicTo(new Un(3)));
        Assert.True(new Un(5).IsIsomorphicTo(new Cn(4)));
        var klein = Group.Create(Product.Group(new Cn(2), new Cn(2)));
        Assert.True(new Un(8).IsIsomorphicTo(klein));
    }

    [Fact]
    public void Test6AutCn()
    {
        for (int k = 2; k < 33; ++k)
        {
            var un = new Un(k);
            var cn = new Cn(k);
            var autCn = Group.Aut(cn[1]);
            Assert.True(un.IsIsomorphicTo(autCn));
        }
    }

    [Fact]
    public void Test6AutCmXCn()
    {
        int N = 10;
        for (int m = 2; m < N; ++m)
        {
            var um = new Un(m);
            var cm = new Cn(m);
            for (int n = 2; n < N; ++n)
            {
                if (IntExt.Gcd(m, n) != 1)
                    continue;

                var un = new Un(n);
                var cn = new Cn(n);
                var gmn = Product.Group(cm, cn);
                var autCmXCn = Group.Aut(gmn[1, 0], gmn[0, 1]);
                var umn = Product.Generate(um, un);
                Assert.Equal(um.Count() * un.Count(), autCmXCn.Count());
                Assert.True(umn.IsIsomorphicTo(autCmXCn));
            }
        }
    }

    [Fact]
    public void Test7GapValidation()
    {
        var g1 = Product.Group(new Cn(3), new Cn(3));
        var autG1 = Group.Aut(g1[1, 0], g1[0, 1]);
        Assert.Equal(48, autG1.Count());

        var g2 = Product.Group(new Cn(2), new Cn(2), new Cn(2));
        var autG2 = Group.Aut(g2[1, 0, 0], g2[0, 1, 0], g2[0, 0, 1]);
        Assert.Equal(168, autG2.Count());

        var g3 = Product.Group(new Cn(2), new Cn(8));
        var autG3 = Group.Aut(g3[1, 0], g3[0, 1]);
        Assert.Equal(16, autG3.Count());

        var g4 = Product.Group(new Cn(2), new Cn(2), new Cn(4));
        var autG4 = Group.Aut(g4[1, 0, 0], g4[0, 1, 0], g4[0, 0, 1]);
        Assert.Equal(192, autG4.Count());
    }
}
using Xunit;
using FastGoat;

namespace Tests;

public class ZnUnitTest
{
    [Fact]
    public void CreateElementTest1()
    {
        var z5 = new Zn(5);
        ZnInt a = (z5, 3);
        ZnInt b = (z5, 8);
        Assert.Equal(a, b);
        Assert.Equal(a ^ 3, z5[4]);
        Assert.Equal(a ^ 5, z5.Neutral());
    }
    [Fact]
    public void CreateElementTest2()
    {
        var z4 = new Zn(4);
        var z5 = new Zn(5);
        var z4xz5 = Group.CartesianProduct(z4, z5);
        var a = z4xz5[3, 0];
        Ep<ZnInt, ZnInt> a0 = ((z4, 3), (z5, 0));
        var b = z4xz5[2, 2];
        Assert.Equal(a, a0);
        Assert.Equal(z4xz5.Invert(a), z4xz5[1, 0]);
        Assert.Equal(z4xz5.Op(a, b), z4xz5[1, 2]);
        Assert.Equal(z4xz5.Op(a, b), z4xz5[1, 2]);
    }
}
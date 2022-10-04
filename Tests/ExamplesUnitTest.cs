using System.Linq;
using FastGoat;
using FastGoat.Examples;
using FastGoat.Gp;
using FastGoat.UserGroup;
using Xunit;

namespace Tests;

public class ExamplesUnitTest
{
    [Fact]
    public void Test1Group21()
    {
        GroupOrder21.Zn21();
        GroupOrder21.ZnDirectProduct();

        GroupOrder21.Symmetric7();
        GroupOrder21.Symmetric7Fast();

        GroupOrder21.SemiDirectProduct();
        Assert.True(true);
    }

    [Fact]
    public void Test1Invariants()
    {
        var c14 = new Cn(14);
        var c21 = new Cn(21);
        var bg = Product.Group(c14, c21);
        var g = Group.Create(bg.Name, bg);
        var decomposition = AbelianInvariantsFactors.Reduce(g);
        Assert.True(decomposition.SequenceEqual(new[] { 7, 42 }));
    }
}
using System;
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
        bool Test()
        {
            try
            {
                GroupOrder21.Zn21();
                GroupOrder21.ZnDirectProduct();
                GroupOrder21.Symmetric7();
                GroupOrder21.Symmetric7Fast();
                GroupOrder21.SemiDirectProduct();
                return true;
            }
            catch (Exception e)
            {
                return false;
            }
        }

        Assert.True(Test());
    }

    [Fact]
    public void Test2Invariants()
    {
        var c14 = new Cn(14);
        var c21 = new Cn(21);
        var bg = Product.Group(c14, c21);
        var g = Group.Create(bg.Name, bg);
        var decomposition = AbelianInvariantsFactors.Reduce(g);
        Assert.True(decomposition.SequenceEqual(new[] { 7, 42 }));

        bool Test()
        {
            try
            {
                AbelianInvariantsFactors.Reduce(g);
                AbelianInvariantsFactors.InvariantFactors294();
                AbelianInvariantsFactors.InvariantFactors600();
                AbelianInvariantsFactors.InvariantFactors4320();
                return true;
            }
            catch (Exception e)
            {
                return false;
            }
        }

        Assert.True(Test());
    }

    [Fact]
    public void Test3GroupAction()
    {
        bool Test()
        {
            try
            {
                GroupAction.C7_C3();
                GroupAction.GroupOrder72();
                return true;
            }
            catch (Exception e)
            {
                return false;
            }
        }

        Assert.True(Test());
    }

    [Fact]
    public void Test4GroupOrder32()
    {
        bool Test()
        {
            try
            {
                GroupOrder32C4C8.Isomorphic();
                GroupOrder32C4C8.IsomorphicAnotherExample();
                GroupOrder32C4C8.NonIsomorphic();
                GroupOrder32C4C8.NonIsomorphicAnotherExample();
                return true;
            }
            catch (Exception e)
            {
                return false;
            }
        }

        Assert.True(Test());
    }

    [Fact]
    public void Test5DihedralAut()
    {
        bool Test()
        {
            try
            {
                DihedralAutomorphisms.Dn();
                DihedralAutomorphisms.AutDn();
                return true;
            }
            catch (Exception e)
            {
                return false;
            }
        }

        Assert.True(Test());
    }
}
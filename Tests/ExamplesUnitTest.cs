using System;
using System.Linq;
using FastGoat;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Examples;
using FastGoat.Structures;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Integers;
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
                DihedralAutomorphisms.AutDnPerm();
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
    public void Test6HolomorphC7()
    {
        bool Test()
        {
            try
            {
                HolomorphC7.Sample();
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
    public void Test7Group18()
    {
        bool Test()
        {
            try
            {
                GroupOrder18.Abelians18();
                GroupOrder18.D18();
                GroupOrder18.C3sdpS3();
                GroupOrder18.C3xS3();
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
    public void Test8WordGroup()
    {
        bool Test()
        {
            try
            {
                ToddCoxeter.Frobenius20();
                ToddCoxeter.Quaternion();
                ToddCoxeter.Quaternion32();
                ToddCoxeter.CyclicGroup();
                ToddCoxeter.DiCyclic3();
                ToddCoxeter.DiCyclic12();
                ToddCoxeter.KleinGroup();
                ToddCoxeter.MoreExamples();
                ToddCoxeter.SmallGroup32_32();
                ToddCoxeter.Symm3Group();
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
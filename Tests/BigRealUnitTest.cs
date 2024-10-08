using System;
using System.Linq;
using System.Numerics;
using FastGoat.Commons;
using FastGoat.UserGroup.Floats;
using Xunit;

namespace Tests;

public class BigRealUnitTest
{
    BigReal GenBR((int M0, int M1) m, (int E0, int E1) e, int O)
    {
        var a0 = IntExt.Rng.Next(m.M0, m.M1 + 1);
        var e0 = IntExt.Rng.Next(e.E0, e.E1 + 1);
        return BigReal.FromBigInteger(a0, O).Mul10PowN(e0);
    }

    string fmt(double a, int O) => String.Format($"{{0:E{O}}}", a);

    [Fact]
    public void TestFromBigInteger()
    {
        var k = 2000;
        var M = 1000;
        var E = 6;
        var O = 12;
        var lt = k.Range().Select(_ => IntExt.Rng.Next(-M * M, M * M) * BigInteger.Pow(10, IntExt.Rng.Next(0, E)))
            .ToArray();
        var err = Double.Pow(10, -O + 2);
        foreach (var a in lt)
        {
            var a_b = BigReal.FromBigInteger(a, O);
            var a_b2 = Double.Parse($"{a}");
            var diff = Double.Abs(a_b.ToDouble - a_b2) * Double.Pow(10, -a_b.V);
            var str = new
            {
                a, a_b.Details, a_b, a_bd = $"{fmt(a_b.ToDouble, O)}", a_bd2 = $"{fmt(a_b2, O)}", diff, epsilon = err
            };
            Assert.False(diff > err, $"{str}");
        }
    }

    [Fact]
    public void TestFromString()
    {
        var k = 2000;
        var M = 1000;
        var E = 6;
        var O = 12;
        var lt = k.Range().Select(_ => GenBR((-M * M, M * M), (-2 * E, 2 * E), O)).ToArray();
        var err = Double.Pow(10, -O + 2);
        foreach (var a in lt)
        {
            var a_b = a;
            var a_b2 = Double.Parse($"{a}");
            var diff = Double.Abs(a_b.ToDouble - a_b2) * Double.Pow(10, -a_b.V);
            var str1 = new { a, ad = $"{fmt(a.ToDouble, O)}", a0 = a.Details };
            var str2 = new { a_b, a_bd = $"{fmt(a_b.ToDouble, O)}", a_bd2 = $"{fmt(a_b2, O)}", diff, epsilon = err };
            Assert.False(diff > err, $"{str1} ### {str2}");
        }
    }

    [Fact]
    public void TestToRational()
    {
        var k = 2000;
        var M = 1000;
        var E = 6;
        var O = 12;
        var lt = k.Range().Select(_ => GenBR((-M * M, M * M), (-2 * E, 2 * E), O)).ToArray();
        var err = Double.Pow(10, -O + 2);
        foreach (var a in lt)
        {
            var a_b = a;
            var a_b2 = (double)a.ToRational;
            var diff = Double.Abs(a_b.ToDouble - a_b2) * Double.Pow(10, -a_b.V);
            var str = new
            {
                a, a_b.Details, a_b, a_bd = $"{fmt(a_b.ToDouble, O)}", a_bd2 = $"{fmt(a_b2, O)}", diff, epsilon = err
            };
            Assert.False(diff > err, $"{str}");
        }
    }

    [Fact]
    public void TestAdditionSubstraction()
    {
        var k = 2000;
        var M = 1000;
        var E = 6;
        var O = 12;
        var lt = k.Range()
            .Select(_ => (GenBR((-M * M, M * M), (-2 * E, 2 * E), O), GenBR((-M, M), (-E, E), O)))
            .Select(e =>
                IntExt.Rng.NextDouble() < 0.025 ? (e.Item1, e.Item1) :
                IntExt.Rng.NextDouble() < 0.5 ? e : (e.Item2, e.Item1))
            .ToArray();
        var err = Double.Pow(10, -O + 2);
        foreach (var (a, b) in lt)
        {
            var a_b = a + b;
            var a_b2 = a.ToDouble + b.ToDouble;
            var diff = Double.Abs(a_b.ToDouble - a_b2) * Double.Pow(10, -a_b.V);
            var str1 = new
            {
                a, ad = $"{fmt(a.ToDouble, O)}", a0 = a.Details, b, bd = $"{fmt(b.ToDouble, O)}", b0 = b.Details
            };
            var str2 = new { a_b, a_bd = $"{fmt(a_b.ToDouble, O)}", a_bd2 = $"{fmt(a_b2, O)}", diff, epsilon = err };
            Assert.False(diff > err, $"{str1} ### {str2}");
        }
    }

    [Fact]
    public void TestMultiplication()
    {
        var k = 2000;
        var M = 1000;
        var E = 6;
        var O = 12;
        var lt = k.Range()
            .Select(_ => (GenBR((-M * M, M * M), (-2 * E, 2 * E), O), GenBR((-M, M), (-E, E), O)))
            .Select(e =>
                IntExt.Rng.NextDouble() < 0.025 ? (e.Item1, e.Item1) :
                IntExt.Rng.NextDouble() < 0.5 ? e : (e.Item2, e.Item1))
            .ToArray();
        var err = Double.Pow(10, -O + 2);
        foreach (var (a, b) in lt)
        {
            var a_b = a * b;
            var a_b2 = a.ToDouble * b.ToDouble;
            var diff = Double.Abs(a_b.ToDouble - a_b2) * Double.Pow(10, -a_b.V);
            var str1 = new
                { a, ad = $"{fmt(a.ToDouble, O)}", a0 = a.Details, b, bd = $"{fmt(b.ToDouble, O)}", b0 = b.Details };
            var str2 = new { a_b, a_bd = $"{fmt(a_b.ToDouble, O)}", a_bd2 = $"{fmt(a_b2, O)}", diff, epsilon = err };
            Assert.False(diff > err, $"{str1} ### {str2}");
        }
    }

    [Fact]
    public void TestDivision()
    {
        var k = 2000;
        var M = 1000;
        var E = 6;
        var O = 12;
        var lt = k.Range()
            .Select(_ => (GenBR((-M * M, M * M), (-2 * E, 2 * E), O), GenBR((-M, M), (-E, E), O)))
            .Select(e =>
                IntExt.Rng.NextDouble() < 0.025 ? (e.Item1, e.Item1) :
                IntExt.Rng.NextDouble() < 0.5 ? e : (e.Item2, e.Item1))
            .ToArray();
        GlobalStopWatch.AddLap();
        var err = Double.Pow(10, -O + 2);
        foreach (var (a, b) in lt)
        {
            if (b.IsZero())
                continue;

            var a_b = a / b;
            var a_b2 = a.ToDouble / b.ToDouble;
            var diff = Double.Abs(a_b.ToDouble - a_b2) * Double.Pow(10, -a_b.V);
            var str1 = new { a, ad = $"{fmt(a.ToDouble, O)}", a0 = a.Details, b };
            var str2 = new { a_b, a_bd = $"{fmt(a_b.ToDouble, O)}", a_bd2 = $"{fmt(a_b2, O)}", diff, epsilon = err };
            Assert.False(diff > err, $"{str1} ### {str2}");
        }
    }

    [Fact]
    public void TestRoundSciForm()
    {
        var form = Rounding.SciForm;
        var a = BigReal.FromString("3.67849275e-5", 15);

        Assert.Equal(BigReal.FromString("4e-5", 15), BigReal.Round(a, 0, MidpointRounding.ToEven, form));
        Assert.Equal(BigReal.FromString("3.7e-5", 15), BigReal.Round(a, 1, MidpointRounding.ToEven, form));
        Assert.Equal(BigReal.FromString("3.68e-5", 15), BigReal.Round(a, 2, MidpointRounding.ToEven, form));
        Assert.Equal(BigReal.FromString("3.678e-5", 15), BigReal.Round(a, 3, MidpointRounding.ToEven, form));
        Assert.Equal(BigReal.FromString("3.6785e-5", 15), BigReal.Round(a, 4, MidpointRounding.ToEven, form));
        Assert.Equal(BigReal.FromString("3.67849e-5", 15), BigReal.Round(a, 5, MidpointRounding.ToEven, form));
        Assert.Equal(BigReal.FromString("3.678493e-5", 15), BigReal.Round(a, 6, MidpointRounding.ToEven, form));
        Assert.Equal(BigReal.FromString("3.67849275e-5", 15), BigReal.Round(a, 8, MidpointRounding.ToEven, form));
        
        // last digit
        Assert.Equal(BigReal.FromString("3.6784928e-5", 15), BigReal.Round(a, 7, MidpointRounding.ToEven, form));
        Assert.Equal(BigReal.FromString("3.6784927e-5", 15), BigReal.Round(a, 7, MidpointRounding.ToZero, form));
        Assert.Equal(BigReal.FromString("3.6784928e-5", 15), BigReal.Round(a, 7, MidpointRounding.AwayFromZero, form));
        Assert.Equal(BigReal.FromString("3.6784927e-5", 15), BigReal.Round(a, 7, MidpointRounding.ToNegativeInfinity, form));
        Assert.Equal(BigReal.FromString("3.6784928e-5", 15), BigReal.Round(a, 7, MidpointRounding.ToPositiveInfinity, form));
        
        Assert.Equal(BigReal.FromString("-3.6784928e-5", 15), BigReal.Round(-a, 7, MidpointRounding.ToEven, form));
        Assert.Equal(BigReal.FromString("-3.6784927e-5", 15), BigReal.Round(-a, 7, MidpointRounding.ToZero, form));
        Assert.Equal(BigReal.FromString("-3.6784928e-5", 15), BigReal.Round(-a, 7, MidpointRounding.AwayFromZero, form));
        Assert.Equal(BigReal.FromString("-3.6784928e-5", 15), BigReal.Round(-a, 7, MidpointRounding.ToNegativeInfinity, form));
        Assert.Equal(BigReal.FromString("-3.6784927e-5", 15), BigReal.Round(-a, 7, MidpointRounding.ToPositiveInfinity, form));
    }
    
    [Fact]
    public void TestRoundFixForm()
    {
        var form = Rounding.FixForm;
        var a = BigReal.FromString("3.67849275e-5", 15);

        Assert.Equal(BigReal.BrZero(15), BigReal.Round(a, 0, MidpointRounding.ToEven, form));
        Assert.Equal(BigReal.BrZero(15), BigReal.Round(a, 1, MidpointRounding.ToEven, form));
        Assert.Equal(BigReal.FromString("4e-5", 15), BigReal.Round(a, 5, MidpointRounding.ToEven, form));
        Assert.Equal(BigReal.FromString("3.7e-5", 15), BigReal.Round(a, 6, MidpointRounding.ToEven, form));
        Assert.Equal(BigReal.FromString("3.68e-5", 15), BigReal.Round(a, 7, MidpointRounding.ToEven, form));
        Assert.Equal(BigReal.FromString("3.678e-5", 15), BigReal.Round(a, 8, MidpointRounding.ToEven, form));
        Assert.Equal(BigReal.FromString("3.6785e-5", 15), BigReal.Round(a, 9, MidpointRounding.ToEven, form));
        Assert.Equal(BigReal.FromString("3.67849e-5", 15), BigReal.Round(a, 10, MidpointRounding.ToEven, form));
        Assert.Equal(BigReal.FromString("3.678493e-5", 15), BigReal.Round(a, 11, MidpointRounding.ToEven, form));
        Assert.Equal(BigReal.FromString("3.67849275e-5", 15), BigReal.Round(a, 13, MidpointRounding.ToEven, form));
        
        // last digit
        Assert.Equal(BigReal.FromString("3.6784928e-5", 15), BigReal.Round(a, 12, MidpointRounding.ToEven, form));
        Assert.Equal(BigReal.FromString("3.6784927e-5", 15), BigReal.Round(a, 12, MidpointRounding.ToZero, form));
        Assert.Equal(BigReal.FromString("3.6784928e-5", 15), BigReal.Round(a, 12, MidpointRounding.AwayFromZero, form));
        Assert.Equal(BigReal.FromString("3.6784927e-5", 15), BigReal.Round(a, 12, MidpointRounding.ToNegativeInfinity, form));
        Assert.Equal(BigReal.FromString("3.6784928e-5", 15), BigReal.Round(a, 12, MidpointRounding.ToPositiveInfinity, form));
        
        Assert.Equal(BigReal.FromString("-3.6784928e-5", 15), BigReal.Round(-a, 12, MidpointRounding.ToEven, form));
        Assert.Equal(BigReal.FromString("-3.6784927e-5", 15), BigReal.Round(-a, 12, MidpointRounding.ToZero, form));
        Assert.Equal(BigReal.FromString("-3.6784928e-5", 15), BigReal.Round(-a, 12, MidpointRounding.AwayFromZero, form));
        Assert.Equal(BigReal.FromString("-3.6784928e-5", 15), BigReal.Round(-a, 12, MidpointRounding.ToNegativeInfinity, form));
        Assert.Equal(BigReal.FromString("-3.6784927e-5", 15), BigReal.Round(-a, 12, MidpointRounding.ToPositiveInfinity, form));
    }
}
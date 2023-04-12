using System;
using FastGoat.Structures;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Padic;
using Xunit;

namespace Tests;

public class PadicUnitTest
{
    [Fact]
    public void Test1ThrowException()
    {
        Assert.Throws<ArgumentException>(() => Valuation.Infinity.V);
        Assert.Throws<ArgumentException>(() => ZnBInt.ZnZero(1));
        Assert.Throws<ArgumentException>(() => new PadicGAP(6, 5));
        Assert.Throws<ArgumentException>(() => new PadicGAP(3, -5));
        Assert.Throws<ArgumentException>(() => new PadicZealous(-3, -5));
    }

    [Fact]
    public void Test2ZnBInt()
    {
        var mod1 = 17;

        var z00 = new ZnInt(mod1, 12);
        var z01 = new ZnInt(mod1, 7);
        var z02 = new ZnInt(mod1, 8);
        var z10 = new ZnBInt(mod1, 12);
        var z11 = new ZnBInt(mod1, 7);
        var z12 = new ZnBInt(mod1, 8);

        Assert.Equal(z00.K, z10.K);
        Assert.Equal(-z00.K, -z10.K);
        Assert.Equal((z00 + z01).K, (z10 + z11).K);
        Assert.Equal((z00 - z01).K, (z10 - z11).K);
        Assert.Equal(z00.Pow(3).K, z10.Pow(3).K);
        Assert.Equal((z00 * z01).K, (z10 * z11).K);
        Assert.Equal((1 / z00).K, (1 / z10).K);
        Assert.Equal((z00 / z01).K, (z10 / z11).K);

        Assert.Equal(Ring.Gcd(z01, z02).K, Ring.Gcd(z11, z12).K);
        var bez0 = Ring.Bezout(z01, z02);
        var bez1 = Ring.Bezout(z11, z12);
        Assert.Equal((bez0.x.K, bez0.y.K), (bez1.x.K, bez1.y.K));
    }

    [Fact]
    public void Test3PadicIntegers()
    {
        var a0 = new PadicGAP(5, 3);
        var a1 = new PadicGAP(5, 3, 101);
        var a2 = new PadicGAP(5, 3, 75);
        var a3 = new PadicGAP(5, 3, 15);
        var a4 = a3 / 125;

        Assert.Equal(0, a0.K);
        Assert.Equal(101, a1.K);
        Assert.Equal(3, a2.K);
        Assert.Equal(3, a3.K);
        Assert.Equal(3, a4.K);
        Assert.True(a0.Val.IsInfinity);
        Assert.Equal(new Valuation(0), a1.Val);
        Assert.Equal(new Valuation(2), a2.Val);
        Assert.Equal(new Valuation(1), a3.Val);
        Assert.Equal(new Valuation(-2), a4.Val);

        Assert.Equal(a1.Zero, (a1 + a1.Opp()));
        Assert.Equal(a1.One, (a1 * a1.Inv()));
        Assert.Equal(new PadicGAP(5, 3, 51), a1 + a2);
        Assert.Equal(new PadicGAP(5, 3, 26), a1 - a2);
        Assert.Equal(new PadicGAP(5, 3, 53, 2), a1 * a2); // val(75) = 2 and 101*3 mod 125 = 303 mod 125 = 53 
        Assert.Equal(new PadicGAP(5, 3, 117, -2),
            a1 / a2); // val(75) = 2 501=1+125*4=3*167 and 101*167=16867 mod 125=117
        Assert.Equal(new PadicGAP(5, 3, 76), a1.Pow(2)); // 101^2=10201 mod 125 = 76
    }

    [Fact]
    public void Test4PadicRational()
    {
        var a1 = new PadicGAP(5, 7, 171);
        var a2 = new PadicGAP(5, 7, new Rational(171, 9));
        var a3 = new PadicGAP(5, 7, new Rational(21, 9));
        var a4 = new PadicGAP(5, 7, new Rational(15, 11));
        var a5 = new PadicGAP(5, 7, new Rational(21, 9) + new Rational(15, 11));
        var a6 = new PadicGAP(5, 7, new Rational(21, 9) - new Rational(15, 11));
        var a7 = new PadicGAP(5, 7, new Rational(21, 9) * new Rational(15, 11));
        var a8 = new PadicGAP(5, 7, new Rational(21, 9) / new Rational(15, 11));

        Assert.Equal(a1 / 9, a2);
        Assert.Equal(a2.Normalized, a2 * a2.Norm);
        Assert.Equal(a5, a3 + a4);
        Assert.Equal(a6, a3 - a4);
        Assert.Equal(a7, a3 * a4);
        Assert.Equal(a8, a3 / a4);
    }

    [Fact]
    public void Test5PadicZealousIntegers()
    {
        var a0 = new PadicZealous(5, 3);
        var a1 = new PadicZealous(5, 3, 101);
        var a2 = new PadicZealous(5, 3, 75);
        var a3 = new PadicZealous(5, 3, 15);
        var a4 = a3 / 125;

        Assert.Equal(0, a0.S);
        Assert.Equal(101, a1.S);
        Assert.Equal(3, a2.S);
        Assert.Equal(3, a3.S);
        Assert.Equal(0, a4.S);
        Assert.Equal(a0.N, a0.Val);
        Assert.Equal(0, a1.Val);
        Assert.Equal(2, a2.Val);
        Assert.Equal(1, a3.Val);
        Assert.Equal(a4.N, a4.Val);

        Assert.Equal(a1.Zero, (a1 + a1.Opp()));
        Assert.Equal(a1.One, (a1 * a1.Inv()));
        Assert.Equal(new PadicZealous(5, 3, 51), a1 + a2);
        Assert.Equal(new PadicZealous(5, 3, 26), a1 - a2);
        Assert.Equal(new PadicZealous(5, 3, 53) * 25, a1 * a2);
        Assert.Equal(new PadicZealous(5, 3, 117) / 25, a1 / a2);
        Assert.Equal(new PadicZealous(5, 3, 76), a1.Pow(2));
    }

    [Fact]
    public void Test6PadicZealousRational()
    {
        var a1 = new PadicZealous(5, 7, 171);
        var a2 = new PadicZealous(5, 7, new Rational(171, 9));
        var a3 = new PadicZealous(5, 7, new Rational(21, 9));
        var a4 = new PadicZealous(5, 7, new Rational(15, 11));
        var a5 = new PadicZealous(5, 7, new Rational(21, 9) + new Rational(15, 11));
        var a6 = new PadicZealous(5, 7, new Rational(21, 9) - new Rational(15, 11));
        var a7 = new PadicZealous(5, 7, new Rational(21, 9) * new Rational(15, 11));
        var a8 = new PadicZealous(5, 7, new Rational(21, 9) / new Rational(15, 11));

        Assert.Equal(a1 / 9, a2);
        Assert.Equal(a2.Normalized, a2 * a2.Norm);
        Assert.Equal(a5, a3 + a4);
        Assert.Equal(a6, a3 - a4);
        Assert.Equal(a7, a3 * a4);
        Assert.Equal(a8, a3 / a4);

        var f1 = new PadicZealous(5, 17, new Rational(313, 75)) * 625;
        var f2 = new PadicZealous(5, 17, 415) / 25;
        var f3 = f1 * 125 + f2 * 125;
        var f4 = (f1 + f2) * 125;
        Assert.Equal(f3, f4);
        Assert.NotEqual(f3.PNVSM, f4.PNVSM);
        Assert.True((f3 - f4).Norm.IsZero());
    }
}
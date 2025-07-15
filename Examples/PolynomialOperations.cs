using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.UserGroup.Integers;

namespace Examples;

public static class PolynomialOperations
{
    public static void Qi()
    {
        var (x, i) = Ring.Polynomial(Rational.KZero(), "X", "i").Deconstruct();
        var f = x.Pow(2) + 1;
        var fi = i.Pow(2) + 1;
        var p = (x - i) * (x + i);

        Console.WriteLine(p);
        Console.WriteLine(p.Div(f));
        Console.WriteLine(p.Div(fi).rem.Div(f));
    }

    public static void Qsqrt2()
    {
        var (x, a) = Ring.Polynomial(Rational.KZero(), "X", "a").Deconstruct();
        var f = x.Pow(2) - 2;
        var fa = a.Pow(2) - 2;
        var p = (x - a) * (x + a);

        Console.WriteLine(p);
        Console.WriteLine(p.Div(f));
        Console.WriteLine(p.Div(fa).rem.Div(f));
    }

    public static void Qcbrt1()
    {
        var (x, a) = Ring.Polynomial(Rational.KZero(), "X", "a").Deconstruct();
        var f = x.Pow(3) - 1;
        var fa = a.Pow(2) + a + 1;
        var p = (x - 1) * (x - a) * (x + a + 1);

        Console.WriteLine(p);
        Console.WriteLine(p.Div(f));
        Console.WriteLine(p.Div(fa).rem.Div(f));
    }

    public static void Qcbrt5()
    {
        var (x, a, j) = Ring.Polynomial(Rational.KZero(), "X", "a", "j").Deconstruct();
        var f = x.Pow(3) - 5;
        var fa = a.Pow(3) - 5;
        var fj = j.Pow(2) + j + 1;
        var p = (x - a) * (x - j * a) * (x + (j + 1) * a);

        Console.WriteLine(p);
        Console.WriteLine(p.Div(f));
        Console.WriteLine(p.Div(fa).rem.Div(f));
        Console.WriteLine(p.Div(fa).rem.Div(fj).rem.Div(f));
    }

    public static void Qsqrt6()
    {
        var (x, a, b) = Ring.Polynomial(Rational.KZero(), "X", "a", "b").Deconstruct();
        var f = x.Pow(2) - 6;
        var fa = a.Pow(2) - 2;
        var fb = b.Pow(2) - 3;
        var p = (x - a * b) * (x + a * b);

        Console.WriteLine(p);
        Console.WriteLine(p.Div(fa).rem.Div(fb));
        Console.WriteLine(p.Div(f));
        Console.WriteLine(p.Div(fa).rem.Div(f));
        Console.WriteLine(p.Div(fa).rem.Div(fb).rem.Div(f));
    }

    public static void Qisqrt2()
    {
        var (x, a, i) = Ring.Polynomial(Rational.KZero(), "X", "a", "i").Deconstruct();
        var f = x.Pow(4) - x.Pow(2) - 2;
        var fi = i.Pow(2) + 1;
        var fa = a.Pow(2) - 2;
        var p = (x - a) * (x + a) * (x - i) * (x + i);

        Console.WriteLine(p);
        Console.WriteLine(p.Div(f));
        Console.WriteLine("#");

        Console.WriteLine(p.Div(fi).rem.Div(f));
        Console.WriteLine(p.Div(fi).rem.Div(fa).rem);
        Console.WriteLine(p.Div(fi).rem.Div(fa).rem.Div(f));

        Console.WriteLine("#");
        Console.WriteLine(p.Div(fa).rem.Div(f));
        Console.WriteLine(p.Div(fa).rem.Div(fi).rem.Div(f));
    }

    public static void Qexpr()
    {
        var (x, a) = Ring.Polynomial(Rational.KZero(), "X", "a").Deconstruct();
        var f = x.Pow(2) - 2 * x - 2;
        var fa = a.Pow(2) - 3;
        var p = (x - (1 + a)) * (x - (1 - a));

        Console.WriteLine(p);
        Console.WriteLine(p.Div(fa).rem);
        Console.WriteLine(p.Div(fa).rem.Div(f));
    }
}
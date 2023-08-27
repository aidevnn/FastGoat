using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Integers;

namespace FastGoat.Examples;

public static class GroebnerBasis
{
    static GroebnerBasis()
    {
        Ring.DisplayPolynomial = MonomDisplay.StarCaret;
    }

    static Polynomial<Rational, Xi> Simplify(Polynomial<Rational, Xi> f)
    {
        var gcdNum = IntExt.GcdBigInt(f.Coefs.Where(c => !c.Value.IsZero()).Select(c => c.Value.Num).ToArray());
        var gcdDenom = IntExt.GcdBigInt(f.Coefs.Where(c => !c.Value.IsZero()).Select(c => c.Value.Denom).ToArray());
        var r = new Rational(gcdDenom, gcdNum);
        return f * r;
    }

    public static void Example1()
    {
        var (x, y) = Ring.Polynomial(Rational.KZero(), "x", "y").Deconstruct();
        foreach (var order in new[] { MonomOrder.Lex, MonomOrder.GrLex })
        {
            x.Indeterminates.SetOrder(order);
            var (p1, p2) = (x.Pow(2) * y - 1, x * y.Pow(2) - 3);
            Console.WriteLine($"GroebnerBasis[{p1}, {p2}]{p1.Indeterminates}");
            var bs0 = Ring.GroebnerBasis(p1, p2);
            Console.WriteLine(bs0.Select(Simplify).Glue("\n"));
            Console.WriteLine();
            Console.WriteLine($"ReducedGrobnerBasis[{p1}, {p2}]{p1.Indeterminates}");
            var bs1 = Ring.ReducedGrobnerBasis(p1, p2);
            Console.WriteLine(bs1.Glue("\n"));
            Console.WriteLine();
        }
    }

    public static void Example2()
    {
        var (x, y) = Ring.Polynomial(Rational.KZero(), "x", "y").Deconstruct();
        foreach (var order in new[] { MonomOrder.Lex, MonomOrder.GrLex })
        {
            x.Indeterminates.SetOrder(order);
            var (p1, p2) = (x.Pow(3) - 2 * x * y, x.Pow(2) * y - 2 * y.Pow(2) + x);
            Console.WriteLine($"GroebnerBasis[{p1}, {p2}]{p1.Indeterminates}");
            var bs0 = Ring.GroebnerBasis(p1, p2);
            Console.WriteLine(bs0.Select(Simplify).Glue("\n"));
            Console.WriteLine();
            Console.WriteLine($"ReducedGrobnerBasis[{p1}, {p2}]{p1.Indeterminates}");
            var bs1 = Ring.ReducedGrobnerBasis(p1, p2);
            Console.WriteLine(bs1.Glue("\n"));
            Console.WriteLine();
        }
    }

    public static void Example3()
    {
        var (x, y, z) = Ring.Polynomial(Rational.KZero(), "x", "y", "z").Deconstruct();
        foreach (var order in new[] { MonomOrder.Lex, MonomOrder.GrLex })
        {
            x.Indeterminates.SetOrder(order);
            var (p1, p2, p3) = (x + y - z, x.Pow(4) + x.Pow(3) + x.Pow(2) + x + 1, y.Pow(2) - 5);
            Console.WriteLine($"GroebnerBasis[{p1}, {p2}, {p3}]{p1.Indeterminates}");
            var bs0 = Ring.GroebnerBasis(p1, p2, p3);
            Console.WriteLine(bs0.Select(Simplify).Glue("\n"));
            Console.WriteLine();
            Console.WriteLine($"ReducedGrobnerBasis[{p1}, {p2}, {p3}]{p1.Indeterminates}");
            var bs1 = Ring.ReducedGrobnerBasis(p1, p2, p3);
            Console.WriteLine(bs1.Glue("\n"));
            Console.WriteLine();
        }
    }

    public static void Example4()
    {
        var (x, y, z, t) = Ring.Polynomial(Rational.KZero(), "x", "y", "z", "t").Deconstruct();
        foreach (var order in new[] { MonomOrder.Lex, MonomOrder.GrLex })
        {
            x.Indeterminates.SetOrder(order);
            var (p1, p2, p3) = (x + y - z, x * x - 2 * t * t, y * y - 5 * t * t);
            Console.WriteLine($"GroebnerBasis[{p1}, {p2}, {p3}]{p1.Indeterminates}");
            var bs0 = Ring.GroebnerBasis(p1, p2, p3);
            Console.WriteLine(bs0.Select(Simplify).Glue("\n"));
            Console.WriteLine();
            Console.WriteLine($"ReducedGrobnerBasis[{p1}, {p2}, {p3}]{p1.Indeterminates}");
            var bs1 = Ring.ReducedGrobnerBasis(p1, p2, p3);
            Console.WriteLine(bs1.Glue("\n"));
            Console.WriteLine();
        }
    }

    public static void LCMexample1()
    {
        var (t, x, y) = Ring.Polynomial(Rational.KZero(), MonomOrder.Lex, "t", "x", "y").Deconstruct();
        var (e1, e2) = ( x.Pow(2) * y, x * y.Pow(2));
        var (p1, p2) = (t * e1, (1 - t) * e2);
        Console.WriteLine($"ReducedGrobnerBasis[{p1}, {p2}]{p1.Indeterminates}");
        var bs1 = Ring.ReducedGrobnerBasis(p1, p2);
        bs1.Println();
        Console.WriteLine();

        var lcm = Ring.LcmPolynomial(e1, e2);
        var gcd = e1 * e2 / lcm;
        var (f1, f2) = (e1.Div(gcd), e2.Div(gcd));
        Console.WriteLine(new { e1, e2, lcm, gcd, f1, f2 });
    }

    public static void LCMexample2()
    {
        var (t, x, y) = Ring.Polynomial(Rational.KZero(), MonomOrder.Lex, "t", "x", "y").Deconstruct();

        var (e1, e2) = ((x + y).Pow(4) * (x.Pow(2) + y).Pow(2) * (x - 5 * y), (x + y) * (x.Pow(2) + y).Pow(3) * (x + 3 * y));
        var (p1, p2) = (t * e1, (1 - t) * e2);
        Console.WriteLine($"ReducedGrobnerBasis[{p1}, {p2}]{p1.Indeterminates}");
        var bs1 = Ring.ReducedGrobnerBasis(p1, p2);
        bs1.Println();
        Console.WriteLine();

        var lcm = Ring.LcmPolynomial(e1, e2);
        var gcd = e1 * e2 / lcm;
        var (f1, f2) = (e1.Div(gcd), e2.Div(gcd));
        var (g1, g2) = ((x + y).Pow(3) *  (x - 5 * y), (x.Pow(2) + y) * (x + 3 * y));
        Console.WriteLine(new { e1 });
        Console.WriteLine(new { e2 });
        Console.WriteLine(new { lcm });
        Console.WriteLine(new { gcd });
        Console.WriteLine(new { f1 });
        Console.WriteLine(new { g1 });
        Console.WriteLine(new { f2 });
        Console.WriteLine(new { g2 });
    }
}
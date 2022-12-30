using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Integers;

namespace FastGoat.Examples;

public static class GroebnerBasis
{
    static GroebnerBasis()
    {
        Monom.Display = MonomDisplay.StarCaret;
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
        var (x, y) = Ring.Polynomial("x", "y", Rational.KZero());
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
        var (x, y) = Ring.Polynomial("x", "y", Rational.KZero());
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
        var (x, y, z) = Ring.Polynomial("x", "y", "z", Rational.KZero());
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
        var (x, y, z, t) = Ring.Polynomial("x", "y", "z", "t", Rational.KZero());
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
}
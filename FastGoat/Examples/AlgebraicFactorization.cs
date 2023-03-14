using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Integers;
using static FastGoat.UserGroup.Polynoms.IntFactorisation;

namespace FastGoat.Examples;

public static class AlgebraicFactorization
{
    public static void MinimalPolynomials()
    {
        {
            var x = FG.QPoly('a');
            var f = x.Pow(2) - 5;
            var a = FG.EPoly(f);
            var X = FG.KPoly('X', a);
            CharacPoly(a);
            NormDetails(X - a);
            CharacPoly((a + 3) / 2);
            NormDetails(X - (a + 3) / 2);
        }

        {
            var x = FG.QPoly('a');
            var f = x.Pow(2) - 3 * x + 1;
            var a = FG.EPoly(f);
            var X = FG.KPoly('X', a);
            CharacPoly(a);
            NormDetails(X - a);
            CharacPoly(2 * a - 3);
            NormDetails(X - (2 * a - 3));
        }

        {
            var x = FG.QPoly('a');
            var f = x.Pow(4) + x.Pow(3) + x.Pow(2) + x + 1;
            var a = FG.EPoly(f);
            var X = FG.KPoly('X', a);
            CharacPoly(a);
            NormDetails(X - a);
            CharacPoly(a + a.Inv());
            NormDetails(X - (a + a.Inv()));
        }

        {
            var a = FG.EQPoly('a', -2, 0, 0, 1);
            var X = FG.KPoly('X', a);
            NormDetails(X.Pow(3) - (a - 1));
        }
    }

    public static void AlgFactorization()
    {
        Ring.DisplayPolynomial = MonomDisplay.Caret;
        
        {
            var i = FG.EQPoly('i', 1, 0, 1);
            var x = FG.KPoly('x', i);
            var A = x.Pow(4) - x.Pow(2) - 2;
            AlgebraicFactors(A, true);
            // x^4 + -x^2 + -2 = (x + -i) * (x + i) * (x^2 + -2)
        }

        {
            var a = FG.EQPoly('a', -1, -3, 0, 1);
            var x = FG.KPoly('x', a);
            var f = x.Pow(3) - 3 * x - 1;
            AlgebraicFactors(f, true);
            // x^3 + -3·x + -1 = (x + a^2 + -2) * (x + -a) * (x + -a^2 + a + 2)
        }

        {
            var a = FG.EQPoly('a', 1, 1, 1, 1, 1);
            var x = FG.KPoly('x', a);
            var A = x.Pow(2) - 5;
            AlgebraicFactors(A, true);
            // x^2 + -5 = (x + 2·a^3 + 2·a^2 + 1) * (x + -2·a^3 + -2·a^2 + -1)
        }

        {
            var a = FG.EQPoly('a', -5, 0, 1);
            var x = FG.KPoly('x', a);
            var A = x.Pow(4) + x.Pow(3) + x.Pow(2) + x + 1;
            AlgebraicFactors(A, true);
            // x^4 + x^3 + x^2 + x + 1 = (x^2 + (1/2·a + 1/2)·x + 1) * (x^2 + (-1/2·a + 1/2)·x + 1)
        }

        {
            var a = FG.EQPoly('a', -2, 0, 0, 1);
            var x = FG.KPoly('x', a);
            var A = x.Pow(3) + 3 * x.Pow(2) + 3 * x - 1;
            AlgebraicFactors(A, true);
            // x^3 + 3·x^2 + 3·x + -1 = (x + -a + 1) * (x^2 + (a + 2)·x + a^2 + a + 1)
        }

        {
            var i = FG.EQPoly('i', 1, 0, 1);
            var x = FG.KPoly('x', i);
            var A = x.Pow(4) + 25 * x.Pow(2) + 50 * x + 25;
            AlgebraicFactors(A, true);
            // x^4 + 25·x^2 + 50·x + 25 = (x^2 + -5·i·x + -5·i) * (x^2 + 5·i·x + 5·i)
        }

        {
            var a = FG.EQPoly('a', 2, 2, 1);
            var x = FG.KPoly('x', a);
            var A = x.Pow(4) + 4;
            AlgebraicFactors(A, true);
            // x^4 + 4 = (x + -a) * (x + a) * (x + a + 2) * (x + -a + -2)
        }

        {
            var a = FG.EQPoly('a', 3, 0, 1);
            var x = FG.KPoly('x', a);
            var A = x.Pow(6) - 1;
            AlgebraicFactors(A, true);
            // x^6 + -1 = (x + 1/2·a + 1/2) * (x + 1/2·a + -1/2) * (x + -1/2·a + 1/2) * (x + -1/2·a + -1/2) * (x + 1) * (x + -1)
        }
    }

    public static void PrimitiveEltExamples()
    {
        {
            var (b, _) = FG.EPolyXC(Rational.KZero(), 'a', 'b', -2, 0, 1);
            var f = b.Pow(2) - 3;
            var (r0, a0, b0) = PrimitiveElt(f);
            Console.WriteLine($"With {b[0].F} = 0 and {f} = 0");
            Console.WriteLine($"Primitive element y with Q(y) = Q(a,b)");
            Console.WriteLine(new[] { "Trager SqfrNorm", $"{r0} = 0", $"a = {a0}", $"b = {b0}" }.Glue("\n"));
            var (r1, a1, b1) = PrimEltGb(f);
            Console.WriteLine(new[] { "Grobner Basis", $"{r1} = 0", $"a = {a1}", $"b = {b1}" }.Glue("\n"));
            Console.WriteLine();
        }

        {
            var (b, _) = FG.EPolyXC(Rational.KZero(), 'a', 'b', 1, 1, 0, 1);
            var f = b.Pow(3) - b.Pow(2) + 4 * b - 3;
            var (r0, a0, b0) = PrimitiveElt(f);
            Console.WriteLine($"With {b[0].F} = 0 and {f} = 0");
            Console.WriteLine($"Primitive element y with Q(y) = Q(a,b)");
            Console.WriteLine(new[] { "Trager SqfrNorm", $"{r0} = 0", $"a = {a0}", $"b = {b0}" }.Glue("\n"));
            var (r1, a1, b1) = PrimEltGb(f);
            Console.WriteLine(new[] { "Grobner Basis", $"{r1} = 0", $"a = {a1}", $"b = {b1}" }.Glue("\n"));
            Console.WriteLine();
        }

        {
            var (b, _) = FG.EPolyXC(Rational.KZero(), 'a', 'b', 1, 1, 1, 1, 1);
            var f = b.Pow(2) - 5;
            var (r0, a0, b0) = PrimitiveElt(f);
            Console.WriteLine($"With {b[0].F} = 0 and {f} = 0");
            Console.WriteLine($"Primitive element y with Q(y) = Q(a,b)");
            Console.WriteLine(new[] { "Trager SqfrNorm", $"{r0} = 0", $"a = {a0}", $"b = {b0}" }.Glue("\n"));
            var (r1, a1, b1) = PrimEltGb(f);
            Console.WriteLine(new[] { "Grobner Basis", $"{r1} = 0", $"a = {a1}", $"b = {b1}" }.Glue("\n"));
            Console.WriteLine();
        }
    }

    public static void SplittingFieldQuarticPolynomial()
    {
        var x = FG.QPoly();

        SplittingField(x.Pow(2) - 3, true);
        SplittingField(x.Pow(2) - 3 * x - 3, true);
        SplittingField(x.Pow(2) + x + 1, true);
        SplittingField(x.Pow(2) + 2 * x - 5, true);

        SplittingField(x.Pow(3) - 2, true);
        SplittingField(x.Pow(3) - 3, true);
        SplittingField(x.Pow(3) - 3 * x - 1, true);
        SplittingField(x.Pow(3) - 2 * x + 2, true);
        SplittingField(x.Pow(3) + 2 * x.Pow(2) - x - 1, true);
        SplittingField(x.Pow(3) + 4 * x.Pow(2) + 3 * x + 1, true);
        SplittingField(x.Pow(3) - x + 1, true);

        SplittingField(x.Pow(4) + 4 * x.Pow(2) + 2, true);
        SplittingField(x.Pow(4) - 4 * x.Pow(2) + 2, true);
        SplittingField(x.Pow(4) + x.Pow(3) + x.Pow(2) + x + 1, true);
        SplittingField(x.Pow(4) - 2, true);
        SplittingField(x.Pow(4) + 2, true);
        SplittingField(x.Pow(4) + 5, true);
        SplittingField(x.Pow(4) + 3 * x.Pow(2) + 3, true);

        SplittingField(x.Pow(4) + x + 1, true); // Time:32489 ms
        SplittingField(x.Pow(4) + 3 * x.Pow(3) - x.Pow(2) + x + 1, true); // Time:46199 ms
    }

    // ∛(∛2 - 1) = ∛(1/9) + ∛(2/9) - ∛(4/9)
    public static void Ramanujan1()
    {
        var x = FG.QPoly();
        var (X, a0) = FG.EPolyXc(x.Pow(3) - 2, 'b');
        var (R, a1, _) = PrimitiveElt(X.Pow(3) - a0 + 1);
        var (X0, b) = FG.EPolyXc(R, 'b');
        var a = a1.Substitute(b);

        // Y^3 - 1/9 = 0 <=> Y = X/3 and X^3 - 3 = 0
        var c = AlgebraicRoots(X0.Pow(3) - 3, true).First() / 3;

        // Y^3 - 2/9 = 0 <=> Y = X/3 and X^3 - 6 = 0
        var d = AlgebraicRoots(X0.Pow(3) - 6, true).First() / 3;

        // Y^3 - 4/9 = 0 <=> Y = X/3 and X^3 - 4 = 0
        var e = AlgebraicRoots(X0.Pow(3) - 12, true).First() / 3;

        Console.WriteLine($"With {b.F} = 0");
        Console.WriteLine($"a - 1 = {a - 1} and a^3 = {a.Pow(3)}");
        Console.WriteLine($"c = {c} and c^3 = {c.Pow(3)}");
        Console.WriteLine($"d = {d} and d^3 = {d.Pow(3)}");
        Console.WriteLine($"e = {e} and e^3 = {e.Pow(3)}");
        Console.WriteLine($"2c^2 = {2 * c.Pow(2)} and d*e = {d * e}");
        Console.WriteLine($"c - d + e = {c - d + e}");
        
        // ∛(∛2 - 1) = ∛(1/9) + ∛(2/9) - ∛(4/9)
        /***
            With b⁹ + 3·b⁶ + 3·b³ + -1 = 0
            a - 1 = b³ and a^3 = 2
            c = 1/3·b⁴ + 2/3·b and c^3 = 1/9
            d = 1/3·b⁷ + b⁴ + 2/3·b and d^3 = 2/9
            e = 1/3·b⁷ + 2/3·b⁴ + b and e^3 = 4/9
            2c^2 = 2/9·b⁸ + 8/9·b⁵ + 8/9·b² and d*e = 2/9·b⁸ + 8/9·b⁵ + 8/9·b²
            c - d + e = b
         */
    }

    // √(∛5 - ∛4) = (∛2 + ∛20 - ∛25) / 3
    public static void Ramanujan2()
    {
        var x = FG.QPoly();
        var (X, _) = FG.EPolyXc(x.Pow(3) - 5, 'b');
        var (R, a0, b0) = PrimitiveElt(X.Pow(3) - 4);
        var (X0, y) = FG.EPolyXc(R, 'y');
        var a1 = a0.Substitute(y);
        var b1 = b0.Substitute(y);
     
        Console.WriteLine((R, a1, b1, a1.Pow(3), b1.Pow(3)));
        CharacPoly(a1 - b1);
        var R1 = Norm(X0 - a1 + b1);
        var (X1, c) = FG.EPolyXc(R1, 'c');
        var a2 = AlgebraicRoots(X1.Pow(3) - 5, true).First();
        var b2 = AlgebraicRoots(X1.Pow(3) - 4, true).First();

        Console.WriteLine($"With {c.F} = 0");
        Console.WriteLine($"a = {a2} and a^3 = {a2.Pow(3)}");
        Console.WriteLine($"b = {b2} and b^3 = {b2.Pow(3)}");

        Console.WriteLine(Norm(X1.Pow(2) - c));
        var d0 = AlgebraicRoots(X1.Pow(2) - c, true).Last();
        CharacPoly(d0);
        var R2 = Norm(X1 - d0);
        var (X2, d) = FG.EPolyXc(R2, 'd');
        var a = a2.Poly.Substitute(d.Pow(2));
        var b = b2.Poly.Substitute(d.Pow(2));
     
        var e = AlgebraicRoots(X2.Pow(3) - 2, true).First();
        var f = AlgebraicRoots(X2.Pow(3) - 20, true).First();
        var g = AlgebraicRoots(X2.Pow(3) - 25, true).First();
     
        Console.WriteLine($"With {d.F} = 0");
        Console.WriteLine($"a = {a} and a^3 = {a.Pow(3)}");
        Console.WriteLine($"b = {b} and b^3 = {b.Pow(3)}");
        Console.WriteLine($"a - b = {a - b}");
        Console.WriteLine();
     
        Console.WriteLine($"e = {e} and e^3 = {e.Pow(3)}");
        Console.WriteLine($"f = {f} and f^3 = {f.Pow(3)}");
        Console.WriteLine($"g = {g} and g^3 = {g.Pow(3)}");
        Console.WriteLine($"(e + f - g)/3 = {(e + f - g) / 3}");
        
        // √(∛5 - ∛4) = (∛2 + ∛20 - ∛25) / 3
        /***
            With d⁹ + 7·d⁶ + 23·d³ + -1 = 0
            a = 5/9·d⁸ + 4·d⁵ + 124/9·d² and a^3 = 5
            b = 5/9·d⁸ + 4·d⁵ + 115/9·d² and b^3 = 4
            a - b = d²

            e = 1/9·d⁷ + d⁴ + 32/9·d and e^3 = 2
            f = 1/3·d⁷ + 2·d⁴ + 23/3·d and f^3 = 20
            g = 4/9·d⁷ + 3·d⁴ + 74/9·d and g^3 = 25
            (e + f - g)/3 = d
         */
    }
}
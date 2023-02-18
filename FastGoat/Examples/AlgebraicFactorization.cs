using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Polynoms;

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
        Monom.Display = MonomDisplay.Caret;
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
            Console.WriteLine($"f={f}; A = f/(x-a) = {f.Div(x - a)}");
            var A = f / (x - a);
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
            AlgebraicFactors(A);
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
            Console.WriteLine($"Trager Primitive element y with Q(y) = Q(a,b)");
            Console.WriteLine(new[] { "SqfrNorm", $"{r0} = 0", $"a = {a0}", $"b = {b0}" }.Glue("\n"));
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

        // SplittingField(x.Pow(4) + x + 1, true); // Time:70s
        // SplittingField(x.Pow(4) + 3*x.Pow(3) - x.Pow(2) + x + 1, true); // Time:113s
    }
}
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Polynoms;
using static FastGoat.UserGroup.Polynoms.IntFactorisation;

namespace FastGoat.Examples;

public static class CardanFormula
{
    static CardanFormula()
    {
        Ring.DisplayPolynomial = MonomDisplay.StarCaret;
    }

    static (KPoly<Rational> p, Rational t) TransformPoly(KPoly<Rational> P)
    {
        var deg = P.Degree;
        if (!P[deg].Equals(P.KOne))
            throw new("Polynom isnt monic");

        if (deg == 0)
            return (P, P.KZero);

        var t = P[deg - 1] / deg;
        return (P.Substitute(P.X - t), t);
    }

    static void Cubic(KPoly<Rational> P)
    {
        if (P.Degree != 3 || !P.LT.Equals(Rational.KOne()))
            throw new("P must be monic and cubic");

        var (P1, c) = ConstCoef(P, monic: true);
        var (P2, t) = TransformPoly(P1);
        var (p, q) = (P2[1], P2[0]);
        var D = -4 * p.Pow(3) - 27 * q.Pow(2);
        var (numD, _) = D.Decomp();
        var numD0 = numD.Where(e => e.Item2 % 2 != 0).ToArray();
        var D0 = numD0.Select(e => new Rational(e.Item1)).Aggregate(Rational.KOne(), (acc, e) => acc * e);
        var x = P.X;

        var (y1, y2, y3) = SplittingField(P, true).Deconstruct();
        var X = FG.KPoly('X', y1);
        var (R, a0, b0) = PrimitiveElt(X.Pow(2) + X + 1);
        var (X1, x1) = FG.EPolyXc(R, 'y');
        var j = b0.Substitute(x1);
        y1 = y1.Substitute(x1).Substitute(a0);
        y2 = y2.Substitute(x1).Substitute(a0);
        y3 = y3.Substitute(x1).Substitute(a0);
        Console.WriteLine((X1 - y1) * (X1 - y2) * (X1 - y3));
        Console.WriteLine((j * j + j, j));
        var (u, v) = (y1 + j * y2 + j * j * y3, y1 + j * j * y2 + j * y3);
        Console.WriteLine((X1 - u.Pow(3)) * (X1 - v.Pow(3)));
        new[] { u, v, u.Pow(3), v.Pow(3) }.Println();
        AlgebraicFactors(X1.Pow(2) + 27 * q * X1.KOne * X1 - 27 * p.Pow(3) * X1.KOne).Println();
        var sqrtD = AlgebraicRoots(X1.Pow(2) - D * X1.KOne)[0];
        var sqrtm3 = 2 * j + 1;
        Console.WriteLine((sqrtD.Pow(2), sqrtm3.Pow(2)));
        Console.WriteLine(x1.F.SubstituteChar('X'));
        Console.WriteLine(ConstCoef(x1.F, monic: true));
        new[] { y1, y2, y3 }.Println("Roots");
        Console.WriteLine();
    }

    static void Cubic2(KPoly<Rational> P)
    {
        if (P.Degree != 3 || !P.LT.Equals(Rational.KOne()))
            throw new("P must be monic and cubic");

        var (P1, c) = ConstCoef(P, monic: true);
        var (P2, t) = TransformPoly(P1);
        var (p, q) = (P2[1], P2[0]);
        var D = -4 * p.Pow(3) - 27 * q.Pow(2);
        var x = P.X;

        var (u3, v3) = FactorsQuadratic(x.Pow(2) + 27 * q * x - 27 * p.Pow(3), false).roots.Select(e => -e[0]).Deconstruct();
        var D0 = u3.X.Pow(2)[0];

        var r3 = AlgebraicRoots(x.Pow(2) + 3);
        var e3 = r3[1];
        var X3 = FG.KPoly('X', e3);

        var (F, X, sqrtm3, sqrtD0) =
            !D0.Equals(-3 * D0.One) ? PrimitiveElt(x.Pow(2) + 3, x.Pow(2) - D0)[0] : (e3.F, X3, e3, e3);

        u3 = u3.Substitute(sqrtD0);
        v3 = v3.Substitute(sqrtD0);

        var (R, a0, b0) = PrimitiveElt(X.Pow(3) - u3);
        var (X1, x1) = FG.EPolyXc(R, 'y');
        var u = b0.Substitute(x1);
        var v3n = v3.Substitute(x1).Substitute(a0);
        var vr = AlgebraicRoots(X1.Pow(3) - v3n);
        var j = (-1 + sqrtm3.Substitute(x1).Substitute(a0)) / 2;
        var t0 = -t + x1.Zero;
        var v = vr.First(v0 =>
        {
            var (vy1, vy2, vy3) = ((u + v0) / 3 + t0, j * (u + j * v0) / 3 + t0, j * (j * u + v0) / 3 + t0);
            return ((X1 - vy1) * (X1 - vy2) * (X1 - vy3)).Equals(P.Substitute(X1));
        });

        var (y1, y2, y3) = ((u + v) / 3 + t0, j * (u + j * v) / 3 + t0, j * (j * u + v) / 3 + t0);
        Console.WriteLine("P = {0}", (X1 - y1) * (X1 - y2) * (X1 - y3));
        Console.WriteLine("Splitting field Q(a, j) minPoly : {0}", R.SubstituteChar('X'));
        new[] { y1, y2, y3 }.Order().Println("Roots");
        Console.WriteLine();

        var Q = P.Substitute(X1);
        AlgebraicRoots(Q).Order().Println();
        Console.WriteLine();
    }

    static void FactorsQDetails(KPoly<Rational> P)
    {
        var fact = FactorsQ(P);
        Console.WriteLine($"P = {P}");
        var c = 'a';
        foreach (var p in fact)
        {
            Console.WriteLine($"    {p}");
            if (p.Item1.Degree != 2)
                continue;

            FactorsQuadratic(p.Item1, a: c);
            c++;
        }

        Console.WriteLine();
    }

    public static void ExamplesFactorsQ()
    {
        var x = FG.QPoly();
        FactorsQDetails(x.Pow(2) + 2);
        FactorsQDetails(x.Pow(2) - new Rational(4, 9));
        FactorsQDetails(x.Pow(2) - new Rational(5, 9));
        FactorsQDetails(x.Pow(2) + new Rational(5, 8));
        FactorsQDetails(x.Pow(2) + x + 1);
        FactorsQDetails((3 * x + 1) * (x + 2));
        FactorsQDetails((3 * x + 1) * (x / 2 + 2));
        FactorsQDetails((3 * x + 1) * (x / 2 + 2).Pow(2));
        FactorsQDetails(x.Pow(2) + 2 * x - 1);
        FactorsQDetails(x.Pow(2) + 2 * x + 1);
        FactorsQDetails(x.Pow(2) / 3 - 2 * x + 3);
        FactorsQDetails(x.Pow(2) + x / 5 + x.KOne / 3);
        FactorsQDetails(x.Pow(2) + x / 3 - x.KOne / 2);
        FactorsQDetails(2 * x.Pow(2) + x / 3 - 1);
        FactorsQDetails(6 * x.Pow(2) + x - 3);

        FactorsQDetails((x.Pow(2) - 4).Pow(2) * (x + 1).Pow(3));
        FactorsQDetails((x + 3).Pow(2) / 4);

        FactorsQDetails((2 * x.Pow(2) + x - 3) * (5 * x.Pow(2) - 3 * x + 1) * (x.Pow(2) + x + 1));
    }

    public static void ExamplesFactorizationQ()
    {
        var X = FG.QPoly('X');

        Console.WriteLine("Random Z[X] Polynomials");
        Console.WriteLine();
        IntExt.RngSeed(1235);

        for (int j = 0; j < 10; j++)
        {
            var amp = IntExt.Rng.Next(2, 9);
            var n = 2 + IntExt.Rng.Next(13);
            var degrees = IntExt.Partitions32[n].Where(l => l.All(i => i != 1) && l.Count > 1)
                .OrderBy(_ => IntExt.Rng.NextDouble())
                .FirstOrDefault(new[] { 2, 3, 4 }.ToList())
                .ToArray();

            var rg = new[] { 1, 1, 1, 2, 2, 3 };
            var polys = degrees.Select(ni =>
                    PolynomialFactorization.RandPoly(Rational.KZero(), amp, ni, false).Pow(rg[IntExt.Rng.Next(6)]))
                .ToArray();
            var f0 = polys.Aggregate((a, b) => a * b);
            var f = new KPoly<Rational>('X', f0.KZero, f0.Coefs);

            Console.WriteLine($"Random factors : [{polys.Glue(" ;  ")}]");
            Console.WriteLine($"f0 = {f}");
            Console.WriteLine();

            var facts = FactorsQ(f, details: true);
            var f1 = facts.Select(e => e.Item1.Pow(e.Item2)).Aggregate((a0, a1) => a0 * a1).SubstituteChar(f0.x);
            if (!f0.Equals(f1))
                throw new($"f1 = {f1} f0 = {f0} => {f0.Equals(f1)}");
        }
    }

    public static void Examples1()
    {
        var x = FG.QPoly();

        Cubic(x.Pow(3) - 3 * x - 1);
        Cubic(x.Pow(3) + x.Pow(2) - 2 * x - 1);
        Cubic(x.Pow(3) - 3 * x - 3);
        Cubic(x.Pow(3) + x.Pow(2) - 2 * x + 2);
    }

    public static void Examples2()
    {
        var x = FG.QPoly();

        Cubic2(x.Pow(3) - 3 * x - 1);
        Cubic2(x.Pow(3) + x.Pow(2) - 2 * x - 1);
        Cubic2(x.Pow(3) - 3 * x - 3);
        Cubic2(x.Pow(3) + x.Pow(2) - 2 * x + 2);
    }
}
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Padic;
using FastGoat.UserGroup.Polynoms;

namespace FastGoat.Examples;

public static class BivariatePolynomialFactorization
{
    static BivariatePolynomialFactorization()
    {
        Ring.DisplayPolynomial = MonomDisplay.StarCaret;
        GlobalStopWatch.Restart();
    }

    /// <summary>
    /// Truncates a bivariate polynomial by restricting its degree between specified lower and upper bounds.
    /// </summary>
    /// <param name="P">The bivariate polynomial to be truncated.</param>
    /// <param name="t">The degree lower bound.</param>
    /// <param name="o">The degree upper bound.</param>
    /// <returns>The truncated bivariate polynomial.</returns>
    static KPoly<ZnBInt> Truncate(KPoly<ZnBInt> P, int t, int o)
    {
        var mod = P.KZero.Details;
        var z = ZnBInt.ZnZero(mod.P, mod.O);
        var coefs = o.Range().Select(i => i < t ? z : P[i]).TrimSeq().ToArray();
        return new KPoly<ZnBInt>(P.x, z, coefs);
    }

    /// <summary>
    /// The BivariatePolynomialFactorization class provides methods for factorizing bivariate polynomials over a finite field.
    /// </summary>
    public static void Example1_Fp()
    {
        Ring.DisplayPolynomial = MonomDisplay.StarCaret;
        var (X2, X1) = Ring.Polynomial(ZnBInt.ZnZero(101), MonomOrder.Lex, "X2", "X1").Deconstruct();

        // (X2^2 + 100*X1^2 + 100) * (X2^2 +  99*X2 + X1^2)
        var P = (X2.Pow(2) + 100 * X1.Pow(2) + 100) * (X2.Pow(2) + 99 * X2 + X1.Pow(2));
        var k = IntExt.Solve_k_pow_m_equal_one_mod_n_strict(101, 100);
        var P0y = P.Substitute(P.Zero, X1).ToKPoly(X2);
        var firr = IntFactorisation.Firr(P0y, P.KOne * k)
            .Select(f => f.ToPolynomial(P.Indeterminates, X2.ExtractIndeterminate))
            .ToArray();
        var lifts = IntFactorisation.HenselLifting(P, firr);
        var facts = IntFactorisation.Recombinaison(P, lifts);

        var P0X2 = firr.Aggregate(P.One, (acc, c) => acc * c);
        Console.WriteLine($"P(X1,X2) = {P}");
        firr.Println($"P(0,X2) = {P0X2}");
        lifts.Println("Hensel Lifting");
        facts.Println("Factors");

        var P0 = facts.Aggregate(P.One, (acc, c) => acc * c);

        Console.WriteLine($"{facts.Glue(" * ", "({0})")} = {P0}");
        Console.WriteLine($"Check:{P.Equals(P0)}");
        Console.WriteLine();
    }

    /// <summary>
    /// This method demonstrates an example of rational polynomial factorization using the FastGoat library.
    /// It computes the factorization of a given bivariate polynomial over a field of rational numbers.
    /// The example starts by setting the display format for polynomials. Then it defines a list of bivariate polynomials.
    /// Each polynomial in the list is created using operations like addition, subtraction, multiplication, and exponentiation
    /// on the monomials and coefficients of the variables X2 and X1.
    /// Finally, it iterates over the list of polynomials and calls the FactorsFxy method to perform the factorization.
    /// </summary>
    public static void Example2_Rational()
    {
        Ring.DisplayPolynomial = MonomDisplay.StarCaret;
        var (X2, X1) = Ring.Polynomial(Rational.KZero(), MonomOrder.Lex, "X2", "X1").Deconstruct();

        // All Working 
        var polys = new[]
        {
            (X2.Pow(2) - 3 * X2 - 4 * X1) * (X2.Pow(2) - 4 * X1 * X2 - 4),
            (X2 - 3 * X1.Pow(2) - 2 * X1) * (X2.Pow(2) - 3 * X1 * X2 - 4),
            (X2.Pow(2) - 3 * X1 - 3) * (X2.Pow(2) - 4 * X2 - 3 * X1.Pow(2)),
            (X2 - 2 * X1 - 3) * (X2.Pow(2) - X2 - 2 * X1),
            (X2 - 2 * X1 - 4) * (X2.Pow(2) - X1.Pow(2) - 3),
            (X2 - 4 * X1.Pow(2) - 2 * X1) * (X2.Pow(2) - 4 * X1.Pow(2) - 3)
        };

        Logger.Level = LogLevel.Level1;
        foreach (var F in polys)
        {
            Console.WriteLine($"F = {F}");
            IntFactorisation.FactorsFxy(F).Println($"Factors in Q[{X1},{X2}]");
            Console.WriteLine();
        }
    }

    /// <summary>
    /// The purpose of this method is to demonstrate batch factorization of random polynomials.
    /// The method uses a nested foreach loop to iterate over pairs of values for `n` and `m`.
    /// For each pair, the method calls the `FactorisationRandPolFxy` method with the specified parameters.
    /// After the loop finishes, the method calls the `Console.Beep` method and then displays the elapsed time
    /// using the `GlobalStopWatch.Show` method.
    /// </summary>
    public static void Example3_BatchRandomPolynomial()
    {
        IntExt.RngSeed(25413);
        Logger.Level = LogLevel.Level1;
        Ring.DisplayPolynomial = MonomDisplay.StarCaret;

        GlobalStopWatch.AddLap();
        foreach (var (n, m) in 4.Range(1).SelectMany(m => (5 - m).Range(2).Select(n => (n, m))))
        {
            var pols = IntFactorisation.RandRationalPolynomials(amplitude: 10, maxMonomDegree: m,
                    nbFacts: n, nbRandPolys: 100)
                .Where(F => F.CoefMax("Y").Degree == 0 && IntFactorisation.FilterRandPolynomialFxy(F))
                .Distinct().Take(10).ToArray();

            foreach (var (idx, F) in pols.Index())
            {
                Console.WriteLine($"P{idx,-2} = {F}");
                var (y, x) = F.Indeterminates.Deconstruct();
                var factsF = IntFactorisation.FactorsFxy(F, rewrite: true)
                    .Select(f => IntFactorisation.Primitive(f)).Where(f => !(f - 1).IsZero()).Order().ToArray();
                var lt = F / factsF.Aggregate((acc, e) => e * acc);
                if (!(lt - 1).IsZero())
                    factsF = factsF.Prepend(lt).ToArray();

                factsF.Println($"Factors in Q[{x},{y}]");

                var prod = factsF.Aggregate((ai, aj) => ai * aj);
                if (!prod.Equals(F))
                    throw new($"Prod = {prod}");

                Console.WriteLine();
            }
        }

        Console.Beep();
        GlobalStopWatch.Show();
    }

    // AECF p405
    public static void Example4_Fp()
    {
        Ring.DisplayPolynomial = MonomDisplay.StarCaret;
        var (X2, X1) = Ring.Polynomial(ZnBInt.ZnZero(101), MonomOrder.Lex, "X2", "X1").Deconstruct();

        // (X2^2 + 100*X1^2 + 100) * (X2^2 +  99*X2 + X1^2)
        var P = (X2.Pow(2) + 100 * X1.Pow(2) + 100) * (X2.Pow(2) + 99 * X2 + X1.Pow(2));
        var k = IntExt.Solve_k_pow_m_equal_one_mod_n_strict(101, 100);
        var P0y = P.Substitute(P.Zero, X1).ToKPoly(X2);
        var firr = IntFactorisation.Firr(P0y, P.KOne * k)
            .Select(f => f.ToPolynomial(P.Indeterminates, X2.ExtractIndeterminate))
            .ToArray();

        var o = 10;
        var I = X1.Pow(o);
        var lifts = IntFactorisation.HenselLifting(P, firr, o);

        var P0X2 = firr.Aggregate(P.One, (acc, c) => acc * c);
        Console.WriteLine($"P(X1,X2) = {P}");
        firr.Println($"P(0,X2) = {P0X2}");
        lifts.Println("Hensel Lifting");
        Console.WriteLine();
        var nu = P.DegreeOf(X1.ExtractIndeterminate);
        var list = new List<KPoly<ZnBInt>[]>();
        var nbLifts = lifts.Length;
        foreach (var f in lifts)
        {
            var _p = P.Div(f).quo;
            var df = f.D(X2);
            var _f = (_p * df).Div(I).rem.Monic();
            var decomp = Ring.Decompose(_f, X2.ExtractIndeterminate).Item1;
            Console.WriteLine($"f = {f}");
            decomp.Println($"_f = {_f}");
            var col = decomp.Select(e => Truncate(e.Value.ToKPoly(X1), nu + 1, o - 1)).ToArray();
            var gcd = col.Aggregate(Ring.FastGCD).Monic;
            list.Add(col.Select(e => e.Div(gcd).quo).ToArray());
        }

        var maxDeg = list.Max(l => l.Max(p => p.Degree));
        var p0 = IntExt.Primes10000.First(p0 => p0 > maxDeg + 1); // p0 = 7
        var a = FG.EPoly(FG.CyclotomicPolynomial(p0).Substitute(FG.ZbPoly(101, 1)));
        var A = KMatrix<EPoly<ZnBInt>>.MergeSameRows(list.Select(l => l.Select(e => e.Substitute(a)).ToKMatrix(nbLifts))
            .ToArray());
        Console.WriteLine("A");
        Console.WriteLine(A);
        Console.WriteLine("Ker(A)");
        var ker = A.NullSpace().Item2;
        Console.WriteLine(ker);

        ker.Cols.Select(col => col.Select((c, i) => c[0] * lifts[i]).Where(f => !f.IsZero())
            .Aggregate((acc, e) => e * acc).Div(X1.Pow(nu + 1)).rem).Println($"Factors in F101[X1,X2]");
    }

    /// <summary>
    /// Contains example related to the substitution of variables using method Rewrite in a bivariate polynomial.
    /// </summary>
    public static void Example5_Substitution()
    {
        var (X2, X1) = Ring.Polynomial(Rational.KZero(), MonomOrder.Lex, "X2", "X1").Deconstruct();
        var P = (X2 * X1 + 1) * (X2 + X1 - 3);
        var (_, i, F) = IntFactorisation.RewritingPolynomialLeadingTerm(P);
        Console.WriteLine($"P = {P}");
        Console.WriteLine($"Substitute {X1} <- {X1 + i * X2}");
        Console.WriteLine($"F = {F}");
        var factsF0 = IntFactorisation.FactorsFxy(F);
        factsF0.Println($"Factors F in Q[{X1},{X2}]");
        var factsF = factsF0.Select(fi => fi.Substitute(X1 - i * X2, X1)).Order().ToArray();

        Console.WriteLine($"Substitute {X1} <- {X1 - i * X2}");
        factsF.Println($"Factors P in Q[{X1},{X2}]");
        Console.WriteLine();
    }

    public static void Example6_Recombinaison_case()
    {
        var (X2, X1) = Ring.Polynomial(Rational.KZero(), MonomOrder.Lex, "X2", "X1").Deconstruct();
        var F = 5 * X1.Pow(2) * X2.Pow(4) + 9 * X1.Pow(3) * X2.Pow(3) + 14 * X1.Pow(2) * X2.Pow(3) +
            4 * X1.Pow(4) * X2.Pow(2) + 6 * X1.Pow(3) * X2.Pow(2) - 3 * X1.Pow(2) * X2.Pow(2) - 18 * X1 * X2.Pow(2) -
            4 * X1.Pow(4) * X2 + X1.Pow(3) * X2 - 14 * X1.Pow(2) * X2 + 10 * X1 * X2 - 2 * X1.Pow(2) - 8;

        var (_, i0, F0) = IntFactorisation.RewritingPolynomialResultantZero(F);
        Console.WriteLine(new { F });
        Console.WriteLine($"Res(F)(0,0) != 0 : {IntFactorisation.FilterRandPolynomialFxy(F)}");
        Console.WriteLine(new { F0 });
        Console.WriteLine($"Substitute {X1} <- {X1 + i0}");
        Console.WriteLine($"Res(F0)(0,0) != 0 : {IntFactorisation.FilterRandPolynomialFxy(F0)}");
        var (_, i1, F1) = IntFactorisation.RewritingPolynomialLeadingTerm(F0);
        Console.WriteLine(new { F1 });
        Console.WriteLine($"Substitute {X1} <- {X1 + i1 * X2}");
        var facts = IntFactorisation.FactorsFxy(F1);
        var factsF1 = facts.Select(f => f.Substitute(X1 - i0 * X2, X1)).Order().ToArray();
        var factsF = factsF1.Select(f => f.Substitute(X1 - i1, X1))
            .Select(f => IntFactorisation.Primitive(f)).Order().ToArray();
        var prod = factsF.Aggregate((acc, e) => e * acc);
        var lc = F.LeadingDetails.lc / prod.LeadingDetails.lc;
        if (!(lc - 1).IsZero())
            factsF = factsF.Prepend(lc * F.One).ToArray();

        prod = factsF.Aggregate((acc, e) => e * acc);
        Console.WriteLine($"F = {F}");
        factsF.Println($"Factors in Q[{X1},{X2}]");
        Console.WriteLine();
        Console.WriteLine($"Check:{F.Equals(prod)}");
    }

    public static void Example7_FullSubstitutions()
    {
        // IntExt.RngSeed(812665);
        GlobalStopWatch.AddLap();
        Logger.Level = LogLevel.Level1;
        foreach (var (n, m) in 4.Range(1).SelectMany(m => (5 - m).Range(2).Select(n => (n, m))))
        {
            var ct = 0;
            GlobalStopWatch.AddLap();
            while (ct < 5)
            {
                var facts = IntFactorisation.RandRationalPolynomials(10, maxMonomDegree: m, nbMonoms: 6, nbFacts: 1)
                    .Distinct().Take(n).ToArray();
                var F = facts.Aggregate((ai, aj) => ai * aj);
                var (y, x) = F.Indeterminates.Deconstruct();
                if (F.CoefMax(y).Degree == 0 && m > 1)
                    continue;

                if (IntFactorisation.FilterRandPolynomialFxy(F))
                    continue;

                Console.WriteLine($"P{ct} = {F}");
                var factsF = IntFactorisation.FactorsFxy(F, true);
                factsF.Println($"Factors in Q[{x},{y}]");

                var prod = factsF.Aggregate((ai, aj) => ai * aj);
                if (!prod.Equals(F))
                    throw new($"Prod = {prod}");
                Console.WriteLine();
                ++ct;
            }

            GlobalStopWatch.Show($"MaxDegree:{m} NbFactors:{n}");
            Console.WriteLine();
        }

        GlobalStopWatch.Show("End");
        Console.Beep();
    }

    public static void Example8_NonSeparable()
    {
        GlobalStopWatch.AddLap();
        Ring.DisplayPolynomial = MonomDisplay.StarCaret;
        var (X2, X1) = Ring.Polynomial(Rational.KZero(), MonomOrder.Lex, "X2", "X1").Deconstruct();
        var F = X2.Pow(6) - 9 * X1 * X2.Pow(5) + 25 * X1.Pow(2) * X2.Pow(4) - 3 * X1 * X2.Pow(4) -
            44 * X1.Pow(3) * X2.Pow(3) + 27 * X1.Pow(2) * X2.Pow(3) + 72 * X1.Pow(4) * X2.Pow(2) -
            63 * X1.Pow(3) * X2.Pow(2) - 32 * X1.Pow(5) * X2 + 24 * X1.Pow(4) * X2 - 48 * X1.Pow(6) + 36 * X1.Pow(5);

        Console.WriteLine($"F = {F}");
        GlobalStopWatch.AddLap();
        var sff = IntFactorisation.SFF(F);
        sff.Println("SFF");
        GlobalStopWatch.Show("SFF");

        var facts = sff.SelectMany(f => IntFactorisation.FactorsFxy(f.g, rewrite: true).Select(h => (g: h, i: f.m)))
            .ToArray();

        Console.WriteLine();
        Console.WriteLine($"F = {F}");
        facts.Println("Factors with multiplicity");

        var prod = facts.Aggregate(F.One, (acc, e) => acc * e.g.Pow(e.i));
        var check = prod.Equals(F);
        Console.WriteLine(new { check });
        GlobalStopWatch.Show();
        Logger.Level = LogLevel.Level1;
        IntFactorisation.FactorsFxy(F, rewrite: true).Println("Fail");
    }

    public static void Example9_FullFactorisationRational()
    {
        // IntExt.RngSeed(812665);
        GlobalStopWatch.AddLap();
        // Logger.Level = LogLevel.Level1;
        foreach (var (n, m) in 4.Range(1).SelectMany(m => (6 - m).Range(2).Select(n => (n, m))))
        {
            var ct = 0;
            GlobalStopWatch.AddLap();
            while (ct < 5)
            {
                var F = IntFactorisation.RandRationalPolynomials(10, maxMonomDegree: m, nbMonoms: 4, nbFacts: 1)
                    .SelectMany(f => Enumerable.Repeat(f, IntExt.Rng.Next(1, 4))).Take(3 * n).Shuffle().Take(n + 1)
                    .Aggregate((ai, aj) => ai * aj);
                var (y, x) = F.Indeterminates.Deconstruct();
                if (F.DegreeOf(y) == 0 || F.Coefs.Max(e => e.Value.Absolute.Num) > 5000)
                    continue;

                Console.WriteLine("#######################################");
                Console.WriteLine($"P{ct} = {F}");
                var (F1, ctx, cty) = IntFactorisation.CT(F);
                var sff = IntFactorisation.SFF(F1).Select(f => (g: IntFactorisation.Primitive(f.g), f.m)).ToArray();

                if (sff.Sum(f => f.m) > 1)
                {
                    if (!F.Equals(F1))
                    {
                        Console.WriteLine($"CTX = {ctx} CTY = {cty}");
                        Console.WriteLine($"F1 = {F1}");
                    }

                    sff.Println("SFF");
                }

                var ctXY = new List<(Polynomial<Rational, Xi> g, int m)>();
                if (!ctx.IsOne())
                    ctXY.Add((ctx, 1));

                if (!cty.IsOne())
                    ctXY.Add((F.X(y), cty.Degree));

                var factsF = sff.SelectMany(f => IntFactorisation.FactorsFxy(f.g, rewrite: true).Select(g => (g, f.m)))
                    .Concat(ctXY)
                    .Select(e => (g: IntFactorisation.Primitive(e.g), e.m))
                    .GroupBy(e => e.g).Select(e => (g: e.Key, m: e.Sum(f => f.Item2)))
                    .OrderBy(e => e.g.NbIndeterminates)
                    .ThenBy(e => e.g.Coefs.Keys.Count)
                    .ThenBy(e => e.m)
                    .ThenBy(e => e.g)
                    .ToArray();
                var prod = factsF.Aggregate(F.One, (acc, e) => e.Item1.Pow(e.Item2) * acc);
                var lc = F.Div(prod).quo;
                if (!lc.Equals(lc.One))
                    factsF = factsF.Prepend((lc, 1)).ToArray();

                Console.WriteLine();
                factsF.Println($"Factors in Q[{x},{y}]");
                Console.WriteLine();

                Console.WriteLine("## sage export");
                Console.WriteLine($"k.<X,Y>=QQ[]; factor({F})");
                prod = factsF.Aggregate(F.One, (acc, e) => e.Item1.Pow(e.Item2) * acc);
                if (!prod.Equals(F) || factsF.Where(e => e.g.Degree > 0).Sum(e => e.m) < n)
                    throw new();

                Console.WriteLine();
                ++ct;
            }

            GlobalStopWatch.Show($"MaxDegree:{m} NbFactors:{n}");
            Console.WriteLine();
        }

        GlobalStopWatch.Show("End");
        Console.Beep();
    }

    public static void Example10_FullFactorisationFp()
    {
        foreach (var (p, n) in new[] { (2, 4), (3, 4), (5, 3), (11, 3) })
        {
            var (Y, X) = Ring.Polynomial(ZnInt.ZnZero(p), MonomOrder.Lex, "Y", "X").Deconstruct();
            var (x, y) = (X.ExtractIndeterminate, Y.ExtractIndeterminate);

            foreach (var (i0, P0) in IntFactorisation.RandZnIntPolynomials(p, nbFacts: n).Distinct().Take(5).Index())
            {
                var P1 = new Polynomial<ZnInt, Xi>(P0.Indeterminates, P0.KZero,
                    new(P0.Coefs.ToDictionary(e => e.Key, e => e.Value)));
                var i = i0 + 1;
                Console.WriteLine($"P{i} = {P1}");
                var facts = IntFactorisation.FactorFxyFp(P1);

                facts.OrderBy(e => e.g.NbIndeterminates)
                    .ThenBy(e => e.g.Coefs.Keys.Count)
                    .ThenBy(e => e.m)
                    .ThenBy(e => e.g)
                    .Select(e => (e.g.ToKPoly(y, x), e.m))
                    .Println($"Factors P{i} = {P1}");

                Console.WriteLine();
                var prod = facts.Aggregate(P1.One, (acc, fi) => acc * fi.g.Pow(fi.m));
                Console.WriteLine("## sage export");
                Console.WriteLine($"k.<X,Y>=GF({p})[]; factor({P1})");
                if (!P1.Equals(prod) || facts.Sum(e => e.m) < n)
                    throw new($"{P1.Div(prod)}");

                Console.WriteLine();
            }
        }
    }
}
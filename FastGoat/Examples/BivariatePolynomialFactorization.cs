using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Polynoms;

namespace FastGoat.Examples;

public static class BivariatePolynomialFactorization
{
    static BivariatePolynomialFactorization()
    {
        Ring.DisplayPolynomial = MonomDisplay.StarCaret;
    }
    
    /// <summary>
    /// Computes the factorization of a given bivariate polynomial over a finite field.
    /// First, evaluates the bivariate polynomial at x = 0 to obtain an univariate polynomial in y,
    /// then factorizes the resulting univariate polynomial in y using Berlekamp's algorithm.
    /// </summary>
    /// <param name="F">The bivariate polynomial to be factorized.</param>
    /// <returns>An array of polynomials representing the factors of the input polynomial.</returns>
    /// <exception cref="ArgumentException">Thrown when the characteristic of the finite field is zero.</exception>
    static Polynomial<ZnInt, Xi>[] Firr(Polynomial<ZnInt, Xi> F)
    {
        if (F.P == 0)
            throw new();

        var p = F.P;
        var (y, x) = F.Indeterminates.Deconstruct();
        var P0y = F.Substitute(F.Zero, x);
        var _P0y = P0y.ToKPoly(y);
        var k = IntExt.SolveAll_k_pow_m_equal_one_mod_n_strict(p, p - 1).First();
        var a0 = k * _P0y.KOne;
        var firr = IntFactorisation.BerlekampProbabilisticAECF(_P0y, a0).Order().ToArray();
        return firr.Select(fi => fi.Substitute(F.X(y))).ToArray();
    }

    /// <summary>
    /// Performs one step of the Hensel Lifting method to lift the factorization of a bivariate polynomial 
    /// from one modulus to a higher power of the modulus.
    /// </summary>
    /// <param name="F">The original bivariate polynomial to be factorized.</param>
    /// <param name="fi">The array of current factors of the polynomial F.</param>
    /// <param name="I">The ideal used in the lifting process, which is a power of 2 of the variable "x"..</param>
    /// <param name="o">The order of the lifting step.</param>
    /// <returns>An array of polynomials representing the factors of the input polynomial after one lifting step.</returns>
    /// <exception cref="ArgumentException">Thrown when the polynomial does not have exactly two indeterminates.</exception>
    static Polynomial<ZnInt, Xi>[] HenselLiftingStep(Polynomial<ZnInt, Xi> F, Polynomial<ZnInt, Xi>[] fi,
        Polynomial<ZnInt, Xi> I, int o)
    {
        var xis = F.ExtractAllIndeterminates;
        if (xis.Length != 2)
            throw new();

        var (x, y) = xis.Deconstruct();
        var P0 = F;
        var P1 = F.D(y);
        var tmp = new List<Polynomial<ZnInt, Xi>>();
        foreach (var f in fi)
        {
            var df = f.D(y).Div(I).rem;
            var F0 = P0.Div(f).rem.Div(I).rem;
            var F1 = P1.Div(f).rem.Div(I).rem;
            var _F = F0.ToKPoly(x);
            var _dF = F1.ToKPoly(x);
            var _dFi = FG.NewtonInverse(_dF, int.Max(_dF.Degree, o) + 1);
            var _FDFi = (_dFi * _F).ToPolynomial(F.Indeterminates, x);
            var R1 = (df * _FDFi).Div(f).rem.Div(I).rem;
            var fr = f + R1;
            var c = fr[new(f.Indeterminates, y)];
            tmp.Add(fr * c.Inv());
        }

        return tmp.Order().ToArray();
    }

    /// <summary>
    /// Performs Hensel lifting on a bivariate polynomial over a finite field.
    /// </summary>
    /// <param name="F">The bivariate polynomial to be factorized.</param>
    /// <param name="firr">An array of polynomials representing the initial irreducible factors
    /// of the polynomial F(0, y).</param>
    /// <returns>An array of polynomials representing the factors of the input polynomial.</returns>
    static Polynomial<ZnInt, Xi>[] HenselLifting(Polynomial<ZnInt, Xi> F, Polynomial<ZnInt, Xi>[] firr)
    {
        var (y, x) = F.Indeterminates.Deconstruct();
        var all = firr.ToArray();

        var o = F.DegreeOf(x) + 1;
        var I = F.X(x);
        while (I.Degree < o)
        {
            I = I.Pow(2);
            all = HenselLiftingStep(F, all, I, o);
        }

        return all;
    }

    /// <summary>
    /// Recombines the lifted factors of a bivariate polynomial to form the valid factors
    /// over the original finite field.
    /// </summary>
    /// <param name="F">The original bivariate polynomial over the finite field.</param>
    /// <param name="fi">Array of lifted factors obtained from the Hensel lifting process.</param>
    /// <returns>An array of polynomials representing the recombined factors of the input polynomial.</returns>
    /// <exception cref="ArgumentException">Thrown when no valid recombination of factors is found.</exception>
    static Polynomial<ZnInt, Xi>[] Recombinaison(Polynomial<ZnInt, Xi> F, Polynomial<ZnInt, Xi>[] fi)
    {
        var (y, x) = F.Indeterminates.Deconstruct();
        var o = F.DegreeOf(x) + 1;
        var xo = F.X(x).Pow(o);

        var facts = new List<Polynomial<ZnInt, Xi>>();
        var rem = new HashSet<Polynomial<ZnInt, Xi>>(fi);
        while (rem.Count != 0)
        {
            var sz = rem.Count;
            var combs = rem.AllCombinations();
            foreach (var l in combs.Where(l => l.Length != 0))
            {
                var fact = l.Aggregate(F.One, (acc, li) => acc * li).Div(xo).rem;
                if (F.Div(fact).rem.IsZero())
                {
                    facts.Add(fact);
                    rem.ExceptWith(l);
                    break;
                }
            }

            if (rem.Count == sz)
                throw new();
        }

        return facts.Order().ToArray();
    }

    static Polynomial<Rational, Xi> ZPoly2QPoly(Polynomial<ZnInt, Xi> f, ZnInt c)
    {
        var p = f.P;
        var coefs1 = f.Coefs.ToDictionary(e => e.Key, e => e.Value * c)
            .ToDictionary(
                e => e.Key,
                e => new Rational(e.Value.K * 2 <= p ? e.Value.K : e.Value.K - p))
            .Where(e => !e.Value.IsZero())
            .ToDictionary(e => e.Key, e => e.Value);
        return new Polynomial<Rational, Xi>(f.Indeterminates, Rational.KZero(), new(coefs1));
    }

    /// <summary>
    /// Factorizes a given bivariate polynomial over the rational numbers using Hensel lifting
    /// over finite fields and recombination techniques.
    /// </summary>
    /// <param name="F">The bivariate polynomial to be factorized.</param>
    /// <exception cref="ArgumentException">Thrown when the polynomial contains non-integer coefficients.</exception>
    static Polynomial<Rational, Xi>[] FactorsFxy(Polynomial<Rational, Xi> F)
    {
        if (F.Coefs.Any(e => !e.Value.IsInteger()))
            throw new();

        var (y, x) = F.Indeterminates.Deconstruct();
        var disc = Ring.Discriminant(F.Substitute(F.Zero, x).ToKPoly(y));
        var decomp = IntExt.PrimesDecompositionBigInt(disc.Absolute.Num).Distinct();
        Console.WriteLine($"F({x},{y}) = {F}");
        foreach (var p in IntExt.Primes10000.Except(decomp).Where(p => p < 5000))
        {
            var coefs = F.Coefs.ToDictionary(e => e.Key, e => new ZnInt(p, (int)e.Value.Num))
                .Where(e => !e.Value.IsZero())
                .ToDictionary(e => e.Key, e => e.Value);
            var Fp = new Polynomial<ZnInt, Xi>(F.Indeterminates, ZnInt.ZpZero(p), new(coefs));
            var lc = Fp.LeadingDetails.lc;
            Fp *= lc.Inv();
            try
            {
                var firr = Firr(Fp);
                var lifts = HenselLifting(Fp, firr);
                var factsFp = Recombinaison(Fp, lifts);

                var P0X2 = firr.Aggregate(Fp.One, (acc, e) => e * acc);
                var P0 = factsFp.Aggregate(Fp.One, (acc, e) => e * acc);

                var factsQ = factsFp.Select(f =>
                {
                    var Pi1 = ZPoly2QPoly(f, lc.One);
                    if (F.Div(Pi1).rem.IsZero())
                        return Pi1;

                    var Pic = ZPoly2QPoly(f, lc);
                    if (F.Div(Pic).rem.IsZero())
                        return Pic;

                    return Pi1.One;
                }).ToArray();

                var P1 = factsQ.Aggregate(F.One, (acc, e) => e * acc);

                if (!Fp.Div(P0).rem.IsZero())
                    throw new();
                if (F.Degree != P1.Degree || !F.Div(P1).rem.IsZero())
                    throw new();

                var lc1 = F.LeadingDetails.lc / P1.LeadingDetails.lc;
                if (!(lc1 - 1).IsZero())
                    factsQ = factsQ.Prepend(lc1 * P1.One).ToArray();

                Console.WriteLine($"P(X1,X2) = {Fp}");
                firr.Println($"P(0,X2) = {P0X2}");
                lifts.Println("Hensel Lifting");
                factsFp.Println($"Factors in F{p}[{x},{y}]");
                Console.WriteLine($"{factsFp.Glue(" * ", "({0})")} = {P0}");
                factsQ.Println($"Factors in Q[{x},{y}]");
                Console.WriteLine($"{factsQ.Glue(" * ", "({0})")} = {lc1 * P1}");
                Console.WriteLine();
                return factsQ;
            }
            catch (Exception)
            {
                Console.WriteLine($"########### P = {p} wont work");
            }
        }

        return [F];
    }

    /// <summary>
    /// Generates a random bivariate polynomial with a given number of factors and maximum degree.
    /// </summary>
    /// <param name="nbFactors">The number of factors in the polynomial.</param>
    /// <param name="maxDegree">The maximum degree of the polynomial.</param>
    /// <param name="scalarLT">The leading term of the polynomial is scalar.</param>
    /// <returns>A tuple containing the generated polynomial F and an array of the polynomial factors.</returns>
    static (Polynomial<Rational, Xi> F, Polynomial<Rational, Xi>[] facts)
        GenerateRandomPolynomialFxy(int nbFactors, int maxDegree, bool scalarLT = true)
    {
        var (X2, X1) = Ring.Polynomial(Rational.KZero(), MonomOrder.Lex, "X2", "X1").Deconstruct();
        var x1s = (1 + maxDegree).Range().Select(X1.Pow).ToArray();
        var x2s = maxDegree.Range().Select(X2.Pow).ToArray();
        var x1x2s = x1s.Grid2D(x2s).Select(e => e.t1 * e.t2).Where(f => f.Degree <= maxDegree).Order().ToArray();
        var nb = x1x2s.Length;
        var x2 = X2.ExtractIndeterminate;

        Polynomial<Rational, Xi> Choose()
        {
            var f = (1 + IntExt.Rng.Next(nb)).Range().Select(_ => x1x2s[IntExt.Rng.Next(nb)] * IntExt.Rng.Next(-5, 5))
                .Where(e => !e.IsZero())
                .Aggregate(X1.Zero, (acc, e) => acc + e);
            var degLT = scalarLT ? f.DegreeOf(x2) + 1 : int.Min(maxDegree, f.Degree);
            var i = scalarLT ? 0 : IntExt.Rng.Next(1 + degLT / 2);
            var j = degLT - i;
            var g = f + X2.Pow(j) * X1.Pow(i);
            var mn = g.Coefs.Keys.Aggregate(Monom<Xi>.Gcd);
            var coefs = g.Coefs.ToDictionary(e => e.Key.Div(mn).Item2, e => e.Value);
            return new(g.Indeterminates, g.KZero, new(coefs));
        }

        var facts = nbFactors.Range().Select(_ => Choose()).Order().ToArray();
        var F = facts.Aggregate((a0, a1) => a0 * a1);
        return (F, facts);
    }

    /// <summary>
    /// Filters a bivariate polynomial to determine if it can be factorized.
    /// Checks if the polynomial has at least 2 indeterminates.
    /// Evaluates the polynomial at x = 0 to obtain a univariate polynomial in y.
    /// Checks if the discriminant of the univariate polynomial is non-zero and its degree is greater than 1.
    /// Calculates the derivative of the bivariate polynomial with respect to y and evaluates it at x = 0 to obtain
    /// a univariate polynomial in y.
    /// Calculates the resultant of the univariate polynomials and checks if it is non-zero.
    /// </summary>
    /// <param name="F">The bivariate polynomial to be filtered.</param>
    /// <returns>True if the polynomial satisfies the filtering conditions, false otherwise.</returns>
    static bool FilterRandPolynomialFxy(Polynomial<Rational, Xi> F)
    {
        if (F.NbIndeterminates < 2)
            return false;

        var (y, x) = F.Indeterminates.Deconstruct();
        var F0y = F.Substitute(F.Zero, x).ToKPoly(y);
        if (Ring.Discriminant(F0y).IsZero() || F0y.Degree <= 1)
            return false;

        var DF = F.D(y);
        var DF0y = DF.Substitute(F.Zero, x).ToKPoly(y);
        var res = Ring.FastResultant(F0y, DF0y);
        return !res.IsZero();
    }

    /// <summary>
    /// Performs factorization of a given bivariate polynomial over the rational numbers.
    /// The method generates random bivariate polynomials and attempts to factorize them using Berlekamp's algorithm.
    /// The generated polynomials must satisfy the following conditions:
    /// - Maximum degree of each factor is limited by the provided maxDegreeByFactor parameter.
    /// - Constant term of each factor must not be zero.
    /// - The generated polynomial must have exactly 2 indeterminates.
    /// - The generated polynomial must pass the FilterRandPolynomialFxy filter.
    /// </summary>
    /// <param name="nbPoly">The number of random polynomials to generate and factorize.</param>
    /// <param name="nbFactors">The number of factors in each generated polynomial.</param>
    /// <param name="maxDegreeByFactor">The maximum degree allowed for each factor.</param>
    static void FactorisationRandPolFxy(int nbPoly, int nbFactors, int maxDegreeByFactor)
    {
        var ct = 0;
        while (ct < nbPoly)
        {
            var (F, facts) = GenerateRandomPolynomialFxy(nbFactors, maxDegreeByFactor);
            if (facts.All(f => f.Degree <= maxDegreeByFactor
                               && !f.ConstTerm.IsZero()
                               && f.NbIndeterminates == 2)
                && FilterRandPolynomialFxy(F))
            {
                ct++;
                facts.Println($"F = {F}");
                FactorsFxy(F);
                Console.WriteLine();
            }
        }
    }
    
    static (int i, Polynomial<Rational, Xi> F2) RewritingPolynomial(Polynomial<Rational, Xi> F)
    {
        var (X3, X2, X1, X0) = Ring.Polynomial(Rational.KZero(), MonomOrder.Lex, "X3", "X2", "X1", "X0").Deconstruct();
        var (x3, x2, x1, x0) = X0.Indeterminates.Deconstruct();

        var (y, x) = F.Indeterminates.Deconstruct();
        var Fcoefs = Ring.Decompose(F, y).Item1;
        var F0 = Fcoefs.Select(e => (e.Key.ToKPoly(y).Substitute(X1), e.Value.ToKPoly(x).Substitute(X0)))
            .Aggregate(X0.Zero, (acc, e) => acc + e.Item1 * e.Item2);
        var F0coefsX1 = Ring.Decompose(F0, x1).Item1;

        if (F0coefsX1.Last().Value.Degree == 0)
            return (0, F);

        var seq = 21.Range(-10).OrderBy(i => i * i).ThenDescending().ToArray();
        foreach (var i in seq)
        {
            var subs1 = new List<(Xi, Polynomial<Rational, Xi>)>()
            {
                (x0, X2 + i * X3),
                (x1, X3)
            };

            var F1 = F0.Substitute(subs1);
            var F1coefs = Ring.Decompose(F1, x3).Item1;
            if (F1coefs.Last().Value.Degree != 0)
                continue;

            var F2 = F1coefs.Select(e => (e.Key.ToKPoly(X3).Substitute(F.X(y)), e.Value.ToKPoly(X2).Substitute(F.X(x))))
                .Aggregate(F.Zero, (acc, e) => acc + e.Item1 * e.Item2);

            return (i, F2);
        }

        throw new();
    }

    /// <summary>
    /// The BivariatePolynomialFactorization class provides methods for factorizing bivariate polynomials over a finite field.
    /// </summary>
    public static void Example1_Fp()
    {
        Ring.DisplayPolynomial = MonomDisplay.StarCaret;
        var (X2, X1) = Ring.Polynomial(ZnInt.ZnZero(101), MonomOrder.Lex, "X2", "X1").Deconstruct();

        // (X2^2 + 100*X1^2 + 100) * (X2^2 +  99*X2 + X1^2)
        var P = (X2.Pow(2) + 100 * X1.Pow(2) + 100) * (X2.Pow(2) + 99 * X2 + X1.Pow(2));

        var firr = Firr(P);
        var lifts = HenselLifting(P, firr);
        var facts = Recombinaison(P, lifts);

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
            (X2.Pow(2) - X1.Pow(2) - 1) * (X2.Pow(2) - 2 * X2 + X1.Pow(2)),
            (X2.Pow(2) - 3 * X2 - 4 * X1) * (X2.Pow(2) - 4 * X1 * X2 - 4),
            (X2 - 3 * X1.Pow(2) - 2 * X1) * (X2.Pow(2) - 3 * X1 * X2 - 4),
            (X2.Pow(2) - 3 * X1 - 3) * (X2.Pow(2) - 4 * X2 - 3 * X1.Pow(2)),
            (X2 - 2 * X1 - 3) * (X2.Pow(2) - X2 - 2 * X1),
            (X2 - 2 * X1 - 4) * (X2.Pow(2) - X1.Pow(2) - 3),
            (X2 - 4 * X1.Pow(2) - 2 * X1) * (X2.Pow(2) - 4 * X1.Pow(2) - 3)
        };

        foreach (var F in polys)
            FactorsFxy(F);
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
        // IntExt.RngSeed(25413);
        Ring.DisplayPolynomial = MonomDisplay.StarCaret;
        GlobalStopWatch.Restart();

        foreach (var (n, m) in 4.Range(1).SelectMany(m => (5 - m).Range(2).Select(n => (n, m))))
            FactorisationRandPolFxy(nbPoly: 5, nbFactors: n, maxDegreeByFactor: m);

        Console.Beep();
        GlobalStopWatch.Show();
    }

    public static void Example4_F4787_case()
    {
        Ring.DisplayPolynomial = MonomDisplay.StarCaret;
        var (X2, X1) = Ring.Polynomial(Rational.KZero(), MonomOrder.Lex, "X2", "X1").Deconstruct();
        var F0 = X2.Pow(8) + 4 * X1 * X2.Pow(7) + 2 * X2.Pow(7) + 8 * X1 * X2.Pow(6) + 6 * X2.Pow(6) -
            3 * X1.Pow(3) * X2.Pow(5) - X1.Pow(2) * X2.Pow(5) + 24 * X1 * X2.Pow(5) + 7 * X2.Pow(5) -
            9 * X1.Pow(4) * X2.Pow(4) - 3 * X1.Pow(3) * X2.Pow(4) - 2 * X1.Pow(2) * X2.Pow(4) + 25 * X1 * X2.Pow(4) +
            X2.Pow(4) + 12 * X1.Pow(5) * X2.Pow(3) - 6 * X1.Pow(3) * X2.Pow(3) - 2 * X1.Pow(2) * X2.Pow(3) +
            12 * X1 * X2.Pow(3) - 2 * X2.Pow(3) + 3 * X1.Pow(5) * X2.Pow(2) - 21 * X1.Pow(3) * X2.Pow(2) -
            6 * X1.Pow(2) * X2.Pow(2) - 18 * X2.Pow(2) + 6 * X1.Pow(6) * X2 + 3 * X1.Pow(4) * X2 - 7 * X1.Pow(3) * X2 -
            3 * X1.Pow(2) * X2 + X1 * X2 - 21 * X2 - 9 * X1.Pow(7) - 15 * X1.Pow(4) - 9 * X1.Pow(3) - 4 * X1 - 12;

        FactorsFxy(F0);
    }

    public static void Example5_Substitution()
    {
        var (X0, X1, X2, X3) = Ring.Polynomial(Rational.KZero(), MonomOrder.Lex, (4, "X")).Deconstruct();
        var P = (X0 * X1 + 1) * (X0 + X1 - 3);

        var (x0, x1, x2, _) = X0.Indeterminates.Deconstruct();
        var subs = new List<(Xi, Polynomial<Rational, Xi>)>()
        {
            (x0, X2 + 2 * X3 + 1),
            (x1, X2)
        };

        var Q = P.Substitute(subs).Monic();
        Console.WriteLine(new { P, Q });
        Ring.Decompose(P, x0).Item1.Println("P");
        var Qcoefs = Ring.Decompose(Q, x2).Item1;
        Qcoefs.Println("Q");

        var (Y, X) = Ring.Polynomial(Rational.KZero(), MonomOrder.Lex, "Y", "X").Deconstruct();
        var Q1 = Qcoefs.Select(e => (e.Key.ToKPoly(X2).Substitute(Y), e.Value.ToKPoly(X3).Substitute(X)))
            .Aggregate(X.Zero, (acc, e) => acc + e.Item1 * e.Item2);
        Console.WriteLine(Q1);
        Ring.Decompose(Q1, Y.ExtractIndeterminate).Item1.Println();

        var FactsQ1 = FactorsFxy(Q1);
        var FactsP = new List<Polynomial<Rational, Xi>>();
        foreach (var f in FactsQ1)
        {
            var fcoefs = Ring.Decompose(f, Y.ExtractIndeterminate).Item1;
            var g = fcoefs.Select(e =>
                    (e.Key.ToKPoly(Y).Substitute(X1), e.Value.ToKPoly(X).Substitute((X0 - X1 - 1) / 2)))
                .Aggregate(X0.Zero, (acc, e) => 2 * e.Item1 * e.Item2 + acc);
            FactsP.Add(g.Monic());
        }

        var prod = FactsP.Aggregate(P.One, (acc, e) => e * acc);
        FactsP.Println($"Prod = {prod}, P = Prod:{P.Equals(prod)}");
    }

    public static void Example6_Substitution()
    {
        var (X2, X1) = Ring.Polynomial(Rational.KZero(), MonomOrder.Lex, "X2", "X1").Deconstruct();
        var P = (X2 * X1 + 1) * (X2 + X1 - 3);
        var (i, F) = RewritingPolynomial(P);
        Console.WriteLine($"F = {P}");
        Console.WriteLine($"Substitute {X1} <- {X1 + i * X2}");
        var factsF0 = FactorsFxy(F);
        var factsF = factsF0.Select(fi => fi.Substitute(X1 - i * X2, X1)).Order().ToArray();

        Console.WriteLine();
        Console.WriteLine($"Substitute {X1} <- {X1 - i * X2}");
        factsF.Println($"Factors in Q[{X1},{X2}]");
        Console.WriteLine();
    }

    // F = X1*X2^2 + X1^2*X2 - 3*X1*X2 + X2 + X1 - 3
    // Substitute X1 <- X2 + X1
    // F(X1,X2) = 2*X2^3 + 3*X1*X2^2 - 3*X2^2 + X1^2*X2 - 3*X1*X2 + 2*X2 + X1 - 3
    // ########### P = 3 wont work
    // ########### P = 5 wont work
    // ########### P = 7 wont work
    // ########### P = 11 wont work
    // P(X1,X2) = X2^3 + 10*X1*X2^2 +  7*X2^2 +  9*X1^2*X2 +  7*X1*X2 + X2 +  9*X1 +  7
    // P(0,X2) = X2^3 +  7*X2^2 + X2 +  7
    //     X2 + 13
    //     X2 +  7
    //     X2 +  4
    // Hensel Lifting
    //     X2 +  9*X1 +  7
    //     X2 +  8*X1^2 +  9*X1 +  4
    //     X2 +  9*X1^2 +  9*X1 + 13
    // Factors in F17[X1,X2]
    //     X2 +  9*X1 +  7
    //     X2^2 + X1*X2 + 1
    // (X2 +  9*X1 +  7) * (X2^2 + X1*X2 + 1) = X2^3 + 10*X1*X2^2 +  7*X2^2 +  9*X1^2*X2 +  7*X1*X2 + X2 +  9*X1 +  7
    // Factors in Q[X1,X2]
    //     2*X2 + X1 - 3
    //     X2^2 + X1*X2 + 1
    // (2*X2 + X1 - 3) * (X2^2 + X1*X2 + 1) = 2*X2^3 + 3*X1*X2^2 - 3*X2^2 + X1^2*X2 - 3*X1*X2 + 2*X2 + X1 - 3
    // 
    // 
    // Substitute X1 <- -X2 + X1
    // Factors in Q[X1,X2]
    //     X2 + X1 - 3
    //     X1*X2 + 1
}
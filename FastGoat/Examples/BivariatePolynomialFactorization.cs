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
        var firr = IntFactorisation.Firr(_P0y, a0).Order().ToArray();
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
    static Polynomial<ZnInt, Xi>[] Recombinaison(Polynomial<ZnInt, Xi> F, Polynomial<ZnInt, Xi>[] fi)
    {
        var (y, x) = F.Indeterminates.Deconstruct();
        var o = F.DegreeOf(x) + 1;
        var xo = F.X(x).Pow(o);

        var facts = new List<Polynomial<ZnInt, Xi>>();
        var rem = new HashSet<Polynomial<ZnInt, Xi>>(fi);
        var nbCombs = 1;
        while (rem.Count != 0)
        {
            var sz = rem.Count;
            foreach (var comb in rem.AllCombinationsFromM(nbCombs).ToArray())
            {
                var fact = comb.Aggregate(F.One, (acc, li) => acc * li).Div(xo).rem;
                if (F.Div(fact).rem.IsZero())
                {
                    facts.Add(fact);
                    rem.ExceptWith(comb);
                    nbCombs = comb.Length;
                    break;
                }
            }

            if (rem.Count == sz)
                break;
        }

        return facts.Order().ToArray();
    }
    /// <summary>
    /// Recombines factors of a bivariate polynomial F over rational numbers from a set of polynomials fi over a finite field.
    /// </summary>
    /// <typeparam name="Rational">The type representing rational numbers.</typeparam>
    /// <typeparam name="Xi">The type representing the polynomial variable.</typeparam>
    /// <typeparam name="ZnInt">The type representing integers modulo n.</typeparam>
    /// <param name="F">The polynomial over rationals to be factored.</param>
    /// <param name="fi">An array of polynomials over ZnInt which are intermediate factors.</param>
    /// <param name="c">A constant value over ZnInt.</param>
    /// <returns>
    /// An array of polynomials over rationals that are factors of the given polynomial F.
    /// </returns>
    static Polynomial<Rational, Xi>[] Recombinaison(Polynomial<Rational, Xi> F, Polynomial<ZnInt, Xi>[] fi, ZnInt c)
    {
        var factsQ = new List<Polynomial<Rational, Xi>>();
        var rem = new HashSet<Polynomial<ZnInt, Xi>>(fi);
        var nbCombs = 1;
        
        while (rem.Count != 0)
        {
            var sz = rem.Count;
            foreach (var comb in rem.AllCombinationsFromM(nbCombs).ToArray())
            {
                var fact = comb.Aggregate((acc, e) => e * acc);
                var factQ1 = ZPoly2QPoly(fact, c.One);
                if (F.Div(factQ1).rem.IsZero())
                {
                    factsQ.Add(Primitive(factQ1));
                    rem.ExceptWith(comb);
                    nbCombs = comb.Length;
                    break;
                }

                var factQc = ZPoly2QPoly(fact, c);
                if (F.Div(factQc).rem.IsZero())
                {
                    factsQ.Add(Primitive(factQc));
                    rem.ExceptWith(comb);
                    nbCombs = comb.Length;
                    break;
                }
            }

            if (rem.Count == sz)
                break;
        }

        return factsQ.Order().ToArray();
    }

    /// <summary>
    /// Converts a bivariate polynomial with coefficients in ZnInt to a bivariate polynomial with coefficients in Rational.
    /// </summary>
    /// <param name="f">The bivariate polynomial to be converted.</param>
    /// <param name="c">The characteristic of the field.</param>
    /// <returns>A bivariate polynomial with coefficients in Rational.</returns>
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

    static Polynomial<Rational, Xi> Primitive(Polynomial<Rational, Xi> f)
    {
        var l = f.Coefs.Values.Where(e => !e.IsZero()).Select(e => e.Absolute.Num).Distinct().Order().ToArray();
        return f * new Rational(f.LeadingDetails.lc.Sign, IntExt.GcdBigInt(l));
    }

    static Polynomial<Rational, Xi>[] FactorsFxy(Polynomial<Rational, Xi> F)
    {
        Console.WriteLine("################# Start #################");
        Console.WriteLine();
        Console.WriteLine($"F = {F}");
        Console.WriteLine();
        
        var facts = new List<Polynomial<Rational, Xi>>();
        var Frem = F;
        var (y, x) = F.Indeterminates.Deconstruct();
        while (Frem.Degree > 1)
        {
            var factsQ = FactorsFxyStep(Frem);
            facts.AddRange(factsQ);
            var prod = factsQ.Aggregate(F.One, (acc, f) => f * acc);
            Frem = Frem.Div(prod).quo;
        }
        
        if(!(Frem - 1).IsZero())
            facts.Add(Frem);

        var factsFinal = facts.Order().ToArray();
        
        Console.WriteLine();
        Console.WriteLine($"F = {F}");
        factsFinal.Println($"Final Factors in Q[{x},{y}]");
        Console.WriteLine();
        Console.WriteLine("#################  End  #################");
        Console.WriteLine();
        
        return factsFinal;
    }

    /// <summary>
    /// Factorizes a given bivariate polynomial over the rational numbers using Hensel lifting
    /// over finite fields and recombination techniques.
    /// </summary>
    /// <param name="F">The bivariate polynomial to be factorized.</param>
    /// <exception cref="ArgumentException">Thrown when the polynomial contains non-integer coefficients.</exception>
    static Polynomial<Rational, Xi>[] FactorsFxyStep(Polynomial<Rational, Xi> F)
    {
        if (F.Coefs.Any(e => !e.Value.IsInteger()))
            throw new();

        var (y, x) = F.Indeterminates.Deconstruct();
        var disc = Ring.Discriminant(F.Substitute(F.Zero, x).ToKPoly(y));
        var decomp = IntExt.PrimesDecompositionBigInt(disc.Absolute.Num).Distinct();
        Console.WriteLine($"F({x},{y}) = {F}");
        var pMin = 2 + F.DegreeOf(y) * (2 * F.DegreeOf(x) - 1);
        foreach (var p in IntExt.Primes10000.Except(decomp).Where(p => p >= pMin && p < 5000))
        {
            var coefs = F.Coefs.ToDictionary(e => e.Key, e => new ZnInt(p, (int)e.Value.Num))
                .Where(e => !e.Value.IsZero())
                .ToDictionary(e => e.Key, e => e.Value);
            var Fp = new Polynomial<ZnInt, Xi>(F.Indeterminates, ZnInt.ZpZero(p), new(coefs));
            var lc = Fp.LeadingDetails.lc;
            Fp *= lc.Inv();
            bool isSplitting = false;
            try
            {
                var firr = Firr(Fp);
                isSplitting = F.DegreeOf(y) == firr.Length;
                if (firr.All(f => f.Degree != 1))
                    throw new("Message:No degre 1 factor");
                
                var P0X2 = firr.Aggregate(Fp.One, (acc, e) => e * acc);
                
                var lifts = HenselLifting(Fp, firr.Where(f => f.Degree == 1).ToArray());
                var factsFp = Recombinaison(Fp, lifts);
                if (factsFp.Length == 0)
                    throw new("Message:Fp Recombinaison"); // TODO

                var P0 = factsFp.Aggregate(Fp.One, (acc, e) => e * acc);
                var factsQ = Recombinaison(F, factsFp, lc);
                if (factsQ.Length == 0)
                    throw new("Message:Rational Recombinaison"); // TODO

                factsQ = factsQ.Select(f => Primitive(f)).Order().ToArray();
                var P1 = factsQ.Aggregate(F.One, (acc, e) => e * acc);

                var nbFactsQ = factsQ.Count(f => f.Degree != 0);
                var msg1 = factsFp.Length != F.DegreeOf(y) ? "Yes" : "No";
                var msg2 = factsFp.Length != nbFactsQ ? "Yes" : "No";
                Console.WriteLine($"P(X1,X2) = {Fp}");
                firr.Println($"P(0,X2) = {P0X2}");
                lifts.Println("Hensel Lifting");
                factsFp.Println($"Factors in F{p}[{x},{y}], Z-Recombinaison:{msg1}");
                Console.WriteLine($"{factsFp.Glue(" * ", "({0})")} = {P0}");
                factsQ.Println($"Factors in Q[{x},{y}], Q-Recombinaison:{msg2}");
                Console.WriteLine($"{factsQ.Glue(" * ", "({0})")} = {P1}");
                Console.WriteLine();
                return factsQ;
            }
            catch (Exception e)
            {
                Console.WriteLine($"########### P = {p,-5} wont work. Is split:{isSplitting,-5}. {e.Message}");
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
            var i = scalarLT ? 0 : IntExt.Rng.Next(1 + degLT);
            var j = degLT - i;
            var g = f + X2.Pow(j) * X1.Pow(i);
            var mn = g.Coefs.Keys.Aggregate(Monom<Xi>.Gcd);
            var coefs = g.Coefs.ToDictionary(e => e.Key.Div(mn).Item2, e => e.Value);
            return new(g.Indeterminates, g.KZero, new(coefs));
        }

        var facts = nbFactors.Range().Select(_ => Choose()).Order().ToArray();
        var F = facts.Aggregate((a0, a1) => a0 * a1);
        var degX1 = Ring.Decompose(F, X1.ExtractIndeterminate).Item1.Max(e => e.Key.Degree);
        var degX2 = Ring.Decompose(F, X2.ExtractIndeterminate).Item1.Max(e => e.Key.Degree);
        if (degX1 > degX2)
        {
            var facts0 = new List<Polynomial<Rational, Xi>>();
            var sub = new List<(Xi, Polynomial<Rational, Xi>)>()
            {
                (X1.ExtractIndeterminate, X2),
                (X2.ExtractIndeterminate, X1)
            };
            foreach (var f in facts)
                facts0.Add(f.Substitute(sub));

            facts = facts0.Order().ToArray();
            F = facts.Aggregate((a0, a1) => a0 * a1);
        }

        return (F, facts);
    }

    /// <summary>
    /// Filters a bivariate polynomial to determine if it can be factorized.
    /// Checks if the polynomial has at least 2 indeterminates.
    /// Evaluates the polynomial at x = 0 to obtain a univariate polynomial in y.
    /// Checks if the discriminant of the univariate polynomial is non-zero and its degree is greater than 1.
    /// </summary>
    /// <param name="F">The bivariate polynomial to be filtered.</param>
    /// <returns>True if the polynomial satisfies the filtering conditions, false otherwise.</returns>
    static bool FilterRandPolynomialFxy(Polynomial<Rational, Xi> F)
    {
        if (F.NbIndeterminates < 2)
            return false;

        var (y, x) = F.Indeterminates.Deconstruct();
        var F0y = F.Substitute(F.Zero, x).ToKPoly(y);
        var disc = Ring.Discriminant(F0y);
        return F0y.Degree > 1 && !disc.IsZero();
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
            if (F.LeadingDetails.lm.ContentIndeterminates.Count() == 1 &&
                facts.All(f => f.Degree <= maxDegreeByFactor && f.Degree > 0 && f.NbIndeterminates == 2))
            {
                var (i, F0) = RewritingPolynomial2(F);
                if (F0.Degree == 0)
                    continue;

                ct++;
                var factsF0 = FactorsFxy(F0);
                var (y, x) = F.Indeterminates.Deconstruct();
                var factsF = factsF0.Select(fi => fi.Substitute(F.X(x) - i, x))
                    .Select(f => Primitive(f)).Where(f => !(f - 1).IsZero()).Order().ToArray();

                facts.Println($"F = {F}");
                factsF.Println($"Factors in Q[{x},{y}]");
                var nbFacts = factsF.Count(f => f.Degree != 0);
                if (facts.Length > nbFacts)
                    Console.WriteLine("######## Problem new factors");
                else if (facts.Length < nbFacts)
                    Console.WriteLine("######## Problem missing factors");

                Console.WriteLine();
            }
        }
    }

    /// <summary>
    /// Rewrites the given polynomial <paramref name="F"/> to simplify its leading term.
    /// The method transforms the input polynomial such that its leading term is in the form <c>c * X2^n</c>,
    /// where <c>c</c> is a rational coefficient and <c>n</c> is the degree.
    /// </summary>
    /// <param name="F">The input polynomial of type <see cref="Polynomial{Rational, Xi}"/> with indeterminates
    /// of type <see cref="Xi"/>.</param>
    /// <returns>
    /// A tuple containing:
    /// <list type="bullet">
    /// <item>
    /// <description>An integer <c>i</c> that indicates the index used for substitution.</description>
    /// </item>
    /// <item>
    /// <description>The rewritten polynomial <c>F2</c> of type <see cref="Polynomial{Rational, Xi}"/>.</description>
    /// </item>
    /// </list>
    /// </returns>
    /// <exception cref="Exception">Thrown when no suitable substitution is found to simplify the leading term.</exception>
    /// <remarks>
    /// This method is part of the polynomial factorization process.
    /// Simplifying the leading term of the polynomial is crucial for applying factorization techniques like
    /// Hensel lifting for bivariate polynomials.
    /// </remarks>
    static (int i, Polynomial<Rational, Xi> F2) RewritingPolynomial1(Polynomial<Rational, Xi> F)
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

        var seq = 21.Range(-10).Where(i => i != 0).OrderBy(i => i * i).ThenDescending().ToArray();
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

    static (int i, Polynomial<Rational, Xi> F2) RewritingPolynomial2(Polynomial<Rational, Xi> F)
    {
        var (x, X) = F.IndeterminatesAndVariables.First();
        var seq = 21.Range(-10).OrderBy(i => int.Abs(i)).ThenDescending().ToArray();
        foreach (var i in seq)
        {
            var F1 = F.Substitute(X + i, x);
            if (!FilterRandPolynomialFxy(F1))
                continue;

            return (i, F1);
        }

        return (0, F.Zero); // TODO
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
            FactorisationRandPolFxy(nbPoly: 10, nbFactors: n, maxDegreeByFactor: m);

        Console.Beep();
        GlobalStopWatch.Show(); // Time:12.633s
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
        // Avg Time:318 ms Dev:8.173
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

    /// <summary>
    /// Contains example related to the substitution of variables using method Rewrite in a bivariate polynomial.
    /// </summary>
    public static void Example6_Substitution()
    {
        var (X2, X1) = Ring.Polynomial(Rational.KZero(), MonomOrder.Lex, "X2", "X1").Deconstruct();
        var P = (X2 * X1 + 1) * (X2 + X1 - 3);
        var (i, F) = RewritingPolynomial1(P);
        Console.WriteLine($"F = {P}");
        Console.WriteLine($"Substitute {X1} <- {X1 + i * X2}");
        var factsF0 = FactorsFxy(F);
        var factsF = factsF0.Select(fi => fi.Substitute(X1 - i * X2, X1)).Order().ToArray();

        Console.WriteLine();
        Console.WriteLine($"Substitute {X1} <- {X1 - i * X2}");
        factsF.Println($"Factors in Q[{X1},{X2}]");
        Console.WriteLine();
    }

    /// <summary>
    /// Runs a batch of polynomial substitutions and factorizations.
    /// </summary>
    public static void Example7_SubstitutionBatch()
    {
        // IntExt.RngSeed(2519);
        GlobalStopWatch.Restart();
        for (int m = 2; m <= 4; ++m)
        {
            var ct = 0;
            GlobalStopWatch.AddLap();
            while (ct < 5)
            {
                var (F, facts) = GenerateRandomPolynomialFxy(nbFactors: 2, maxDegree: m, scalarLT: false);
                var (y, x) = F.Indeterminates.Deconstruct();
                if (F.CoefMax(y).Degree == 0 || !FilterRandPolynomialFxy(F))
                    continue;

                var (i, F0) = RewritingPolynomial1(F);
                if (facts.Any(fi => fi.NbIndeterminates != 2) || !FilterRandPolynomialFxy(F0))
                    continue;

                var line = Enumerable.Repeat('#', 30).Glue();
                Console.WriteLine(line);
                facts.Println($"F = {F}");
                Console.WriteLine();
                Console.WriteLine($"Substitute {x} <- {F.X(x) + i * F.X(y)}");
                Console.WriteLine($"F0 = {F0}");
                Console.WriteLine(line);
                Console.WriteLine();

                var factsF0 = FactorsFxy(F0);
                var factsF = factsF0.Select(fi => fi.Substitute(F.X(x) - i * F.X(y), x))
                    .Select(f => Primitive(f)).Where(f => !(f - 1).IsZero()).Order().ToArray();

                var P1 = factsF.Aggregate(F.One, (acc, e) => e * acc);
                var c0 = F.LeadingDetails.lc / P1.LeadingDetails.lc;
                if (!(c0 - 1).IsZero())
                    factsF = factsF.Prepend(c0 * P1.One).ToArray();

                Console.WriteLine();
                facts.Println($"F = {F}");
                factsF.Println($"Factors in Q[{x},{y}]");
                var nbFacts = factsF.Count(f => f.Degree != 0);
                if (facts.Length > nbFacts)
                    Console.WriteLine("######## Problem new factors");
                else if (facts.Length < nbFacts)
                    Console.WriteLine("######## Problem missing factors");

                Console.WriteLine();
                ++ct;
            }

            GlobalStopWatch.Show($"MaxDegree:{m}");
        }

        GlobalStopWatch.Show("End");
        Console.Beep();
    }
    // F = -X1^3*X2^3 - 4*X1^2*X2^3 + 6*X1*X2^3 + 5*X2^3 - 3*X1^3*X2^2 - 22*X1^2*X2^2 + 8*X1*X2^2 + 12*X2^2 +
    // 3*X1^4*X2 + X1^3*X2 - 12*X1^2*X2 + 4*X2 + 12*X1^4 - 3*X1^3 - 6*X1^2
    // Factors in Q[X1,X2]
    //     -1
    //     X1^2*X2 - X1*X2 - X2 + 4*X1^2 - X1 - 2
    //     X1*X2^2 + 5*X2^2 - X1*X2 + 2*X2 - 3*X1^2
    
    public static void Example8_ResultantZero()
    {
        var (X2, X1) = Ring.Polynomial(Rational.KZero(), MonomOrder.Lex, "X2", "X1").Deconstruct();

        var F = X2.Pow(4) + 2 * X1 * X2.Pow(3) + X2.Pow(3) - 2 * X1.Pow(2) * X2.Pow(2) - 3 * X1 * X2.Pow(2) -
            2 * X2.Pow(2) - 12 * X1.Pow(2) * X2 - 4 * X1 * X2 + 10 * X1.Pow(3) + 4 * X1.Pow(2);
        var (i, F0) = RewritingPolynomial2(F);
        Console.WriteLine(new { F });
        Console.WriteLine($"Res(F)(0,0) != 0 : {FilterRandPolynomialFxy(F)}");
        Console.WriteLine(new { i, F0 });
        Console.WriteLine($"Substitute {X1} <- {X1 + i}");
        Console.WriteLine($"Res(F0)(0,0) != 0 : {FilterRandPolynomialFxy(F0)}");
        var facts = FactorsFxy(F0);
        var factsF = facts.Select(f => f.Substitute(X1 - i, X1)).Order().ToArray();
        Console.WriteLine($"F = {F}");
        factsF.Println($"Factors in Q[{X1},{X2}]");
    }

    public static void Example9_ResultantZeroBatch()
    {
        // IntExt.RngSeed(2519);
        var ct = 0;
        while (ct < 50)
        {
            var (F, facts) = GenerateRandomPolynomialFxy(nbFactors: 2, maxDegree: 3, scalarLT: true);
            var (y, x) = F.Indeterminates.Deconstruct();
            if (facts.Any(fi => fi.NbIndeterminates != 2))
                continue;

            var line = Enumerable.Repeat('#', 30).Glue();
            Console.WriteLine(line);
            facts.Println($"F = {F}");
            Console.WriteLine();

            var (i, F0) = RewritingPolynomial2(F);
            if (!FilterRandPolynomialFxy(F0))
                throw new("Res[F](0,0) = 0");
            if (i != 0)
            {
                Console.WriteLine($"Substitute {x} <- {F.X(x) + i}");
                Console.WriteLine($"F0 = {F0}");
            }
            else
                Console.WriteLine("No substitution");

            Console.WriteLine(line);
            Console.WriteLine();

            var factsF0 = FactorsFxy(F0);
            var factsF = factsF0.Select(fi => fi.Substitute(F.X(x) - i, x))
                .Select(f => Primitive(f)).Where(f => !(f - 1).IsZero()).Order().ToArray();

            var P1 = factsF.Aggregate(F.One, (acc, e) => e * acc);
            var c0 = F.LeadingDetails.lc / P1.LeadingDetails.lc;
            if (!(c0 - 1).IsZero())
                factsF = factsF.Prepend(c0 * P1.One).ToArray();

            if (i != 0)
            {
                Console.WriteLine();
                facts.Println($"F = {F}");
                factsF.Println($"Factors in Q[{x},{y}]");
                var nbFacts = factsF.Count(f => f.Degree != 0);
                if (facts.Length > nbFacts)
                    Console.WriteLine("######## Problem new factors");
                else if (facts.Length < nbFacts)
                    Console.WriteLine("######## Problem missing factors");

                Console.WriteLine();
            }

            ++ct;
        }
    }
}
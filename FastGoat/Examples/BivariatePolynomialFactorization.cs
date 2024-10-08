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

    /// <summary>
    /// Compute the primitive polynomial of a bivariate polynomial with coefficients in Rational.
    /// </summary>
    /// <param name="f">The bivariate polynomial.</param>
    /// <returns>The primitive form of <paramref name="f"/>.</returns>
    static Polynomial<Rational, Xi> Primitive(Polynomial<Rational, Xi> f)
    {
        if (f.IsZero())
            return f;

        var arrGcd = f.Coefs.Values.Where(e => !e.IsZero()).Select(e => e.Absolute.Num).Distinct().Order().ToArray();
        var arrLcm = f.Coefs.Values.Select(e => e.Absolute.Denom).Distinct().Order().ToArray();
        return f * new Rational(f.LeadingDetails.lc.Sign * IntExt.LcmBigInt(arrLcm), IntExt.GcdBigInt(arrGcd));
    }

    /// <summary>
    /// Converts a bivariate polynomial represented by a dictionary-based structure into 
    /// a bivariate polynomial using a jagged array structure.
    /// </summary>
    /// <param name="F">The bivariate polynomial to be converted, represented by a dictionary-based structure.</param>
    /// <returns>
    /// A bivariate polynomial represented by a jagged array, where the coefficients are fractional polynomials
    /// with rational coefficients.
    /// </returns>
    static KPoly<FracPoly<Rational>> ConvertToJaggedArrayPoly(Polynomial<Rational, Xi> F)
    {
        var Y = FG.KPoly('Y', FG.QFracPoly('X'));
        var X = Y.KOne.X * Y.One;
        var (y, x) = F.Indeterminates.Deconstruct();
        var Fcoefs = Ring.Decompose(F, y).Item1;
        return Fcoefs.Select(e => e.Key.ToKPoly(y).Substitute(Y) * e.Value.ToKPoly(x).Substitute(X))
            .Aggregate((acc, e) => e + acc);
    }

    /// <summary>
    /// Converts a bivariate polynomial from a jagged array-based structure into a dictionary-based structure.
    /// </summary>
    /// <param name="F">The bivariate polynomial in a jagged array-based structure.</param>
    /// <param name="X">The polynomial used as the first indeterminate.</param>
    /// <param name="Y">The polynomial used as the second indeterminate.</param>
    /// <returns>
    /// A bivariate polynomial represented by a dictionary-based structure, where the coefficients are rational numbers.
    /// </returns>
    static Polynomial<Rational, Xi> ConvertToDictionaryPoly(KPoly<FracPoly<Rational>> F, Polynomial<Rational, Xi> X,
        Polynomial<Rational, Xi> Y)
    {
        var lcm = F.Coefs.Select(f => f.Denom).Aggregate(Ring.Lcm);
        var F0 = new FracPoly<Rational>(lcm) * F;
        if (F0.Coefs.Any(c => !(c.Denom - 1).IsZero()))
            throw new();
        return Primitive(F0.Coefs.Select((c0, i) => c0.Num.Substitute(X) * Y.Pow(i)).Aggregate((acc, e) => e + acc));
    }

    /// <summary>
    /// Factorizes a bivariate polynomial with coefficients in Rational into Padic Integer polynomial factor.
    /// </summary>
    /// <param name="F0y">The bivariate polynomial to be factorized.</param>
    /// <param name="lc">The leading coefficient of the polynomial.</param>
    /// <param name="p">The characteristic of the field.</param>
    /// <param name="o">The degree of the finite field.</param>
    /// <returns>A tuple consisting of a modulus, an array of factors
    static (Modulus po, KPoly<ZnBInt>[] firr, KPoly<ZnBInt>[] firr1) Firr(KPoly<Rational> F0y, Rational lc, int p,
        int o)
    {
        var k = IntExt.SolveAll_k_pow_m_equal_one_mod_n_strict(p, p - 1).First();
        var a0 = ZnBInt.ZnZero(p) + k;
        var po = a0.Details;
        var f0y = IntFactorisation.QPoly2ZnInt(F0y, po);
        var firr = IntFactorisation.Firr(f0y.Monic, a0).ToArray();

        var all = new List<KPoly<ZnBInt>>(firr);
        while (po.O < o && all.Count > 1)
        {
            // ++po;
            po *= 2;
            var tmp = new List<KPoly<ZnBInt>>();
            var fa = IntFactorisation.QPoly2ZnInt(F0y, po) * lc.ToZnBInt(po).Inv();
            foreach (var g in all)
            {
                var gi = IntFactorisation.ZPoly2ZnInt(g, po);
                var y = FG.EPoly(gi);
                var dgi = gi.Derivative.Substitute(y);
                var fi = fa.Substitute(y);
                var dfi = fa.Derivative.Substitute(y);
                var ri = (dgi * fi * dfi.Inv()).Poly;
                tmp.Add(gi + ri);
            }

            all.Clear();
            all = tmp.ToList();
        }

        return (po, firr, all.ToArray());
    }

    /// <summary>
    /// Performs Hensel lifting on a bivariate polynomial over a finite field.
    /// </summary>
    /// <param name="F">The bivariate polynomial to be factorized.</param>
    /// <param name="firr">An array of polynomials representing the initial irreducible factors
    /// of the polynomial F(0, y).</param>
    /// <param name="o">The precision of the lifting.</param>
    /// <returns>An array of polynomials representing the factors of the input polynomial.</returns>
    static Polynomial<ZnBInt, Xi>[] HenselLifting(Polynomial<ZnBInt, Xi> F, Polynomial<ZnBInt, Xi>[] firr, int o)
    {
        var (y, x) = F.Indeterminates.Deconstruct();
        var all = firr.ToArray();

        var I = F.X(x);
        while (I.Degree < o)
        {
            I = I.Pow(2);
            all = HenselLiftingStep(F, all, I, o);
        }

        return all;
    }

    /// <summary>
    /// Performs one step of the Hensel Lifting method to lift the factorization of a bivariate polynomial 
    /// from one modulus to a higher power of the modulus.
    /// </summary>
    /// <param name="F">The original bivariate polynomial to be factorized.</param>
    /// <param name="firr">The array of current factors of the polynomial F.</param>
    /// <returns>An array of polynomials representing the factors of the input polynomial after one lifting step.</returns>
    static Polynomial<ZnBInt, Xi>[] HenselLifting(Polynomial<ZnBInt, Xi> F, Polynomial<ZnBInt, Xi>[] firr)
    {
        var (_, x) = F.Indeterminates.Deconstruct();
        var o = F.DegreeOf(x) + 1;
        return HenselLifting(F, firr, o);
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
    static Polynomial<ZnBInt, Xi>[] HenselLiftingStep(Polynomial<ZnBInt, Xi> F, Polynomial<ZnBInt, Xi>[] fi,
        Polynomial<ZnBInt, Xi> I, int o)
    {
        var xis = F.ExtractAllIndeterminates;
        if (xis.Length != 2)
            throw new();

        var (x, y) = xis.Deconstruct();
        var P0 = F;
        var P1 = F.D(y);
        var tmp = new List<Polynomial<ZnBInt, Xi>>();
        foreach (var f in fi)
        {
            var df = f.D(y).Div(I).rem;
            var F0 = P0.Div(f).rem.Div(I).rem;
            var F1 = P1.Div(f).rem.Div(I).rem;
            var _F = F0.ToKPoly(x);
            var _dF = F1.ToKPoly(x);
            var _dFi = Ring.NewtonInverse(_dF, int.Max(_dF.Degree, o) + 1);
            var _FDFi = (_dFi * _F).ToPolynomial(F.Indeterminates, x);
            var R1 = (df * _FDFi).Div(f).rem.Div(I).rem;
            var fr = f + R1;
            var c = fr[new(f.Indeterminates, y)];
            tmp.Add(fr * c.Inv());
        }

        return tmp.Order().ToArray();
    }

    /// <summary>
    /// Recombines the lifted factors of a bivariate polynomial to form the valid factors
    /// over the original finite field.
    /// </summary>
    /// <param name="F">The original bivariate polynomial over the finite field.</param>
    /// <param name="fi">Array of lifted factors obtained from the Hensel lifting process.</param>
    /// <returns>An array of polynomials representing the recombined factors of the input polynomial.</returns>
    static Polynomial<ZnBInt, Xi>[] Recombinaison(Polynomial<ZnBInt, Xi> F, Polynomial<ZnBInt, Xi>[] fi)
    {
        var (y, x) = F.Indeterminates.Deconstruct();
        var o = F.DegreeOf(x) + 1;
        var xo = F.X(x).Pow(o);

        var facts = new List<Polynomial<ZnBInt, Xi>>();
        var rem = new HashSet<Polynomial<ZnBInt, Xi>>(fi);
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
    /// <param name="F">The polynomial over rationals to be factored.</param>
    /// <param name="fi">An array of polynomials over ZnInt which are intermediate factors.</param>
    /// <param name="c">A constant value over ZnInt.</param>
    /// <returns>
    /// An array of polynomials over rationals that are factors of the given polynomial F.
    /// </returns>
    static Polynomial<Rational, Xi>[] Recombinaison(Polynomial<Rational, Xi> F, Polynomial<ZnBInt, Xi>[] fi, ZnBInt c)
    {
        var factsQ = new List<Polynomial<Rational, Xi>>();
        var rem = new HashSet<Polynomial<ZnBInt, Xi>>(fi);
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
    /// Converts a bivariate polynomial with coefficients in P-adic integers to a bivariate polynomial
    /// with coefficients in Rational.
    /// </summary>
    /// <param name="f">The bivariate polynomial to be converted.</param>
    /// <param name="c">The characteristic of the field.</param>
    /// <returns>A bivariate polynomial with coefficients in Rational.</returns>
    static Polynomial<Rational, Xi> ZPoly2QPoly(Polynomial<ZnBInt, Xi> f, ZnBInt c)
    {
        var coefs1 = f.Coefs.ToDictionary(e => e.Key, e => e.Value * c)
            .ToDictionary(
                e => e.Key,
                e => new Rational(e.Value.K * 2 <= e.Value.Mod ? e.Value.K : e.Value.K - e.Value.Mod))
            .Where(e => !e.Value.IsZero())
            .ToDictionary(e => e.Key, e => e.Value);
        return new Polynomial<Rational, Xi>(f.Indeterminates, Rational.KZero(), new(coefs1));
    }

    /// <summary>
    /// Performs a step in the factorization of a bivariate polynomial with coefficients in Rational.
    /// </summary>
    /// <param name="F">The bivariate polynomial to be factorized.</param>
    /// <param name="name">The bivariate polynomial name to print.</param>
    /// <returns>An array of bivariate polynomials with coefficients in Rational.</returns>
    static Polynomial<Rational, Xi>[] FactorsFxyStep(Polynomial<Rational, Xi> F, string name = "F")
    {
        if (F.Coefs.Any(e => !e.Value.IsInteger()))
            throw new($"F = {F}");

        var (y, x) = F.Indeterminates.Deconstruct();
        Console.WriteLine($"{name}({x},{y}) = {F}");
        var F0y = F.Substitute(F.Zero, x).ToKPoly(y);
        var lc0 = F.LeadingDetails.lc;
        var pMin = 2 + F.DegreeOf(y) * (2 * F.DegreeOf(x) - 1);
        foreach (var (p, o) in IntFactorisation.PSigma(F0y, 600).Where(e => e.p >= pMin))
        {
            try
            {
                var (po, _, firr1) = Firr(F0y, lc0, p, o);
                var firrPAdic = firr1.Select(f => f.Monic.ToPolynomial(F.Indeterminates, y)).ToArray();
                var P0X2 = firrPAdic.Aggregate((acc, e) => e * acc);
                var coefs = F.Coefs.ToDictionary(e => e.Key, e => e.Value.ToZnBInt(po))
                    .Where(e => !e.Value.IsZero())
                    .ToDictionary(e => e.Key, e => e.Value);
                var Fp = new Polynomial<ZnBInt, Xi>(F.Indeterminates, ZnBInt.ZnZero(p, po.O), new(coefs));
                var lc1 = Fp.LeadingDetails.lc;
                Fp *= lc1.Inv();

                var lifts = HenselLifting(Fp, firrPAdic);
                var factsFp = Recombinaison(Fp, lifts);

                var P0 = factsFp.Aggregate(Fp.One, (acc, e) => e * acc);
                var factsQ = Recombinaison(F, factsFp, lc1);
                if (factsQ.Length == 0)
                    throw new("Message:Rational Recombinaison"); // TODO

                factsQ = factsQ.Select(f => Primitive(f)).Order().ToArray();
                var nbFactsQ = factsQ.Count(f => f.Degree != 0);
                var msg1 = factsFp.Length != F.DegreeOf(y) ? "Yes" : "No";
                var msg2 = factsFp.Length != nbFactsQ ? "Yes" : "No";
                Console.WriteLine($"{name}(X1,X2) = {Fp}");
                firrPAdic.Println($"{name}(0,X2) = {P0X2}");
                lifts.Println("Hensel Lifting");
                factsFp.Println($"Factors in F{p}[{x},{y}], Z-Recombinaison:{msg1}");
                Console.WriteLine($"{factsFp.Glue(" * ", "({0})")} = {P0}");
                factsQ.Println($"Factors in Q[{x},{y}], Q-Recombinaison:{msg2}");
                Console.WriteLine();
                return factsQ;
            }
            catch (Exception e)
            {
                Console.WriteLine($"########### P = {p,-5} and O = {o} wont work. {e.Message}");
            }
        }

        return [F];
    }

    /// <summary>
    /// Computes the factors of a bivariate polynomial with coefficients in Rational.
    /// </summary>
    /// <param name="F">The bivariate polynomial to compute factors for.</param>
    /// <returns>An array of bivariate polynomials with coefficients in Rational that are the factors of F.</returns>
    static Polynomial<Rational, Xi>[] FactorsFxy(Polynomial<Rational, Xi> F, bool rewrite = false)
    {
        Console.WriteLine("################# Start #################");
        Console.WriteLine();
        Console.WriteLine($"F = {F}");
        Console.WriteLine();
        var (i0, i1, F1, F2) = (0, 0, F, F);
        var ((x, X), (y, Y)) = F2.IndeterminatesAndVariables.Deconstruct();
        var stack = new Stack<string>();
        stack.Push("F");
        if (rewrite)
        {
            (i0, F1) = RewritingPolynomial1(F);
            if (i0 != 0)
            {
                Console.WriteLine($"Substitute {X}  <-  {X + i0 * Y}");
                Console.WriteLine($"F_1 = {F1}");
                Console.WriteLine();
                stack.Push("F_1");
            }

            (i1, F2) = RewritingPolynomial2(F1);
            if (i1 != 0)
            {
                Console.WriteLine($"Substitute {X}  <-  {X + i1}");
                Console.WriteLine($"F_2 = {F1}");
                Console.WriteLine();
                stack.Push("F_2");
            }
        }

        var facts = new List<Polynomial<Rational, Xi>>();
        var Frem = F2;
        var name = stack.Peek();
        while (Frem.Degree > 0)
        {
            var factsQ = FactorsFxyStep(Frem, name);
            facts.AddRange(factsQ);
            var prod = factsQ.Aggregate(F2.One, (acc, f) => f * acc);
            Frem = Primitive(Frem.Div(prod).quo);
            name = $"sub{name}";
        }

        if (rewrite)
        {
            if (i1 != 0)
            {
                facts.Println($"{stack.Pop()} factors");
                Console.WriteLine($"Substitute {X}  <-  {X - i1}");
                facts = facts.Select(f => Primitive(f.Substitute(X - i1, X))).Order().ToList();
                Console.WriteLine();
            }

            if (i0 != 0)
            {
                facts.Println($"{stack.Pop()} factors");
                Console.WriteLine($"Substitute {X}  <-  {X - i0 * Y}");
                facts = facts.Select(f => Primitive(f.Substitute(X - i0 * Y, X))).Order().ToList();
                Console.WriteLine();
            }
        }

        var lt = F / facts.Aggregate(F.One, (acc, f) => f * acc);
        if (!(lt - 1).IsZero())
            facts.Add(lt);

        var factsFinal = facts.Order().ToArray();

        Console.WriteLine();
        Console.WriteLine($"F = {F}");
        factsFinal.Println($"Final Factors in Q[{x},{y}]");
        Console.WriteLine();
        var check = F.Equals(factsFinal.Aggregate((acc, e) => e * acc));
        Console.WriteLine($"Check:{check}");
        if (!check)
            throw new("Wrong factorisation");

        Console.WriteLine("#################  End  #################");
        Console.WriteLine();

        return factsFinal;
    }

    /// <summary>
    /// Generates a random bivariate polynomial with a given number of factors and maximum degree.
    /// </summary>
    /// <param name="nbFactors">The number of factors in the polynomial.</param>
    /// <param name="maxDegree">The maximum degree of the polynomial.</param>
    /// <param name="scalarLT">The leading term of the polynomial is scalar.</param>
    /// <returns>A tuple containing the generated polynomial F and an array of the polynomial factors.</returns>
    /// <remarks>The result polynomial dont always respects the parameters specifications</remarks>
    static (Polynomial<Rational, Xi> F, Polynomial<Rational, Xi>[] facts)
        GenerateRandomPolynomialFxy(int nbFactors, int maxDegree, bool scalarLT = true, bool multiplicity = false)
    {
        // TODO
        var (X2, X1) = Ring.Polynomial(Rational.KZero(), MonomOrder.Lex, "X2", "X1").Deconstruct();
        var x1s = (1 + maxDegree).Range().Select(X1.Pow).ToArray();
        var x2s = maxDegree.Range().Select(X2.Pow).ToArray();
        var x1x2s = x1s.Grid2D(x2s).Select(e => e.t1 * e.t2).Where(f => f.Degree <= maxDegree).Order().ToArray();
        var nb = x1x2s.Length;

        Polynomial<Rational, Xi> Choose()
        {
            var f = (1 + IntExt.Rng.Next(nb)).Range().Select(_ => x1x2s[IntExt.Rng.Next(nb)] * IntExt.Rng.Next(-5, 5))
                .Where(e => !e.IsZero())
                .Aggregate(X1.Zero, (acc, e) => acc + e);
            var degLT = int.Min(maxDegree, f.Degree + 1);
            var i = scalarLT ? 0 : IntExt.Rng.Next(degLT);
            var j = degLT - i;
            var g = f + X2.Pow(j) * X1.Pow(i);
            if (!g.ConstTerm.IsZero() || g.Coefs.Keys.Aggregate(Monom<Xi>.Gcd).Degree == 0)
                return Primitive(g);

            return Primitive(g + IntExt.RngSign * IntExt.Rng.Next(1, 5));
        }

        var c = multiplicity ? (maxDegree == 2 ? 3 : 2) : 1;
        var facts = (10 * nbFactors).Range().Select(_ => Choose().Pow(IntExt.Rng.Next(1, 1 + c)))
            .DistinctBy(f => f.Monic()).Take(nbFactors).Order().ToArray();
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
            var (y, x) = F.Indeterminates.Deconstruct();
            if (F.LeadingDetails.lm.ContentIndeterminates.Count() == 1 &&
                facts.All(f => f.CoefMax(y).Degree == 0 &&
                               f.Degree <= maxDegreeByFactor && 
                               f.Degree > 0 && 
                               f.NbIndeterminates == 2))
            {
                var (i, F0) = RewritingPolynomial2(F);
                if (F0.Degree == 0)
                    continue;

                ct++;
                var factsF0 = FactorsFxy(F0);
                var factsF = factsF0.Select(fi => fi.Substitute(F.X(x) - i, x))
                    .Select(f => Primitive(f)).Where(f => !(f - 1).IsZero()).Order().ToArray();
                var lt = F / factsF.Aggregate((acc, e) => e * acc);
                if (!(lt - 1).IsZero())
                    factsF = factsF.Prepend(lt).ToArray();

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
        var ((x, X), (y, Y)) = F.IndeterminatesAndVariables.Deconstruct();
        if (F.CoefMax(y).Degree == 0)
            return (0, F);

        var degY = F.DegreeOf(y);
        var seq = (2 * degY + 1).Range(-degY).OrderBy(i => int.Abs(i)).ThenDescending().ToArray();
        foreach (var i in seq)
        {
            var subs1 = new List<(Xi, Polynomial<Rational, Xi>)>()
            {
                (x, X + i * Y),
                (y, Y)
            };

            var F1 = F.Substitute(subs1);
            if (F1.CoefMax(y).Degree != 0)
                continue;

            return (i, Primitive(F1));
        }

        throw new();
    }

    /// <summary>
    /// Rewrites a bivariate polynomial with not null resultant
    /// </summary>
    /// <param name="F">The bivariate polynomial to be rewritten.</param>
    /// <returns>A tuple containing the rewriting constant and the rewritten polynomial.</returns>
    static (int i, Polynomial<Rational, Xi> F2) RewritingPolynomial2(Polynomial<Rational, Xi> F)
    {
        var ((x, X), (y, Y)) = F.IndeterminatesAndVariables.Deconstruct();
        var degY = F.DegreeOf(y);
        var seq = (2 * degY + 1).Range(-degY).OrderBy(i => int.Abs(i)).ThenDescending().ToArray();
        foreach (var i in seq)
        {
            var F1 = F.Substitute(X + i, x);
            if (!FilterRandPolynomialFxy(F1))
                continue;

            return (i, Primitive(F1));
        }

        throw new("Non separable polynomial");
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
        IntExt.RngSeed(25413);
        Ring.DisplayPolynomial = MonomDisplay.StarCaret;
        GlobalStopWatch.AddLap();

        foreach (var (n, m) in 4.Range(1).SelectMany(m => (5 - m).Range(2).Select(n => (n, m))))
            FactorisationRandPolFxy(nbPoly: 10, nbFactors: n, maxDegreeByFactor: m);

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
        GlobalStopWatch.AddLap();
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
    // F = -8*X1^2*X2^5 - 4*X1^3*X2^4 - 6*X1*X2^4 + 20*X1^3*X2^3 - 3*X1^2*X2^3 + 6*X1*X2^3 + 24*X1^2*X2^2 + 3*X2^2 - 5*X1^2*X2 + 6*X1*X2 - X2 - 2*X1
    //     -4*X1*X2^2 - 3*X2 + 1
    //     2*X1*X2^3 + X1^2*X2^2 - 5*X1^2*X2 - X2 - 2*X1
    // Factors in Q[X1,X2]
    //     -1
    //     4*X1*X2^2 + 3*X2 - 1
    //     2*X1*X2^3 + X1^2*X2^2 - 5*X1^2*X2 - X2 - 2*X1
    // 

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
        // IntExt.RngSeed(7219);
        var ct = 0;
        while (ct < 50)
        {
            var (F, facts) = GenerateRandomPolynomialFxy(nbFactors: 2, maxDegree: 3, scalarLT: false);
            var (y, x) = F.Indeterminates.Deconstruct();
            if (facts.Any(fi => fi.NbIndeterminates != 2 || fi.CoefMax(y).Degree > 0))
                continue;

            var line = Enumerable.Repeat('#', 30).Glue();
            Console.WriteLine(line);
            facts.Println($"F = {F}");
            Console.WriteLine();

            var (i, F0) = RewritingPolynomial2(F);
            if (F0.IsZero() || !FilterRandPolynomialFxy(F0))
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

    public static void Example10_Recombinaison_case()
    {
        var (X2, X1) = Ring.Polynomial(Rational.KZero(), MonomOrder.Lex, "X2", "X1").Deconstruct();
        var F = 5 * X1.Pow(2) * X2.Pow(4) + 9 * X1.Pow(3) * X2.Pow(3) + 14 * X1.Pow(2) * X2.Pow(3) +
            4 * X1.Pow(4) * X2.Pow(2) + 6 * X1.Pow(3) * X2.Pow(2) - 3 * X1.Pow(2) * X2.Pow(2) - 18 * X1 * X2.Pow(2) -
            4 * X1.Pow(4) * X2 + X1.Pow(3) * X2 - 14 * X1.Pow(2) * X2 + 10 * X1 * X2 - 2 * X1.Pow(2) - 8;

        var (i0, F0) = RewritingPolynomial2(F);
        Console.WriteLine(new { F });
        Console.WriteLine($"Res(F)(0,0) != 0 : {FilterRandPolynomialFxy(F)}");
        Console.WriteLine(new { F0 });
        Console.WriteLine($"Substitute {X1} <- {X1 + i0}");
        Console.WriteLine($"Res(F0)(0,0) != 0 : {FilterRandPolynomialFxy(F0)}");
        var (i1, F1) = RewritingPolynomial1(F0);
        Console.WriteLine(new { F1 });
        Console.WriteLine($"Substitute {X1} <- {X1 + i1 * X2}");
        var facts = FactorsFxy(F1);
        var factsF1 = facts.Select(f => f.Substitute(X1 - i0 * X2, X1)).Order().ToArray();
        var factsF = factsF1.Select(f => f.Substitute(X1 - i1, X1))
            .Select(f => Primitive(f)).Order().ToArray();
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

    // AECF p405
    public static void Example11_Fp()
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
        var lifts = HenselLifting(P, firr, o);

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

    public static void Example12_FullSubstitutions()
    {
        // IntExt.RngSeed(812665);
        GlobalStopWatch.AddLap();
        foreach (var (n, m) in 4.Range(1).SelectMany(m => (5 - m).Range(2).Select(n => (n, m))))
        {
            var ct = 0;
            GlobalStopWatch.AddLap();
            while (ct < 5)
            {
                var (F, facts) = GenerateRandomPolynomialFxy(nbFactors: n, maxDegree: m, scalarLT: false);
                var (y, x) = F.Indeterminates.Deconstruct();
                if (F.DegreeOf(y) == 0 || facts.Any(f => f.ExtractAllIndeterminates.Length != 2))
                    continue;

                var factsF = FactorsFxy(F, true);
                facts.Println($"F = {F}");
                factsF.Println($"Factors in Q[{x},{y}]");
                var nbfacts = facts.Count(f => f.Degree != 0);
                var nbFacts = factsF.Count(f => f.Degree != 0);
                if (nbfacts > nbFacts)
                    Console.WriteLine("######## Problem new factors");
                else if (nbfacts < nbFacts)
                    Console.WriteLine("######## Problem missing factors");

                Console.WriteLine();
                ++ct;
            }

            GlobalStopWatch.Show($"MaxDegree:{m} NbFactors:{n}");
        }

        GlobalStopWatch.Show("End");
        Console.Beep();
    }

    public static void Example13()
    {
        // F = -5*X1*X2^4 - 2*X2^4 + X1^2*X2^3 + 25*X1^3*X2^2 + 10*X1^2*X2^2 - 5*X1^4*X2 - 20*X1^3*X2 - 3*X1^2*X2 + 12*X1*X2 + 4*X2 + 4*X1^4 - X1^3 - 2*X1^2
        // F = -5*X1*X2.Pow(4)  - 2*X2.Pow(4)  + X1.Pow(2) *X2.Pow(3)  + 25*X1.Pow(3) *X2.Pow(2)  + 10*X1.Pow(2) *X2.Pow(2)  - 5*X1.Pow(4) *X2 - 20*X1.Pow(3) *X2 - 3*X1.Pow(2) *X2 + 12*X1*X2 + 4*X2 + 4*X1.Pow(4)  - X1.Pow(3)  - 2*X1.Pow(2)

        Ring.DisplayPolynomial = MonomDisplay.StarCaret;
        var (X2, X1) = Ring.Polynomial(Rational.KZero(), MonomOrder.Lex, "X2", "X1").Deconstruct();
        var F = -5 * X1 * X2.Pow(4) - 2 * X2.Pow(4) + X1.Pow(2) * X2.Pow(3) + 25 * X1.Pow(3) * X2.Pow(2) +
            10 * X1.Pow(2) * X2.Pow(2) - 5 * X1.Pow(4) * X2 - 20 * X1.Pow(3) * X2 - 3 * X1.Pow(2) * X2 + 12 * X1 * X2 +
            4 * X2 + 4 * X1.Pow(4) - X1.Pow(3) - 2 * X1.Pow(2);

        FactorsFxy(F, true);
    }

    public static void Example14_NonSeparable()
    {
        GlobalStopWatch.AddLap();
        Ring.DisplayPolynomial = MonomDisplay.StarCaret;
        var (X2, X1) = Ring.Polynomial(Rational.KZero(), MonomOrder.Lex, "X2", "X1").Deconstruct();
        var F = X2.Pow(6) - 9 * X1 * X2.Pow(5) + 25 * X1.Pow(2) * X2.Pow(4) - 3 * X1 * X2.Pow(4) -
            44 * X1.Pow(3) * X2.Pow(3) + 27 * X1.Pow(2) * X2.Pow(3) + 72 * X1.Pow(4) * X2.Pow(2) -
            63 * X1.Pow(3) * X2.Pow(2) - 32 * X1.Pow(5) * X2 + 24 * X1.Pow(4) * X2 - 48 * X1.Pow(6) + 36 * X1.Pow(5);

        Console.WriteLine($"F = {F}");
        GlobalStopWatch.AddLap();
        var sff = IntFactorisation.MusserSFF(ConvertToJaggedArrayPoly(F)).Where(e => e.g.Degree > 0)
            .Select(e => (g: ConvertToDictionaryPoly(e.g, X1, X2), e.i)).OrderBy(e => e.i).ToArray();
        sff.Println("SFF");
        GlobalStopWatch.Show("SFF");

        var ((f0, _), (f1, i)) = sff.Deconstruct();
        var facts = FactorsFxy(f0, rewrite: true).Select(f => (g: f, i: 1)).Prepend((g: f1, i)).ToArray();

        Console.WriteLine();
        Console.WriteLine($"F = {F}");
        facts.Println("Factors with multiplicity");

        var prod = facts.Aggregate(F.One, (acc, e) => acc * e.g.Pow(e.i));
        var check = prod.Equals(F);
        Console.WriteLine(new { check });
        GlobalStopWatch.Show();
    }

    public static void Example15_FullFactorisation()
    {
        // IntExt.RngSeed(812665);
        GlobalStopWatch.AddLap();
        foreach (var (n, m) in 4.Range(1).SelectMany(m => (6 - m).Range(2).Select(n => (n, m))))
        {
            var ct = 0;
            GlobalStopWatch.AddLap();
            while (ct < 5)
            {
                var (F, facts) =
                    GenerateRandomPolynomialFxy(nbFactors: n, maxDegree: m, scalarLT: false, multiplicity: true);
                var (y, x) = F.Indeterminates.Deconstruct();
                if (F.DegreeOf(y) == 0
                    || facts.Any(f => f.ExtractAllIndeterminates.Length != 2)
                    || F.Coefs.Max(e => e.Value.Absolute.Num) > 5000)
                    continue;

                var sff = IntFactorisation.MusserSFF(ConvertToJaggedArrayPoly(F)).Where(e => e.g.Degree > 0)
                    .Select(e => (g: ConvertToDictionaryPoly(e.g, F.X(x), F.X(y)), e.i)).OrderBy(e => e.i).ToArray();

                var factorsWithMultiplicity = new List<(Polynomial<Rational, Xi>, int)>();
                var tmpSff = new List<(Polynomial<Rational, Xi>, int)>();
                foreach (var (g, i) in sff)
                {
                    if (g.DegreeOf(y) == 1)
                        factorsWithMultiplicity.Add((g, i));
                    else
                        tmpSff.Add((g, i));
                }

                foreach (var (F0, i) in tmpSff)
                {
                    var factsF0 = FactorsFxy(F0, true);
                    foreach (var g in factsF0)
                    {
                        if (g.Degree != 0)
                            factorsWithMultiplicity.Add((g, i));
                    }
                }

                var factsF = factorsWithMultiplicity.Order().ToArray();
                var prod = factsF.Aggregate(F.One, (acc, e) => e.Item1.Pow(e.Item2) * acc);
                var lc = F.Div(prod).quo;
                if (!lc.Equals(lc.One))
                    factsF = factsF.Prepend((lc, 1)).ToArray();

                facts.Println($"F = {F}");
                sff.Println("SFF");
                factsF.Println($"Factors in Q[{x},{y}]");
                var nbfacts = facts.Count(f => f.Degree != 0);
                var nbFacts = factsF.Count(f => f.Item1.Degree != 0);
                if (nbfacts > nbFacts)
                    Console.WriteLine("######## Problem new factors");
                else if (nbfacts < nbFacts)
                    Console.WriteLine("######## Problem missing factors");
                
                if (nbFacts >= 2 * sff.Length + 1 && factsF.Count(f => f.Item1.CoefMax(y).Degree != 0) > 0)
                    Console.WriteLine("Interesting case");

                prod = factsF.Aggregate(F.One, (acc, e) => e.Item1.Pow(e.Item2) * acc);
                Console.WriteLine($"Check2:{prod.Equals(F)}");

                Console.WriteLine();
                ++ct;
            }

            GlobalStopWatch.Show($"MaxDegree:{m} NbFactors:{n}");
        }

        GlobalStopWatch.Show("End");
        Console.Beep();
    }
    // F = 9*X1^2*X2^9 - 36*X1*X2^9 + 36*X2^9 - 99*X1^3*X2^8 + 393*X1^2*X2^8 - 372*X1*X2^8 - 36*X2^8 + 162*X1^4*X2^7
    // - 690*X1^3*X2^7 + 727*X1^2*X2^7 - 60*X1*X2^7 + 144*X2^7 - 63*X1^5*X2^6 + 345*X1^4*X2^6 - 567*X1^3*X2^6
    // + 431*X1^2*X2^6 - 480*X1*X2^6 + 192*X2^6 - 9*X1^6*X2^5 - 210*X1^5*X2^5 + 647*X1^4*X2^5 - 420*X1^3*X2^5
    // + 516*X1^2*X2^5 - 824*X1*X2^5 - 48*X2^5 + 180*X1^6*X2^4 - 492*X1^5*X2^4 + 166*X1^4*X2^4 - 288*X1^3*X2^4
    // + 608*X1^2*X2^4 + 56*X1*X2^4 + 192*X2^4 + 18*X1^7*X2^3 + 39*X1^6*X2^3 - 228*X1^5*X2^3 + 184*X1^4*X2^3
    // - 64*X1^3*X2^3 + 132*X1^2*X2^3 - 208*X1*X2^3 + 208*X2^3 - 153*X1^7*X2^2 + 159*X1^6*X2^2 + 84*X1^5*X2^2
    // + 380*X1^4*X2^2 - 188*X1^3*X2^2 - 164*X1^2*X2^2 - 384*X1*X2^2 - 16*X2^2 - 9*X1^8*X2 + 48*X1^6*X2 + 24*X1^5*X2
    // - 52*X1^4*X2 - 112*X1^3*X2 + 64*X1*X2 + 64*X2 + 36*X1^8 - 12*X1^6 - 96*X1^5 - 32*X1^4 - 32*X1^3 + 80*X1^2 + 64*X1 + 64
    // Factors in Q[X1,X2]
    //     (X2^3 - 9*X1*X2^2 - X2^2 - X1^2*X2 + 4*X2 + 4*X1^2 + 4, 1)
    //     (3*X1*X2^3 - 6*X2^3 - 3*X1^2*X2^2 + 7*X1*X2^2 + 3*X1^3 - 2*X1 - 4, 2)
    //
    // Sympy check
    // >>> from sympy import *
    // >>> var('X1,X2')
    // >>> F = 9*X1**2 *X2**9  - 36*X1*X2**9  + 36*X2**9  - 99*X1**3 *X2**8  + 393*X1**2 *X2**8  - 372*X1*X2**8  - 36*X2**8
    // + 162*X1**4 *X2**7  - 690*X1**3 *X2**7  + 727*X1**2 *X2**7  - 60*X1*X2**7  + 144*X2**7  - 63*X1**5 *X2**6
    // + 345*X1**4 *X2**6  - 567*X1**3 *X2**6  + 431*X1**2 *X2**6  - 480*X1*X2**6  + 192*X2**6  - 9*X1**6 *X2**5
    // - 210*X1**5 *X2**5  + 647*X1**4 *X2**5  - 420*X1**3 *X2**5  + 516*X1**2 *X2**5  - 824*X1*X2**5  - 48*X2**5
    // + 180*X1**6 *X2**4  - 492*X1**5 *X2**4  + 166*X1**4 *X2**4  - 288*X1**3 *X2**4  + 608*X1**2 *X2**4  + 56*X1*X2**4
    // + 192*X2**4  + 18*X1**7 *X2**3  + 39*X1**6 *X2**3  - 228*X1**5 *X2**3  + 184*X1**4 *X2**3  - 64*X1**3 *X2**3
    // + 132*X1**2 *X2**3  - 208*X1*X2**3  + 208*X2**3  - 153*X1**7 *X2**2  + 159*X1**6 *X2**2  + 84*X1**5 *X2**2
    // + 380*X1**4 *X2**2  - 188*X1**3 *X2**2  - 164*X1**2 *X2**2  - 384*X1*X2**2  - 16*X2**2  - 9*X1**8 *X2 + 48*X1**6 *X2
    // + 24*X1**5 *X2 - 52*X1**4 *X2 - 112*X1**3 *X2 + 64*X1*X2 + 64*X2 + 36*X1**8  - 12*X1**6  - 96*X1**5  - 32*X1**4
    // - 32*X1**3  + 80*X1**2  + 64*X1 + 64
    // >>> factor(F)
    // -(3*X1**3 - 3*X1**2*X2**2 + 3*X1*X2**3 + 7*X1*X2**2 - 2*X1 - 6*X2**3 - 4)**2
    // *(X1**2*X2 - 4*X1**2 + 9*X1*X2**2 - X2**3 + X2**2 - 4*X2 - 4)
    // 
}
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

    static void FactorsFxy(Polynomial<Rational, Xi> F)
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
            try
            {
                var firr = Firr(Fp);
                var lifts = HenselLifting(Fp, firr);
                var factsFp = Recombinaison(Fp, lifts);

                var P0X2 = firr.Aggregate(Fp.One, (acc, e) => e * acc);
                var P0 = factsFp.Aggregate(Fp.One, (acc, e) => e * acc);
                if (!Fp.Equals(P0))
                    throw new();

                var factsQ = factsFp.Select(f =>
                {
                    var coefs1 = f.Coefs.ToDictionary(
                            e => e.Key,
                            e => new Rational(e.Value.K * 2 <= p ? e.Value.K : e.Value.K - p))
                        .Where(e => !e.Value.IsZero())
                        .ToDictionary(e => e.Key, e => e.Value);
                    return new Polynomial<Rational, Xi>(F.Indeterminates, Rational.KZero(), new(coefs1));
                }).ToArray();

                var P1 = factsQ.Aggregate(F.One, (acc, e) => e * acc);
                if (!F.Equals(P1))
                    throw new();

                Console.WriteLine($"P(X1,X2) = {Fp}");
                firr.Println($"P(0,X2) = {P0X2}");
                lifts.Println("Hensel Lifting");
                factsFp.Println($"Factors in F{p}[{x},{y}]");
                Console.WriteLine($"{factsFp.Glue(" * ", "({0})")} = {P0}");
                factsQ.Println($"Factors in Q[{x},{y}]");
                Console.WriteLine($"{factsQ.Glue(" * ", "({0})")} = {P1}");
                Console.WriteLine();
                break;
            }
            catch (Exception)
            {
                Console.WriteLine($"########### P = {p} wont work");
            }
        }
    }
    
    static (Polynomial<Rational, Xi> F, Polynomial<Rational, Xi>[] facts) 
        GenerateRandomPolynomialFxy(int nbFactors, int maxDegree)
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
            var degX2 = f.DegreeOf(x2);
            return f + X2.Pow(degX2 + 1);
        }

        var facts = nbFactors.Range().Select(_ => Choose()).Order().ToArray();
        var F = facts.Aggregate((a0, a1) => a0 * a1);
        return (F, facts);
    }
    
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
}
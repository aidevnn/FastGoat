using System.Numerics;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Padic;

namespace FastGoat.Examples;

public static class PolynomialFactorizationPart2
{
    static double Nu(KPoly<Rational> f)
    {
        var n = f.Degree;
        var norm = f.Coefs.Select(e => Double.Abs(e)).Max();
        return Double.Sqrt(n + 1) * Double.Pow(2, n) * norm;
    }

    static IEnumerable<(int p, int s)> PSigma(KPoly<Rational> f, bool details = false)
    {
        var discPrimes = IntExt.PrimesDecompositionBigInt(BigInteger.Abs(Ring.Discriminant(f).Num))
            .Distinct().ToArray();

        if (details)
            Console.WriteLine($"discPrimes = {discPrimes.Glue("; ")}");
        var primes = IntExt.Primes10000.Except(discPrimes).ToArray();
        var nu = Nu(f);

        var deg = f.Degree + 3;
        var rg = deg.Range(1).ToArray();

        var all = primes.Take(50).Select(p => (p, s: (int)(Double.Log(2 * nu) / Double.Log(p)) + 1))
            .OrderByDescending(e => e.s)
            .ToArray();

        if (details)
            Console.WriteLine($"nu={nu} prime first = {primes.First()} all = {all.Length}");
        
        foreach (var (p, s) in all)
        {
            yield return (p, s);
        }

        // throw new ArgumentException("Prime not found");
    }

    static void PSigma2(KPoly<Rational> f)
    {
        var discPrimes = IntExt.PrimesDecompositionBigInt(Ring.Discriminant(f).Num).Distinct();
        var p = IntExt.Primes10000.Except(discPrimes).First();
    }

    static KPoly<Rational> ZPoly2QPoly(KPoly<ZnBInt> f) =>
        new(f.x, Rational.KZero(), f.Coefs.Select(e => new Rational(e.K * 2 <= e.Mod ? e.K : e.K - e.Mod)).ToArray());

    static KPoly<ZnBInt> QPoly2ZnInt(KPoly<Rational> f, Modulus po)
    {
        var coefs = f.Coefs.Select(e => new ZnBInt(po, e.Num)).ToArray();
        return new(f.x, po.Zero, coefs);
    }

    static KPoly<ZnBInt> ZPoly2ZnInt(KPoly<ZnBInt> f, Modulus po) =>
        new(f.x, po.Zero, f.Coefs.Select(e => new ZnBInt(po, e.K)).ToArray());

    static KPoly<Rational>[] HenselLifting(KPoly<Rational> f, int p, int o, bool details = false)
    {
        if (!f.Coefs.Last().Equals(f.KOne) || f.Coefs.Any(c => !c.Denom.IsOne))
            throw new ArgumentException();

        var nu = Nu(f);

        var ai0 = (new Un(p)).GetGenerators().First()[new(p, 1)];
        var a0 = ZnBInt.KZero(p) + ai0.K;
        var po = a0.Details;
        var f0 = QPoly2ZnInt(f, po);

        var firr0 = PolynomialFactorization.Firr(f0, a0).ToArray();
        var all = new List<KPoly<ZnBInt>>(firr0);
        var o0 = 1;

        while (po.O < o && all.Count > 1)
        {
            ++o0;
            ++po;
            var tmp = new List<KPoly<ZnBInt>>();
            var fa = QPoly2ZnInt(f, po);
            foreach (var g in all)
            {
                var gi = ZPoly2ZnInt(g, po);
                var dgi = new EPoly<ZnBInt>(gi, gi.Derivative);
                var fi = new EPoly<ZnBInt>(gi, fa.Div(gi).rem);
                var dfi = new EPoly<ZnBInt>(gi, fa.Derivative.Div(gi).rem);
                var ri = (dgi * fi / dfi).Poly;
                tmp.Add(gi + ri);
            }

            all.Clear();
            all = tmp.ToList();
        }

        var xp = FG.KPoly(f.x, new ZnBInt(new Modulus(p, o0), 0));
        var F = new KPoly<Rational>(f.x, f.KZero, f.Coefs);

        var listIrr = new List<KPoly<Rational>>();
        var allS = all.ToArray();
        var sz = 0;
        while (allS.Length != sz)
        {
            sz = allS.Length;
            foreach (var combs in allS.AllCombinations().AscendingByCount())
            {
                if (!combs.Any())
                    continue;

                var H = combs.Aggregate(xp.One, (acc, a) => a * acc);
                var G = ZPoly2QPoly(H);
                if (F.Div(G).rem.IsZero())
                {
                    listIrr.Add(G);
                    F /= G;
                    allS = allS.Except(combs).ToArray();
                    break;
                }
            }
        }

        if (allS.Length > 0)
        {
            var Hlast = allS.Aggregate(xp.One, (acc, a) => a * acc);
            var Glast = ZPoly2QPoly(Hlast);
            listIrr.Add(Glast);
        }

        var check = f.Equals(listIrr.Aggregate((a, b) => a * b));
        if (!check)
        {
            throw new Exception();
        }
        
        var listIrr0 = listIrr.OrderBy(q => new KPoly<Rational>(q.x, q.KZero, q.Coefs.Select(Rational.Absolute).ToArray())).ToArray();
        if (details)
        {
            Console.WriteLine($"Prime P = {p}; Sigma = {o}; P^o={BigInteger.Pow(p, o)} 2*Nu={2*nu,0:0.00}");
            Console.WriteLine($"f = {f0} mod {p}");
            Console.WriteLine("Fact(f) = {0} mod {1}", firr0.Order().Glue("*", "({0})"), p);
            Console.WriteLine("Fact(f) = {0} in Z[X]", listIrr0.Glue("*", "({0})"));
            Console.WriteLine();
        }

        return listIrr0;
    }

    public static KPoly<Rational>[] FirrQ(KPoly<Rational> f, bool details = false)
    {
        var m0 = f[f.Degree];
        var f0 = f.Monic;
        var denoms = f.Monic.Coefs.Select(e => e.Denom).Distinct().Where(e => !e.IsOne).ToArray();
        if (denoms.Length == 0)
            denoms = new[] { BigInteger.One, };

        var gcd = IntExt.GcdBigInt(denoms);
        var lcm = denoms.Select(e => e / gcd).Aggregate(BigInteger.One, (acc, a) => acc * a);
        var ct = new Rational(lcm) * new Rational(gcd);
        var f1 = f0 * ct;
        var c = f1[f1.Degree];
        var f2 = c.Pow(f1.Degree - 1) * f1.Substitute(f1.X / c);

        if (!c.Equals(c.One) && details)
        {
            Console.WriteLine($"f0 = {f}");
            Console.WriteLine($"Monic f = {f2}");
        }

        var firr = FirrZ(f2);

        if (!c.Equals(c.One))
        {
            var firr0 = firr.Select(g => g.Substitute(g.X * c).Monic).ToArray();
            var m = m0 * ct.Inv() * c.Pow(f1.Degree - 1).Inv() *
                    firr.Select(g => g.Substitute(g.X * c)[g.Degree]).Aggregate(new Rational(1), (acc, a) => a * acc);
            if(details)
            {
                Console.WriteLine("Finalize");
                if (!m.Equals(m.One))
                    Console.WriteLine("Fact(f0) = {0}*{1} in Z[X]", m, firr0.Order().Glue("*", "({0})"));
                else
                    Console.WriteLine("Fact(f0) = {0} in Z[X]", firr0.Order().Glue("*", "({0})"));

                Console.WriteLine();
            }

            return firr0;
        }

        return firr;
    }

    public static bool EisensteinCriterion(KPoly<Rational> f, bool details = false)
    {
        var a0 = f[0];
        var an = f.Coefs.Last();
        var ai = f.Coefs.SkipLast(1).Where(ai => !ai.IsZero()).ToArray();
        var pi = ai.Length == 0
            ? Array.Empty<int>()
            : ai.Select(a => IntExt.PrimesDecompositionBigInt(BigInteger.Abs(a.Num)).Distinct())
                .Aggregate((a, b) => a.Intersect(b)).ToArray();
        
        if (details)
            Console.WriteLine($"Common Primes {pi.Glue("; ")}");
        
        if (pi.Length == 0)
            return false;
        else
        {
            foreach (var p in pi)
            {
                if (an % p != 0 && a0 % (p * p) != 0)
                {
                    if (details)
                    {
                        Console.WriteLine($"Irreductibility for Eisenstein Criterion for p = {p}");
                        Console.WriteLine($"Fact(f) = {f}");
                        Console.WriteLine();
                    }
                    
                    return true;
                }
            }
        }

        return false;
    }

    public static KPoly<Rational>[] FirrZ(KPoly<Rational> f, bool details = false)
    {
        var discQ = Ring.Discriminant(f).Num;
        var discDecomp = IntExt.PrimesDec(discQ);
        if (details)
        {
            Console.WriteLine($"f = {f}");
            Console.WriteLine($"Disc(f) = {discQ} ~ {discDecomp.AscendingByKey().GlueMap(" * ", "{0}^{1}")}");
        }

        if (EisensteinCriterion(f, details))
        {
            return new[] { f };
        }
        
        foreach (var (p, o) in PSigma(f))
        {
            try
            {
                return HenselLifting(f, p, o, details);
            }
            catch (Exception)
            {
                if (details)
                    Console.WriteLine($"#### Prime {p} and Sigma {o} wont work ####");
            }
        }

        return new[] { f };
    }

    public static void IrreductibleFactorizationZ()
    {
        Console.WriteLine();
        Console.WriteLine("Irreductible Factorization in Z[X]");
        Console.WriteLine();
        Monom.Display = MonomDisplay.StarCaret;
        var X = FG.QPoly('X');

        // AECF example 21.2, page 387
        FirrZ(X.Pow(4) - 1);
        FirrZ(X.Pow(15) - 1);
        FirrZ(X.Pow(11) - 1);
        FirrZ(X.Pow(8) - 40 * X.Pow(6) + 352 * X.Pow(4) - 960 * X.Pow(2) + 576);
        FirrZ(X.Pow(12) - 50 * X.Pow(10) + 753 * X.Pow(8) - 4520 * X.Pow(6) + 10528 * X.Pow(4) - 6720 * X.Pow(2) + 576);

        FirrZ((X + 3) * (X - 5) * (X + 11) * (X - 17));
        FirrZ(X.Pow(6) + 2 * X.Pow(4) - 1);
        FirrZ((X.Pow(3) + 3 * X.Pow(2) + -2) * (X.Pow(3) + -3 * X.Pow(2) + 2));
        FirrZ((X.Pow(3) + 3 * X.Pow(2) + -2) * (X.Pow(3) + -5 * X.Pow(2) + 2 * X + -4));

        FirrZ(X.Pow(6) + 5 * X.Pow(5) + 3 * X.Pow(4) + -7 * X.Pow(3) + -3 * X.Pow(2) + 7 * X + -2);

        // FirrZ(X.Pow(12) + 5 * X.Pow(11) + -202 * X.Pow(10) + -155 * X.Pow(9) + 11626 * X.Pow(8) + -37275 * X.Pow(7)
        //       + -33479 * X.Pow(6) + 547100 * X.Pow(5) + -560012 * X.Pow(4) + -616520 * X.Pow(3) + 351876 * X.Pow(2)
        //       + 146520 * X + -56160); // Bug

        Console.WriteLine("Random Z[X] Polynomials");
        Console.WriteLine();
        for (int j = 0; j < 20; j++)
        {
            var amp = PolynomialFactorization.rnd.Next(2, 20);
            var n = 2 + PolynomialFactorization.rnd.Next(11);
            var degrees = IntExt.Partitions32[n].Where(l => l.All(i => i != 1) && l.Count > 1)
                .OrderBy(i => PolynomialFactorization.rnd.NextDouble())
                .FirstOrDefault(new int[] { 2, 3, 4 }.ToList())
                .ToArray();
            var f0 = degrees.Select(ni => PolynomialFactorization.RandPolySep(Rational.KZero(), amp, ni))
                .Aggregate((a, b) => a * b);
            var f = new KPoly<Rational>('X', f0.KZero, f0.Coefs);
            if (Ring.Discriminant(f).IsZero()) // Separable polynomial
                continue;

            FirrZ(f, details: true);

            /***
             *  Examples of outputs
             * 
                f = X^6 + 12*X^5 + 38*X^4 + 8*X^3 + -9*X^2 + 73*X + -42
                Disc(f) = -95855189967891 = -1^1 * 3^20 * 37^1 * 743^1
                Prime P = 11; Sigma = 5; P^o=161051 Nu=12360.95
                f = X^6 + X^5 + 5*X^4 + -3*X^3 + 2*X^2 + -4*X + 2 mod 11
                Fact(f) = (X + 2)*(X + -4)*(X + -2)*(X^3 + 5*X^2 + -4*X + -4) mod 11
                Fact(f) = (X + 2)*(X^2 + 5*X + -3)*(X^3 + 5*X^2 + -4*X + 7) in Z[X]

                f = X^10 + 9*X^9 + 27*X^8 + 55*X^7 + -29*X^6 + -123*X^5 + -291*X^4 + 370*X^3 + 144*X^2 + -57*X + -126
                Disc(f) = 4406013110764543812064828211149450560 = 2^6 * 3^6 * 5^1 * 19^1 * 23^2 * 47^2 * 53^1 * 271^2
                #### Prime 7 and Sigma 8 wont work ####
                #### Prime 7 and Sigma 9 wont work ####
                Prime P = 13; Sigma = 6; P^o=4826809 Nu=1256602.80
                f = X^10 + -4*X^9 + X^8 + 3*X^7 + -3*X^6 + -6*X^5 + -5*X^4 + 6*X^3 + X^2 + -5*X + 4 mod 13
                Fact(f) = (X + 2)*(X + 3)*(X^2 + 3*X + 6)*(X^2 + -2*X + -4)*(X^4 + 3*X^3 + 2*X^2 + -5*X + -4) mod 13
                Fact(f) = (X^2 + 3*X + 6)*(X^2 + 5*X + -7)*(X^6 + X^5 + 5*X^4 + -8*X^3 + -2*X^2 + 2*X + 3) in Z[X]


             */
        }
    }
}
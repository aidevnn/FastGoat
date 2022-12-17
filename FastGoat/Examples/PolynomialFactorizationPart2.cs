using System.Numerics;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Integers;

namespace FastGoat.Examples;

public static class PolynomialFactorizationPart2
{
    static double Nu(KPoly<Rational> f)
    {
        var n = f.Degree;
        var norm = f.Coefs.Select(e => Double.Abs(e)).Max();
        return Double.Sqrt(n + 1) * Double.Pow(2, n) * norm;
    }

    static IEnumerable<(int p, int s)> PSigma(KPoly<Rational> f)
    {
        var discPrimes = IntExt.PrimesDecompositionBigInt(Ring.Discriminant(f).Num)
            .Concat(IntExt.PrimesDecomposition(f.Degree))
            // .Concat(f.Coefs.SelectMany(e => IntExt.PrimesDecompositionBigInt(e.Num)))
            .Distinct();
        var primes = IntExt.Primes10000.Except(discPrimes).ToArray();
        var nu = Nu(f);

        var rg = f.Degree.Range(2).SkipLast(2).ToArray();
        var all = primes.Take(50).Grid2D(rg)
            .Select(e => (s: e.t2, p: e.t1, pow: float.Pow(e.t1, e.t2)))
            .Where(e => e.pow > 2 * nu).OrderBy(e => e.pow)
            .ToArray();

        foreach (var (s, p, _) in all)
        {
            yield return (p, s);
        }

        throw new ArgumentException("Prime not found");
    }

    static KPoly<Rational> ZPoly2QPoly(KPoly<ZnBInt> f) =>
        new(f.x, Rational.KZero(), f.Coefs.Select(e => new Rational(e.K * 2 <= e.Mod ? e.K : e.K - e.Mod)).ToArray());

    static KPoly<ZnBInt> QPoly2ZnInt(KPoly<Rational> f, ZnBInt.Infos po)
    {
        var coefs = f.Coefs.Select(e => new ZnBInt(po, e.Num)).ToArray();
        return new(f.x, po.Zero, coefs);
    }

    static KPoly<ZnBInt> ZPoly2ZnInt(KPoly<ZnBInt> f, ZnBInt.Infos po) =>
        new(f.x, po.Zero, f.Coefs.Select(e => new ZnBInt(po, e.K)).ToArray());

    static KPoly<Rational>[] HenselLifting(KPoly<Rational> f, int p, int o)
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

        var xp = FG.KPoly(f.x, ZnBInt.KZero(p.Pow(o0)));
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

        Console.WriteLine($"Prime P = {p}; Sigma = {o}; P^o={BigInteger.Pow(p, o)} Nu={nu,0:0.00}");
        Console.WriteLine($"f = {f0} mod {p}");
        Console.WriteLine("Fact(f) = {0} mod {1}", firr0.Order().Glue("*", "({0})"), p);
        Console.WriteLine("Fact(f) = {0} in Z[X]", listIrr.Order().Glue("*", "({0})"));
        Console.WriteLine();

        return listIrr.ToArray();
    }


    static KPoly<Rational>[] FirrZ(KPoly<Rational> f)
    {
        Console.WriteLine($"f = {f}");
        var discQ = Ring.Discriminant(f).Num;
        var discDecomp = IntExt.PrimesDec(discQ);
        Console.WriteLine($"Disc(f) = {discQ} = {discDecomp.AscendingByKey().GlueMap(" * ", "{0}^{1}")}");
        foreach (var (p, o) in PSigma(f))
        {
            try
            {
                return HenselLifting(f, p, o);
            }
            catch (Exception e)
            {
                Console.WriteLine($"#### Prime {p} and Sigma {o} wont work ####");
            }
        }

        return Array.Empty<KPoly<Rational>>();
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

            FirrZ(f);

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
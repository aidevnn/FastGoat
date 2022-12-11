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
        var discPrimes = IntExt.PrimesDecompositionBigInt(Ring.Discriminant(f).Num).Distinct();
        var primes = IntExt.Primes10000.Except(discPrimes);
        var nu = Nu(f);

        foreach (var p in primes.Take(20))
        {
            for (int s = 1; s < f.Degree; ++s)
            {
                if (double.Pow(p, s) > 2 * nu)
                    yield return (p, s);
            }
        }

        throw new ArgumentException("Prime not found");
    }

    static KPoly<Rational> Padic2QPoly(KPoly<Padic> f, bool signed = true) =>
        new(f.x, Rational.KZero(), f.Coefs.Select(e => new Rational(signed ? e.ToSignedBigInt : e.ToBigInt)).ToArray());

    static KPoly<Padic> QPoly2Padic(KPoly<Rational> f, int p, int o)
    {
        var coefs = f.Coefs.Select(e => Padic.Convert(p, o, e.Num)).ToArray();
        return new(f.x, new Padic(p, o), coefs);
    }

    static KPoly<Padic> Resize(KPoly<Padic> f, int o0) =>
        new(f.x, f.KZero.Resize(o0), f.Coefs.Select(e => e.Resize(o0)).ToArray());

    static KPoly<Rational>[] HenselLifting(KPoly<Rational> f, int p, int o)
    {
        if (!f.Coefs.Last().Equals(f.KOne) || f.Coefs.Any(c => !c.Denom.IsOne))
            throw new ArgumentException();

        var nu = Nu(f);

        var padicZero = new Padic(p, 1);
        var a0 = (new Un(p)).GetGenerators().First()[new(p, 1)];
        var f0 = QPoly2Padic(f, p, 1);

        var firr0 = PolynomialFactorization.Firr(f0, a0 + padicZero).ToArray();
        var all = new List<KPoly<Padic>>(firr0);
        var o0 = 1;

        while (o0 < o && all.Count > 1)
        {
            var tmp = new List<KPoly<Padic>>();
            o0 += 1;
            var fa = QPoly2Padic(f, p, o0);
            foreach (var g in all)
            {
                var gi = Resize(g, o0);
                var dgi = new EPoly<Padic>(gi, gi.Derivative);
                var fi = new EPoly<Padic>(gi, fa.Div(gi).rem);
                var dfi = new EPoly<Padic>(gi, fa.Derivative.Div(gi).rem);
                var ri = (dgi * fi / dfi).Poly;
                tmp.Add(gi + ri);
            }

            all.Clear();
            all = tmp.ToList();
        }


        var xp = FG.PadicPoly(p, o0);
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
                var G = Padic2QPoly(H);
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
            var Glast = Padic2QPoly(Hlast);
            listIrr.Add(Glast);
        }

        Console.WriteLine($"Prime P = {p}; Sigma = {o}; P^o={BigInteger.Pow(p, o)} Nu={nu,0:0.00}");
        Console.WriteLine($"f = {f0} mod {p}");
        Console.WriteLine("Fact(f) = {0} mod {1}", firr0.Order().Glue("*", "({0})"), p);
        Console.WriteLine("Fact(f) = {0} in Z[X]", listIrr.Order().Glue("*", "({0})"));
        var check = f.Equals(listIrr.Aggregate((a, b) => a * b));
        if (!check)
            Console.WriteLine("Validity Product {0}", check);

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

    public static void IrrductibleFactorizationZ()
    {
        Console.WriteLine();
        Console.WriteLine("Irreductible Factorization in Z[X]");
        Console.WriteLine();
        Monom.Display = MonomDisplay.StarCaret;
        var x = FG.QPoly('X');

        // AECF example 21.2, page 387
        FirrZ(x.Pow(4) - 1);
        FirrZ(x.Pow(8) - 40 * x.Pow(6) + 352 * x.Pow(4) - 960 * x.Pow(2) + 576);
        FirrZ(x.Pow(12) - 50 * x.Pow(10) + 753 * x.Pow(8) - 4520 * x.Pow(6) + 10528 * x.Pow(4) - 6720 * x.Pow(2) + 576);

        FirrZ((x + 3) * (x - 5) * (x + 11) * (x - 17));
        FirrZ(x.Pow(6) + 2 * x.Pow(4) - 1);
        FirrZ((x.Pow(3) + 3 * x.Pow(2) + -2) * (x.Pow(3) + -3 * x.Pow(2) + 2));
        FirrZ((x.Pow(3) + 3 * x.Pow(2) + -2) * (x.Pow(3) + -5 * x.Pow(2) + 2 * x + -4));

        FirrZ(x.Pow(6) + 5 * x.Pow(5) + 3 * x.Pow(4) + -7 * x.Pow(3) + -3 * x.Pow(2) + 7 * x + -2);

        Console.WriteLine("Random Z[X] Polynomials");
        Console.WriteLine();
        for (int j = 0; j < 20; j++)
        {
            var p = IntExt.Primes10000[PolynomialFactorization.rnd.Next(5)]; // 2, 3, 5, 7, 11
            var n = 2 + PolynomialFactorization.rnd.Next(9);
            var degrees = IntExt.Partitions32[n].Where(l => l.All(i => i != 1) && l.Count > 1)
                .OrderBy(i => PolynomialFactorization.rnd.NextDouble())
                .FirstOrDefault(new int[] { 2, 3, 4 }.ToList())
                .ToArray();
            var f = degrees.Select(ni => PolynomialFactorization.RandPolySep(Rational.KZero(), p, ni))
                .Aggregate((a, b) => a * b);
            if (Ring.Discriminant(f).IsZero())
                continue;

            FirrZ(f);
        }
    }
}
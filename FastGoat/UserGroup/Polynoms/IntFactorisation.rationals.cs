using System.Numerics;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Padic;

namespace FastGoat.UserGroup.Polynoms;

public static partial class IntFactorisation
{
    /// <summary>
    /// Computes the square norm of a vector row or column
    /// </summary>
    /// <param name="v">KMatrix<Rational></param>
    /// <returns>Rational</returns>
    /// <exception cref="ArgumentException"></exception>
    public static Rational SquareNorm2(KMatrix<Rational> v)
    {
        if (v.M == 1)
            return (v * v.T)[0, 0];
        else if (v.N == 1)
            return (v.T * v)[0, 0];

        throw new ArgumentException();
    }

    /// <summary>
    /// Swap two rows of a matrix
    /// </summary>
    /// <param name="i"></param>
    /// <param name="j"></param>
    /// <param name="A"></param>
    /// <typeparam name="T"></typeparam>
    static void SwapRows<T>(int i, int j, T[,] A)
    {
        var cols = A.GetLength(1);
        for (int k = 0; k < cols; k++)
            (A[i, k], A[j, k]) = (A[j, k], A[i, k]);
    }

    static Rational Round(Rational e)
    {
        var (num, denom) = e;
        var (q, r) = BigInteger.DivRem(num, denom);
        var rs = r.Sign;
        var r0 = r * rs * 2;
        if (r0 < denom || (r0 == denom && BigInteger.IsEvenInteger(q)))
            return new(q, 1);

        return new(q + rs, 1);
    }

    public static KMatrix<Rational> LLL(KMatrix<Rational> v)
    {
        var n = v.N;
        var w = v.Cols;
        var (Ws, M) = Ring.GramSchmidt(v);
        var ws = Ws.Cols;
        var N = M.Coefs;
        int i = 1;
        while (i < n)
        {
            for (int j = i - 1; j >= 0; j--)
            {
                var ruij = Round(N[i, j]);
                w[i] -= ruij * w[j];
                for (int k = 0; k <= j; k++)
                {
                    N[i, k] -= ruij * N[j, k];
                }
            }

            if (i >= 1)
            {
                var wsip2 = SquareNorm2(ws[i - 1]);
                var wsi2 = SquareNorm2(ws[i]);
                if (wsip2.CompareTo(2 * wsi2) > 0)
                {
                    var a = N[i, i - 1];
                    var b = a * wsip2 / (wsi2 + a.Pow(2) * wsip2);
                    (ws[i - 1], ws[i]) = (ws[i] + a * ws[i - 1], ws[i - 1] - b * (ws[i] + a * ws[i - 1]));
                    (w[i - 1], w[i]) = (w[i], w[i - 1]);
                    SwapRows(i - 1, i, N);
                    for (int k = i - 1; k < n; k++)
                    {
                        (N[k, i - 1], N[k, i]) = (b * N[k, i - 1] + (1 - a * b) * N[k, i], N[k, i - 1] - a * N[k, i]);
                    }

                    i--;
                }
                else
                {
                    i++;
                }
            }
            else
            {
                i++;
            }
        }

        return KMatrix<Rational>.MergeSameRows(w);
    }

    static double Nu(KPoly<Rational> f)
    {
        var n = f.Degree;
        var norm = f.Coefs.Select(e => Double.Abs(e)).Max();
        return Double.Sqrt(n + 1) * Double.Pow(2, n) * norm;
    }

    static IEnumerable<(int p, int s)> PSigma(KPoly<Rational> f, bool details = false)
    {
        var disc = Ring.Discriminant(f).Num;
        var nu = Nu(f);
        var all = IntExt.Primes10000.Where(p => !BigInteger.Remainder(disc, p).IsZero).Take(150)
            .Select(p => (p, s: (int)(Double.Log(2 * nu) / Double.Log(p)) + 1))
            .OrderByDescending(e => e.s).ToArray();

        foreach (var (p, s) in all)
        {
            yield return (p, s);
        }

        // throw new ArgumentException("Prime not found");
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

    public static (KPoly<ZnBInt> f0, KPoly<ZnBInt>[] firr0, KPoly<ZnBInt>[] firr) HenselLifting(KPoly<Rational> f, int p, int o)
    {
        if (!f.Coefs.Last().Equals(f.KOne) || f.Coefs.Any(c => !c.Denom.IsOne))
            throw new ArgumentException();

        var ai0 = (new Un(p)).GetGenerators().First()[new(p, 1)];
        var a0 = ZnBInt.KZero(p) + ai0.K;
        var po = a0.Details;
        var f0 = QPoly2ZnInt(f, po);

        // GlobalStopWatch.AddLap();
        // var firr0 = Firr(f0, a0).ToArray();
        var firr0 = BerlekampProbabilisticVShoup(f0, a0).ToArray();
        // var firr0 = BerlekampProbabilisticAECF(f0, a0).ToArray();
        // var firr0 = CantorZassenhausVShoup(f0, a0).ToArray();
        // GlobalStopWatch.Show("Firr");
        
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
                var ri = ((dgi * fi) / dfi).Poly;
                tmp.Add(gi + ri);
            }

            all.Clear();
            all = tmp.ToList();
        }

        return (f0, firr0, all.ToArray());
    }

    public static KPoly<Rational>[] HenselLiftingNaive(KPoly<Rational> f, int p, int o, bool details = false)
    {
        if (!f.Coefs.Last().Equals(f.KOne) || f.Coefs.Any(c => c.Denom!=1))
            throw new ArgumentException();

        var (f0, firr0, allS) = HenselLifting(f, p, o);
        var o0 = allS.Max(pl => pl.KZero.Details.O);
        var xp = FG.KPoly(f.x, new ZnBInt(new Modulus(p, o0), 0));
        var F = new KPoly<Rational>(f.x, f.KZero, f.Coefs);

        var listIrr = new List<KPoly<Rational>>();
        var sz = 0;
        while (allS.Length != sz)
        {
            sz = allS.Length;
            foreach (var combs in allS.AllCombinations())
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
            var nu = Nu(f);
            Console.WriteLine($"Prime P = {p}; Sigma = {o}; P^o={BigInteger.Pow(p, o)} 2*Nu={2 * nu,0:0.00}");
            Console.WriteLine($"f = {f0} mod {p}");
            Console.WriteLine("Fact(f) = {0} mod {1}", firr0.Order().Glue("*", "({0})"), p);
            Console.WriteLine("Fact(f) = {0} in Z[X]", listIrr0.Glue("*", "({0})"));
            Console.WriteLine();
        }

        return listIrr0;
    }


    static KMatrix<Rational> Lattice(KPoly<Rational> F, KPoly<ZnBInt>[] irrs, int n, int p, int sigma, int tau)
    {
        var s = irrs.Length;
        var m = n + s;
        var mat = new KMatrix<Rational>(Rational.KZero(), m, m);

        var Fp = QPoly2ZnInt(F, new(p, sigma));
        for (int i = 0; i < s; i++)
            mat.Coefs[i, i] = Rational.KOne();

        for (int i = 0; i < n; i++)
            mat.Coefs[i + s, i + s] = new Rational(p).Pow(sigma - tau);

        for (int i = 0; i < s; i++)
        {
            var fi = irrs[i];
            var gi = (Fp * fi.Derivative) / fi;
            // Console.WriteLine((Fp * fi.Derivative).Div(fi).rem);
            for (int j = 0; j < n; j++)
            {
                var z = ZnBInt.Truncate(gi[j], tau);
                mat.Coefs[i, j + s] = new Rational(z.trunc.ToSignedBigInt);
            }
        }

        return LLL(mat.T);
    }

    public static KPoly<Rational>[] VanHoeijFactorization(KPoly<Rational> f, int max = 2)
    {
        var discQ = Ring.Discriminant(f).Num;
        Console.WriteLine($"f = {f}");
        var discDecomp = IntExt.PrimesDec(discQ);
        Console.WriteLine($"Disc(f) = {discQ} ~ {discDecomp.AscendingByKey().GlueMap(" * ", "{0}^{1}")}");

        foreach (var p in IntExt.Primes10000.Where(p => !BigInteger.Remainder(discQ, p).IsZero))
        {
            try
            {
                return SearchVanHoeij(f, p, max);
            }
            catch (Exception e)
            {
                Console.WriteLine($"#### Prime {p} wont work ####");
            }
        }

        throw new();
    }

    static KPoly<Rational>[] SearchVanHoeij(KPoly<Rational> P, int p, int max = 2)
    {
        var n = P.Degree;
        var x = FG.ZPoly(p);
        var P0 = P.Coefs.Select((e, i) => (int)BigInteger.Remainder(e.Num, p) * x.Pow(i)).Aggregate(x.Zero, (sum, xi) => sum + xi);
        var p0 = Un.FirstGen(p);
        var irrs = FirrFsep(P0, p0);
        var s = irrs.Count;
        var m = s + n;
        var nu = Double.Sqrt(s + n * Double.Pow(s / 2.0, 2));

        var norm2 = Double.Sqrt(P.Coefs.Sum(f => Double.Pow(f.Abs(), 2)));
        var boundTau = n * Double.Pow(2, n) * norm2;
        var tau = (int)(Double.Log(boundTau) / Double.Log(p)) + 1;
        var pPowTau = Double.Pow(p, tau);
        var sqrtN = Double.Sqrt(n);
        var theta = sqrtN * ((1 + sqrtN / 2.0) * Double.Pow(2, (m - 1) / 2.0) * (1 + nu) * nu + 0.5) * pPowTau;
        var boundSigma = 2 * Double.Pow(theta, n) * Double.Pow(norm2, n - 1);
        var sigma = (int)((Double.Log(2) + n * Double.Log(theta) + (n - 1) * Double.Log(norm2)) / Double.Log(p)) + 1;
        var nb = Int32.Min(Int32.Max(max, Int32.Abs(sigma / tau)), max);
        Console.WriteLine(
            $"Search Van-Hoeij Prime P = {p} Tau = {tau} nb = {nb} Theta = {theta} BoundSigma = {boundSigma} MaxSigma = {sigma}; Nu = {nu}");
        for (int i = 2; i <= nb; i++)
        {
            var sigma2 = tau * i;
            try
            {
                var (f0, firr0, irrs2) = HenselLifting(P, p, sigma2);
                var lll = Lattice(P, irrs2, n, p, sigma2, tau);
                var cols = lll.Cols;
                var rgM = m.Range();
                var rgS = s.Range();
                var t = rgM.Where(i0 => rgM.Where(j => j > i0).All(j => Double.Sqrt(SquareNorm2(cols[j])) > nu)).Min();

                var bs = (t + 1).Range().Select(i0 => cols[i0].T.Extract(0, 1, 0, s)).ToArray();
                var coefs = bs.Select(mat => mat.Select(r => (int)r.Num).ToArray()).ToArray();
                var one0 = irrs2[0].One;
                var polys = coefs.Select(l => rgS.Aggregate(one0, (prod, i0) => prod * irrs2[i0].Pow(l[i0])))
                    .Select(ZPoly2QPoly).ToArray();

                var one1 = polys[0].One;
                Console.WriteLine($"       Found Sigma = {sigma2} = {i}*Tau");
                Console.WriteLine($"f = {f0} mod {p}");
                Console.WriteLine("Fact(f) = {0} mod {1}", firr0.Glue("*", "({0})"), p);
                bs.Println("LLL Combinaisons");
                Console.WriteLine("Fact(f) = {0} in Z[X]", polys.Glue("*", "({0})"));
                Console.WriteLine($"f = Prod[Fact(f)] : {P.Equals(polys.Aggregate(one1, (prod, e) => prod * e))}");
                Console.WriteLine();
                return polys;
            }
            catch (Exception e)
            {
                Console.WriteLine($"#### Prime {p} Sigma {sigma2} = {i}*Tau wont work ####");
            }
        }

        throw new();
    }

    public static (KPoly<K> nf, KPoly<K> nx) Deflate<K>(KPoly<K> f) 
        where K : struct, IFieldElt<K>, IRingElt<K>, IElt<K>
    {
        var pows = f.Coefs.Select((c, i) => (c, i)).Where(e => e.i != 0 && !e.c.IsZero()).ToArray();
        var gcd = IntExt.Gcd(pows.Select(e => e.i).ToArray());
        var coefs = f.Coefs.Where((c, i) => i % gcd == 0).ToArray();
        return (new(f.x, f.KZero, coefs), f.X.Pow(gcd));
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

        var firr = FirrZ(f2, details);

        if (!c.Equals(c.One))
        {
            var firr0 = firr.Select(g => g.Substitute(g.X * c).Monic).ToArray();
            var m = m0 * ct.Inv() * c.Pow(f1.Degree - 1).Inv() *
                    firr.Select(g => g.Substitute(g.X * c)[g.Degree]).Aggregate(new Rational(1), (acc, a) => a * acc);
            if (details)
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

        if (details && pi.Length > 0)
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

    public static KPoly<Rational>[] FirrZ2(KPoly<Rational> f, bool details = false)
    {
        if (f.Degree == 1)
            return new[] { f };
        
        var (f0, x0) = Deflate(f);
        if (x0.Degree == 1)
        {
            return FirrZ(f0, details);
        }
        else if(f0.Degree == 1)
        {
            return FirrZ(f, details);
        }
        else
        {
            var facts = FirrZ(f0, details);
            if (facts.Length == 1)
                return FirrZ(f, details);

            var facts2 = facts.SelectMany(f1 => FirrZ2(f1.Substitute(x0), details)).ToArray();
            if (details)
            {
                Console.WriteLine($"f = {f}");
                Console.WriteLine("Fact(f) = {0} in Z[X]", facts2.Glue("*", "({0})"));
                Console.WriteLine();
            }

            return facts2;
        }
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

        if (f.Degree == 1)
            return new[] { f };

        // if (EisensteinCriterion(f, details))
        // {
        //     return new[] { f };
        // }

        foreach (var (p, o) in PSigma(f))
        {
            try
            {
                return HenselLiftingNaive(f, p, o, details);
            }
            catch (Exception)
            {
                if (details)
                    Console.WriteLine($"#### Prime {p} and Sigma {o} wont work ####");
            }
        }

        return new[] { f };
    }
}
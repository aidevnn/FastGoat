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
    /// <param name="v">Matrix of rationals</param>
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
    public static void SwapRows<T>(int i, int j, T[,] A)
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

    public static KMatrix<K> LLL<K>(KMatrix<K> A) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>, IFloatElt<K>
    {
        var n = A.N;
        var w = A.Cols;
        var (Ws, M) = Ring.GramSchmidt2(A);
        var ws = Ws.Cols;
        var N = new KMatrix<K>(M.Coefs);
        int i = 1;
        while (i < n)
        {
            var wi = w[i];
            for (int j = i - 1; j >= 0; j--)
            {
                var ruij = N[i, j].RoundEven;
                // w[i] -= ruij * w[j];
                var wj = w[j];
                for (int k = 0; k < n; k++)
                    wi.Coefs[k, 0] -= ruij * wj.Coefs[k, 0];

                for (int k = 0; k <= j; k++)
                    N.Coefs[i, k] -= ruij * N[j, k];
            }

            if (i >= 1)
            {
                var wsip2 = A.KZero;
                var wsi2 = A.KZero;
                var wsip = ws[i - 1];
                var wsi = ws[i];
                var wip = w[i - 1];
                for (int k = 0; k < n; k++)
                {
                    wsip2 += wsip.Coefs[k, 0].Pow(2);
                    wsi2 += wsi.Coefs[k, 0].Pow(2);
                }

                if (wsip2.CompareTo(2 * wsi2) > 0)
                {
                    var a = N[i, i - 1];
                    var b = a * wsip2 / (wsi2 + a.Pow(2) * wsip2);
                    for (int k = 0; k < n; k++)
                    {
                        var tmp = wsi.Coefs[k, 0] + a * wsip.Coefs[k, 0];
                        (wsip.Coefs[k, 0], wsi.Coefs[k, 0]) = (tmp, wsip.Coefs[k, 0] - b * tmp);
                        (wip.Coefs[k, 0], wi.Coefs[k, 0]) = (wi.Coefs[k, 0], wip.Coefs[k, 0]);
                    }

                    for (int k = 0; k < n; k++)
                        (N.Coefs[i - 1, k], N.Coefs[i, k]) = (N[i, k], N[i - 1, k]);

                    var coef = 1 - a * b;
                    for (int k = i - 1; k < n; k++)
                        (N.Coefs[k, i - 1], N.Coefs[k, i]) = (b * N[k, i - 1] + coef * N[k, i], N[k, i - 1] - a * N[k, i]);

                    i--;
                }
                else
                    i++;
            }
            else
                i++;
        }

        return KMatrix<K>.MergeSameRows(w);
    }

    public static KMatrix<Rational> LLLtheoric(KMatrix<Rational> v)
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

    static IEnumerable<(int p, int s)> PSigma(KPoly<Rational> f)
    {
        var disc = (Ring.Discriminant(f) * f.LT).Num;

        var n = f.Degree;
        var norm = f.Coefs.Select(e => Double.Abs(e)).Max();
        var normb = f.Coefs.Select(e => BigInteger.Log(BigInteger.Abs(e.Num)) - BigInteger.Log(e.Denom)).Max();
        // return Double.Sqrt(n + 1) * Double.Pow(2, n) * norm;

        var nu = Nu(f);
        if (Logger.Level == LogLevel.Level1)
            Console.WriteLine($"nu = {nu} => {Double.Log(2 * nu)} ~ {Double.Log(n + 1) / 2 + n + normb}");

        var all = IntExt.Primes10000.Where(p => !BigInteger.Remainder(disc, p).IsZero).Take(150)
            .Select(p => (p, s: (int)((Double.Log(n + 1) / 2 + n * Double.Log(2) + normb) / Double.Log(p)) + 1))
            .OrderByDescending(e => e.s).ToArray();

        foreach (var (p, s) in all)
        {
            yield return (p, s);
        }

        // throw new ArgumentException("Prime not found");
    }

    public static KPoly<Rational> ZPoly2QPoly(KPoly<ZnBInt> f) =>
        new(f.x, Rational.KZero(), f.Coefs.Select(e => new Rational(e.K * 2 <= e.Mod ? e.K : e.K - e.Mod)).ToArray());

    public static KPoly<ZnBInt> QPoly2ZnInt(KPoly<Rational> f, Modulus po)
    {
        var coefs = f.Coefs.Select(e => e.ToZnBInt(po)).ToArray();
        return new(f.x, po.Zero, coefs);
    }

    static KPoly<ZnInt> QPoly2ZnInt(KPoly<Rational> f, int p)
    {
        var coefs = f.Coefs.Select(e => e.ToZnInt(p)).ToArray();
        return new(f.x, ZnInt.ZnZero(p), coefs);
    }

    static KPoly<ZnBInt> ZPoly2ZnInt(KPoly<ZnBInt> f, Modulus po) =>
        new(f.x, po.Zero, f.Coefs.Select(e => new ZnBInt(po, e.K)).ToArray());

    public static (KPoly<ZnBInt> f0, KPoly<ZnBInt>[] firr0, KPoly<ZnBInt>[] firr) HenselLifting(KPoly<Rational> f, int p, int o)
    {
        if (f.Coefs.Any(c => !c.Denom.IsOne))
            throw new ArgumentException();

        var ai0 = (new Un(p)).GetGenerators().First()[new(p, 1)];
        var a0 = ZnBInt.ZnZero(p) + ai0.K;
        var po = a0.Details;
        var f0 = QPoly2ZnInt(f, po);

        // GlobalStopWatch.AddLap();
        // var firr0 = Firr(f0, a0).ToArray();
        var firr0 = BerlekampProbabilisticVShoup(f0.Monic, a0).ToArray();
        // var firr0 = BerlekampProbabilisticAECF(f0, a0).ToArray();
        // var firr0 = CantorZassenhausVShoup(f0, a0).ToArray();
        // GlobalStopWatch.Show("Firr");

        var all = new List<KPoly<ZnBInt>>(firr0);
        while (po.O < o && all.Count > 1)
        {
            // ++po;
            po *= 2;
            var tmp = new List<KPoly<ZnBInt>>();

            var fa = QPoly2ZnInt(f, po);
            foreach (var g in all)
            {
                var gi = ZPoly2ZnInt(g, po);
                var y = FG.EPoly(gi);
                var dgi = gi.Derivative.Substitute(y); // new EPoly<ZnBInt>(gi, gi.Derivative.Div(gi).rem);
                var fi = fa.Substitute(y); // new EPoly<ZnBInt>(gi, fa.Div(gi).rem);
                var dfi = fa.Derivative.Substitute(y); // new EPoly<ZnBInt>(gi, fa.Derivative.Div(gi).rem);
                var ri = ((dgi * fi) / dfi).Poly;
                tmp.Add(gi + ri);
            }

            all.Clear();
            all = tmp.ToList();
        }

        return (f0, firr0, all.ToArray());
    }

    public static KPoly<Rational>[] HenselLiftingNaive(KPoly<Rational> f, int p, int o)
    {
        if (f.Coefs.Any(c => !c.Denom.IsOne))
            throw new ArgumentException();

        var (f0, firr0, allS) = HenselLifting(f, p, o);
        var c0 = f.LT * allS[0].KOne;
        var o0 = allS.Max(pl => pl.KZero.Details.O);
        var xp = FG.KPoly(f.x, new ZnBInt(new Modulus(p, o0), 0));
        var F = new KPoly<Rational>(f.x, f.KZero, f.Coefs);

        var listIrr = new List<KPoly<Rational>>();
        var sz = 0;
        var nbCombs = 1;
        var itr = 0;

        if (Logger.Level != LogLevel.Off)
            Console.WriteLine($"######## Candidate Prime P = {p}; Sigma = {o}; P^o={BigInteger.Pow(p, o)}");

        while (allS.Length != sz)
        {
            sz = allS.Length;
            if (Logger.Level == LogLevel.Level2)
                Console.WriteLine($"######## Nb Combinaisons : 2^{sz} = {BigInteger.Pow(2, sz)} ########");

            foreach (var combs in allS.AllCombinationsFromM(nbCombs))
            {
                itr++;
                var H = c0 * combs.Aggregate(xp.One, (acc, a) => a * acc);
                var G = ZPoly2QPoly(H).Primitive();
                if (F.Div(G).rem.IsZero())
                {
                    listIrr.Add(G);
                    F /= G;
                    allS = allS.Except(combs).ToArray();

                    nbCombs = combs.Length;
                    if (Logger.Level == LogLevel.Level2)
                        Console.WriteLine($"@@@@@@@@ itr:{itr} nbCombs:{nbCombs}");
                    
                    break;
                }
            }
        }

        if (allS.Length > 0)
        {
            var Hlast = c0 * allS.Aggregate(xp.One, (acc, a) => a * acc);
            var Glast = ZPoly2QPoly(Hlast).Primitive();
            listIrr.Add(Glast);
        }

        var f1 = listIrr.Aggregate((a, b) => a * b);
        if (!f.LT.Equals(Rational.KOne()))
        {
            firr0 = firr0.Prepend(f.LT * f0.KOne * f0.One).ToArray();
            var a = f.LT / f1.LT;
            if (!a.Equals(a.One))
            {
                listIrr.Insert(0, a * f.One);
            }
        }

        var check = f.Monic.Equals(f1.Monic);
        if (!check)
        {
            if (Logger.Level != LogLevel.Off)
                Console.WriteLine($"@@@@@@@@@@ {f.Monic} <> {f1.Monic} @@@@@@@@@@");
            throw new Exception();
        }

        var listIrr0 = listIrr.OrderBy(q => new KPoly<Rational>(q.x, q.KZero, q.Coefs.Select(Rational.Absolute).ToArray())).ToArray();
        if (Logger.Level != LogLevel.Off)
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

        var Fp = QPoly2ZnInt(F, new Modulus(p, sigma));
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

        Console.WriteLine(mat.T);
        return LLL(mat.T);
        // var lll = ExternLibs.Run_fpLLL(mat);
        // Console.WriteLine(lll);
        // Console.WriteLine();
        // Console.WriteLine(LLL(mat.T));
        // Console.ReadLine();
        // return lll;
    }

    public static KPoly<Rational>[] VanHoeijFactorization(KPoly<Rational> f, int max = 2)
    {
        var discQ = Ring.Discriminant(f).Num;
        if (Logger.Level != LogLevel.Off)
            Console.WriteLine($"f = {f}");
        var discDecomp = IntExt.PrimesDec(discQ);
        if (Logger.Level != LogLevel.Off)
            Console.WriteLine($"Disc(f) = {discQ} ~ {discDecomp.AscendingByKey().GlueMap(" * ", "{0}^{1}")}");

        foreach (var p in IntExt.Primes10000.Where(p => !BigInteger.Remainder(discQ, p).IsZero))
        {
            try
            {
                return SearchVanHoeij(f, p, max).Item3;
            }
            catch (Exception e)
            {
                if (Logger.Level != LogLevel.Off)
                    Console.WriteLine($"#### Prime {p} wont work ####");
            }
        }

        throw new();
    }

    public static (int sigma, int tau, KPoly<Rational>[]) SearchVanHoeij(KPoly<Rational> P, int p, int max = 2)
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
        if (Logger.Level != LogLevel.Off)
            Console.WriteLine(
                $"Search Van-Hoeij Prime P = {p} Tau = {tau} nb = {nb} Theta = {theta} BoundSigma = {boundSigma} MaxSigma = {sigma}; Nu = {nu}");
        for (int i = 2; i <= nb; i++)
        {
            var sigma2 = tau * i;
            try
            {
                var (f0, firr0, irrs2) = HenselLifting(P, p, sigma2);
                var lll = Lattice(P, irrs2, n, p, sigma2, tau);
                var cols = lll.Cols.OrderBy(SquareNorm2).ToArray();
                var rgM = m.Range();
                var rgS = s.Range();
                var t = rgM.Where(i0 => rgM.Where(j => j > i0).All(j => Double.Sqrt(SquareNorm2(cols[j])) > nu * nu)).Min();

                var bs = (t + 1).Range().Select(i0 => cols[i0].T.Extract(0, 1, 0, s)).ToArray();
                var coefs = bs.Select(mat => mat.Select(r => (int)r.Num).ToArray()).ToArray();
                var one0 = irrs2[0].One;
                var c0 = one0 * (one0.KOne * P.LT);
                var polys = coefs.Select(l => rgS.Aggregate(c0, (prod, i0) => prod * irrs2[i0].Pow(l[i0])))
                    .Select(e => ZPoly2QPoly(e).Monic).ToArray();

                var one1 = polys[0].One;
                if (!P.LT.Equals(P.KOne))
                {
                    polys = polys.Prepend(P.LT * P.One).ToArray();
                    firr0 = firr0.Prepend(P.LT * f0.KOne * f0.One).ToArray();
                }

                if (Logger.Level != LogLevel.Off)
                {
                    Console.WriteLine($"       Found Sigma = {sigma2} = {i}*Tau");
                    Console.WriteLine($"f = {f0} mod {p}");
                    Console.WriteLine("Fact(f) = {0} mod {1}", firr0.Glue("*", "({0})"), p);
                    bs.Println("LLL Combinaisons");
                    Console.WriteLine("Fact(f) = {0} in Z[X]", polys.Glue("*", "({0})"));
                    Console.WriteLine($"f = Prod[Fact(f)] : {P.Equals(polys.Aggregate(one1, (prod, e) => prod * e))}");
                    Console.WriteLine();
                }
                return (sigma2, tau, polys);
            }
            catch (Exception e)
            {
                if (Logger.Level != LogLevel.Off)
                    Console.WriteLine($"#### Prime {p} Sigma {sigma2} = {i}*Tau wont work ####");
            }
        }

        throw new();
    }

    public static (KPoly<K> nf, KPoly<K> nx) DeflateMin<K>(KPoly<K> f)
        where K : struct, IFieldElt<K>, IRingElt<K>, IElt<K>
    {
        var pows = f.Coefs.Select((c, i) => (c, i)).Where(e => e.i != 0 && !e.c.IsZero()).ToArray();
        if (pows.Length == 0)
            return (f, f.X);

        var gcd = IntExt.Gcd(pows.Select(e => e.i).ToArray());
        var divs = IntExt.Dividors(gcd).Where(k => k != 1).Append(gcd).Distinct().Order().ToArray();
        if (divs.Length == 0)
            return (f, f.X);

        var div = divs.Min();
        var coefs = f.Coefs.Where((c, i) => i % div == 0).ToArray();
        return (new(f.x, f.KZero, coefs), f.X.Pow(div));
    }

    public static (KPoly<K> nf, KPoly<K> nx) DeflateMax<K>(KPoly<K> f)
        where K : struct, IFieldElt<K>, IRingElt<K>, IElt<K>
    {
        var pows = f.Coefs.Select((c, i) => (c, i)).Where(e => e.i != 0 && !e.c.IsZero()).ToArray();
        var gcd = IntExt.Gcd(pows.Select(e => e.i).ToArray());
        var coefs = f.Coefs.Where((c, i) => i % gcd == 0).ToArray();
        return (new(f.x, f.KZero, coefs), f.X.Pow(gcd));
    }

    public static (KPoly<K> nf, KPoly<K> nx) DeflateByN<K>(KPoly<K> f, int n)
        where K : struct, IFieldElt<K>, IRingElt<K>, IElt<K>
    {
        var (f0, x0) = DeflateMax(f);
        if (x0.Degree % n != 0)
            throw new();

        return (f0.Substitute(f.X.Pow(x0.Degree / n)), f.X.Pow(n));
    }

    static (KPoly<Rational> nf, Rational c) CoefDiv(KPoly<Rational> f0)
    {
        var f = f0.ZPoly();
        var coefs = f.Coefs.Select(c => Rational.Absolute(c)).ToArray();

        var dicoNum = coefs.Select((c, i) => (c.Num, i)).Where(e => !e.Num.IsZero && e.i != 0)
            .Select(e => IntExt.PrimesDec(e.Num).Where(kp => kp.Value >= e.i)
                .SelectMany(kp => Enumerable.Repeat(kp.Key, kp.Value / e.i)))
            .ToArray();

        var inter = dicoNum.Aggregate((a0, b0) => a0.IntersectList(b0).ToArray());
        var resNum = inter.Aggregate(BigInteger.One, (prod, k) => k * prod);
        var c = new Rational(resNum);
        var nf = f0.Substitute(f.X / c);
        return (nf, c);
    }

    static (KPoly<Rational> nf, Rational c) CoefMul(KPoly<Rational> f0)
    {
        var f = f0.ZPoly();
        var deg = f.Degree;
        var coefs = f.Coefs.Select(c => Rational.Absolute(c)).ToArray();

        var dicoNum = coefs.Select((c, i) => (c.Num, i)).Where(e => !e.Num.IsZero && e.i != deg)
            .Select(e => IntExt.PrimesDec(e.Num).Where(kp => kp.Value >= deg - e.i)
                .SelectMany(kp => Enumerable.Repeat(kp.Key, kp.Value / (deg - e.i))))
            .ToArray();

        var inter = dicoNum.Aggregate((a0, b0) => a0.IntersectList(b0).ToArray());
        var resNum = inter.Aggregate(BigInteger.One, (prod, k) => k * prod);
        var c = new Rational(resNum);
        var nf = f0.Substitute(f.X * c);

        return (nf, c);
    }


    public static (KPoly<Rational> nf, Rational c) ConstCoef(KPoly<Rational> f, bool monic = false)
    {
        if (f.Degree <= 1)
            return (f, Rational.KOne());

        var (nf1, c1) = CoefDiv(f);
        var (nf2, c2) = CoefMul(nf1);
        var c = c2 / c1;

        if (Logger.Level != LogLevel.Off)
            if (!c.Equals(c.One))
                Console.WriteLine($"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Simplify f={f.Monic} to nf={nf2.Monic} and c={c}");

        if (monic)
            return (nf2.Monic, c);

        return (nf2, c);
    }

    public static (KPoly<Rational> newP, Rational c) QPoly2ZPoly(KPoly<Rational> f)
    {
        var ct = f.Coefs.Last();
        if (!ct.Equals(f.KOne))
        {
            var ans = QPoly2ZPoly(f.Monic);
            return (ct * ans.newP, ans.c);
        }

        var n = f.Degree;
        var denoms = f.Coefs.Select(e => e.Denom).Distinct().ToArray();
        if (denoms.Length == 0)
            denoms = new[] { BigInteger.One, };

        var lcm = IntExt.LcmBigInt(denoms);
        var c0 = IntExt.PrimesDecompositionBigInt(BigInteger.Abs(lcm)).GroupBy(e => e).ToDictionary(e => e.Key, e => e.Count());
        var c1 = c0.ToDictionary(e => e.Key, e => e.Value % n == 0 ? e.Value / n : e.Value)
            .Aggregate(BigInteger.One, (prod, e) => prod * BigInteger.Pow(e.Key, e.Value));
        var c2 = new Rational(c1);

        var f2 = f.Substitute(f.X / c2);
        return (f2, c2);
    }

    public static bool EisensteinCriterion(KPoly<Rational> f)
    {
        var a0 = f[0];
        var an = f.Coefs.Last();
        var ai = f.Coefs.SkipLast(1).Where(ai => !ai.IsZero()).ToArray();
        var pi = ai.Length == 0
            ? Array.Empty<int>()
            : ai.Select(a => IntExt.PrimesDecompositionBigInt(BigInteger.Abs(a.Num)).Distinct())
                .Aggregate((a, b) => a.Intersect(b)).ToArray();

        if (Logger.Level != LogLevel.Off)
            if (pi.Length > 0)
                Console.WriteLine($"Common Primes {pi.Glue("; ")}");

        if (pi.Length == 0)
            return false;
        else
        {
            foreach (var p in pi)
            {
                if (an % p != 0 && a0 % (p * p) != 0)
                {
                    if (Logger.Level != LogLevel.Off)
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

    public static KPoly<Rational>[] FirrZ2(KPoly<Rational> f)
    {
        // if (f.Coefs.Any(c => !c.Denom.IsOne))
        //     throw new($"f isnt in Z[X] : {f}");

        if (f.Degree == 1)
            return new[] { f };

        var (f0, x0) = DeflateMin(f);
        if (x0.Degree == 1)
        {
            return FirrZ(f0);
        }
        else if (f0.Degree == 1)
        {
            return FirrZ(f);
        }
        else
        {
            var facts = FirrZ(f0);
            if (facts.Length == 1)
                return FirrZ(f);

            var facts2 = facts.SelectMany(f1 => FirrZ2(f1.Substitute(x0))).ToArray();
            if (Logger.Level != LogLevel.Off)
            {
                Console.WriteLine($"f = {f}");
                Console.WriteLine("Fact(f) = {0} in Z[X]", facts2.Glue("*", "({0})"));
                Console.WriteLine();
            }

            return facts2;
        }
    }

    public static KPoly<Rational>[] FirrZ(KPoly<Rational> P)
    {
        if (P.Degree <= 1)
            return new[] { P };

        var (f0, c) = ConstCoef(P);
        var f = f0.PrimitiveZPoly();
        var discQ = Ring.Discriminant(f).Num;
        var discDecomp = IntExt.PrimesDec(discQ);
        if (Logger.Level != LogLevel.Off)
        {
            Console.WriteLine($"f = {f}");
            Console.WriteLine($"Disc(f) = {discQ} ~ {discDecomp.AscendingByKey().GlueMap(" * ", "{0}^{1}")}");
        }

        if (f.Degree == 1)
            return new[] { P };

        // if (EisensteinCriterion(f, details))
        // {
        //     return new[] { P };
        // }

        foreach (var (p, o) in PSigma(f))
        {
            try
            {
                var polys = HenselLiftingNaive(f, p, o)
                    .Select(f2 => f2.Degree == 0 ? P.LT * P.One : f2.Substitute(f2.X / c).Monic)
                    .OrderBy(f2 => f2.Degree)
                    .ThenBy(f2 => f2)
                    .ToArray();

                if (polys.All(f2 => f2.Degree != 0) && !P.LT.Equals(P.KOne))
                    polys = polys.Prepend(P.LT * P.One).ToArray();

                if (Logger.Level != LogLevel.Off)
                {
                    Console.WriteLine($"P({P.x}) = {P}");
                    Console.WriteLine($"Fact(P({P.x})) = {{0}} in Q[X]", polys.Glue("*", "({0})"));
                    Console.WriteLine();
                }

                return polys;
                // return SearchVanHoeij(f, p, 2).Item3.Select(f0 => f0.Substitute(f0.X / c).Monic).ToArray();
            }
            catch (Exception)
            {
                if (Logger.Level != LogLevel.Off)
                    Console.WriteLine($"#### Prime {p} and Sigma {o} wont work ####");
            }
        }

        return new[] { P };
    }

    static IEnumerable<(KPoly<Rational>, int)> FactorsMul(KPoly<Rational> P)
    {
        if (P.Degree <= 1)
        {
            yield return (P, 1);
            yield break;
        }

        foreach (var (p0, _, i0) in YunSFF(P))
        {
            foreach (var p2 in FirrZ2(p0))
            {
                yield return (p2.SubstituteChar(P.x), i0);
            }
        }
    }

    public static (KPoly<Rational>, int)[] FactorsQ(KPoly<Rational> P)
    {
        var list0 = FactorsMul(P).ToList();
        var res = list0.Where(e => !e.Item1.Equals(P.One))
            .OrderBy(e => e.Item1.Degree)
            .ThenBy(e => e.Item2)
            .ThenBy(e => e.Item1).ToArray();

        if (Logger.Level != LogLevel.Off)
        {
            string Fmt((KPoly<Rational>, int) e) => e.Item2 == 1 ? $"({e.Item1})" : $"({e.Item1})^{e.Item2} ";
            Console.WriteLine($"f0 = {P}");
            Console.WriteLine("Fact(f0) = {0} in Q[X]", res.Select(Fmt).Glue("*"));
            Console.WriteLine();
        }

        return res;
    }

    public static KPoly<EPoly<ZnInt>>[] FirrFq(KPoly<Rational> f, int q)
    {
        var dec = IntExt.PrimesDec(q);
        if (dec.Count != 1)
            throw new($"q = {q} isnt prime or power of prime");

        var (p, n) = dec.First();
        var fq = new Fq(q, 'a');
        var f0 = new KPoly<EPoly<ZnInt>>(f.x, fq.Zero, f.Coefs.Select(c => c.ToZnInt(p) * fq.One).ToArray());
        return BerlekampProbabilisticAECF(f0, fq.One.X).ToArray();
    }

    public static KPoly<ZnBInt>[] FirrFp(KPoly<Rational> f, int p)
    {
        if (!IntExt.Primes10000.Contains(p))
            throw new($"p = {p} isnt prime");

        var e = IntExt.Solve_k_pow_m_equal_one_mod_n_strict(p, p - 1);
        var a0 = new ZnBInt(p, e, 1);
        var po = a0.Details;
        var f0 = QPoly2ZnInt(f, po);
        return BerlekampProbabilisticAECF(f0, a0).ToArray();
    }
}
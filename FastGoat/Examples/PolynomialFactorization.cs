using System.Numerics;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Polynoms;

namespace FastGoat.Examples;

public static class PolynomialFactorization
{

    static KPoly<K> RandPoly<K>(K scalar, int p, int n) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var coefs = n.Range().Select(i => IntExt.Rng.Next(-p, p + 1) * scalar.One).TrimSeq().Append(scalar.One).ToArray();
        return new KPoly<K>('x', scalar, coefs);
    }

    public static KPoly<K> RandPoly<K>(K scalar, int p, int[] degrees)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        Dictionary<int, int> maxSep = new() { [2] = 3, [3] = 2, [5] = 1 };

        KPoly<K> RandPolyDegSep(K s0, int p0, int n0)
        {
            var f = RandPoly(s0, p0, n0);
            if (maxSep.ContainsKey(p0) && n0 <= 3 && IntExt.Rng.Next(2) == 0)
                f = f.Substitute(f.X.Pow(p * IntExt.Rng.Next(1, maxSep[p0] + 1)));

            return f;
        }

        while (true)
        {
            var f = degrees.Select(n => RandPolyDegSep(scalar, p, n)).Aggregate((a, b) => a * b).Monic;
            if (f.Degree > 1)
                return f;
        }
    }

    public static KPoly<K> RandPolySep<K>(K scalar, int p, int n) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        while (true)
        {
            var f = RandPoly(scalar, p, n);
            if (f.Degree > 1 && !Ring.Discriminant(f).IsZero())
                return f.Monic;
        }
    }

    public static KPoly<K> RandPolySepStrict<K>(K scalar, int p, int n) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        while (true)
        {
            var f = RandPoly(scalar, p, n);
            if (f.Degree == n && !Ring.Discriminant(f).IsZero())
                return f.Monic;
        }
    }

    static KPoly<ZnInt> ProdIrr(int p, int d)
    {
        var x = FG.ZPoly(p);
        return x.Pow(p.Pow(d)) - x;
    }

    static bool IsIrreductibleFp(KPoly<ZnInt> f)
    {
        if (f.Degree < 1)
            return true;

        var p = f.P;
        var n = f.Degree;
        var divs = IntExt.Dividors(n);
        return ProdIrr(p, n).Div(f).rem.IsZero() && divs.All(d => !ProdIrr(p, d).Div(f).rem.IsZero());
    }

    public static void IrreductibleRandPolys()
    {
        var n = 5;
        var p = 5;
        for (int i = 0; i < 4 * n; i++)
        {
            var f = RandPoly(ZnInt.KZero(p), p, n);
            Console.WriteLine($"{{0,{-7 * (n + 1)}}} is irreductible : {{1}}", f, IsIrreductibleFp(f));
        }
    }

    public static List<(KPoly<K> g, int q, int i)> MusserSFF<K>(KPoly<K> f)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var L = new List<(KPoly<K>, int, int)>();
        var c = Ring.Gcd(f, f.Derivative);
        var i = 1;
        var g = f / c;
        while (g.Degree >= 1)
        {
            var p = Ring.Gcd(c, g);
            c = c / p;
            if (g.Degree > p.Degree)
                L.Add(((g / p).Monic, 1, i));

            g = p;
            ++i;
        }

        return L;
    }

    public static List<(KPoly<K> g, int q, int i)> YunSFF<K>(KPoly<K> f)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var L = new List<(KPoly<K>, int, int)>();
        var l = 1;
        var df = f.Derivative;
        var u = Ring.Gcd(f, df);
        var v = f / u;
        var w = f.Derivative / u;
        while (v.Degree >= 1)
        {
            var w_dv = w - v.Derivative;
            var h = Ring.Gcd(v, w_dv);
            w = w_dv / h;
            v = v / h;
            if (h.Degree >= 1)
                L.Add((h.Monic, 1, l));

            ++l;
        }

        return L;
    }

    public static List<(KPoly<K> g, int q, int i)> YunSFFDetails<K>(KPoly<K> f)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var sep = YunSFF(f);
        Console.WriteLine($"Square Free Factors of F = {f}");
        Console.WriteLine(sep.Glue("\n"));
        Console.WriteLine();
        return sep;
    }

    static (KPoly<K> c, int p) DeflateP<K>(KPoly<K> f) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var p = f.KZero.P;
        if (p == 0)
            throw new ArgumentException();

        var d = f.Degree;
        var decomp = IntExt.PrimesDecomposition(d).ToArray();
        if (!decomp.Contains(p))
            return (f, 1);

        var q = decomp.Count(i => i == p);
        var xi = q.Range(1).Select(i => p.Pow(i)).ToArray();
        var dq = (d + 1).Range();
        if (dq.Where(j => j % p != 0).All(j => f[j].IsZero()))
        {
            var coefs = dq.Where(j => j % p == 0).Select(j => f[j]).ToArray();
            return (new KPoly<K>(f.x, f.KZero, coefs), p);
        }

        return (f, 1);
    }

    // Page 334, Book ACFE, 18.4 Factorisation séparable
    static IEnumerable<(KPoly<K> g, int q, int m)> GianniTrager<K>(KPoly<K> f, int q = 1)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        if (f.Degree != 0)
        {
            var sff = MusserSFF(f);
            var l0 = sff.Select(e => (e.g, q, m: e.i)).ToList();
            foreach (var l in l0)
            {
                yield return l;
            }

            var gi = sff.Aggregate(f.One, (acc, a) => acc * a.g.Pow(a.i));
            var c = f / gi;
            var cf = DeflateP(c);
            if (cf.p != 1)
            {
                foreach (var l in GianniTrager(cf.c, cf.p * q))
                {
                    yield return l;
                }
            }
        }
    }

    static void CheckSeparability<K>(KPoly<K> f, (KPoly<K> g, int q, int m)[] fSep)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var x = f.X;
        var fSep2 = fSep.Select(e => (gxp: e.g.Substitute(x.Pow(e.q)), e.q, e.m)).ToArray();
        var f0 = fSep2.Aggregate(x.One, (acc, e) => acc * e.gxp.Pow(e.m));
        if (f.Equals(f0))
        {
            Console.WriteLine("Prop (S1) pass");
        }
        else
        {
            Console.WriteLine("Prop (S1) fail");
            return;
        }

        var tuples = fSep2.Grid2D(fSep2).Where(e => !e.t1.gxp.Equals(e.t2.gxp)).ToArray();
        if (tuples.All(e => Ring.Gcd(e.t1.gxp, e.t2.gxp).Monic.Equals(f.One)))
        {
            Console.WriteLine("Prop (S2) pass");
        }
        else
        {
            Console.WriteLine("Prop (S2) fail");
            var pb = tuples.First(e => !Ring.Gcd(e.t1.gxp, e.t2.gxp).Monic.Equals(f.One));
            Console.WriteLine(pb);
            return;
        }

        var p = f.P;
        if (p == 0 || fSep2.All(e => e.m % p != 0))
        {
            Console.WriteLine("Prop (S3) pass");
        }
        else
        {
            Console.WriteLine("Prop (S3) fail");
            return;
        }

        if (fSep.All(e => e.g.Degree >= 1 && !Ring.Discriminant(e.g).IsZero()))
        {
            Console.WriteLine("Prop (S4) pass");
        }
        else
        {
            Console.WriteLine("Prop (S4) fail");
            var pb = fSep.First(e => e.g.Degree == 0 || Ring.Discriminant(e.g).IsZero());
            Console.WriteLine(pb);
            Console.WriteLine(Ring.Discriminant(pb.g));
            return;
        }

        if (tuples.All(e => !(e.t1.q, e.t1.m).Equals((e.t2.q, e.t2.m))))
        {
            Console.WriteLine("Prop (S5) pass");
        }
        else
        {
            Console.WriteLine("Prop (S5) fail");
            return;
        }

        Console.WriteLine("Successful Factorization");
        Console.WriteLine();
    }

    static void CheckIrreductibility((KPoly<ZnInt> g, int q, int m)[] fSep)
    {
        var digits = fSep.Max(e => e.g.ToString().Length);
        var fmt = $"{{0,-{digits}}} IsIrreductible : {{1}}";
        foreach (var e in fSep)
            Console.WriteLine(fmt, e.g, IsIrreductibleFp(e.g));

        Console.WriteLine();
    }

    static EPoly<K>[] CanonicalBase<K>(KPoly<K> f) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var x = new EPoly<K>(f);
        return f.Degree.Range().Select(i => x.Pow(i)).ToArray();
    }

    static KMatrix<K> BerlekampMatrix<K>(KPoly<K> f,  EPoly<K>[] baseCan, int q)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var n = baseCan.Length;
        var M = new K[n, n];
        var polys = baseCan.Select(g => g.Pow(q) - g).ToArray();
        foreach (var (i, j) in n.Range().Grid2D(n.Range()))
        {
            M[i, j] = polys[j][i];
        }

        return new(M);
    }

    public static KPoly<K>[] FrobeniusKernel<K>(KPoly<K> f, int q)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var baseCan = CanonicalBase(f);
        var bm = BerlekampMatrix(f, baseCan, q);
        var (nt, ns) = bm.NullSpace();
        var (m, n) = ns.Dim;
        var polys = n.Range().Select(j => m.Range().Select(i => ns[i, j] * baseCan[i]).Aggregate((a, b) => a + b).Poly)
            .ToArray();
        return polys;
    }

    static IEnumerable<KPoly<K>> FirrInternal<K>(KPoly<K> f, List<K> allF)where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var polys = FrobeniusKernel(f, allF.Count);
        if (polys.Length > 1)
        {
            foreach (var (g, a) in polys.Where(g => g.Degree > 0).Grid2D(allF.Skip(1)))
            {
                var g_a = g - a;
                var gcd = Ring.Gcd(f, g_a).Monic;
                if (gcd.Degree != 0)
                {
                    foreach (var f1 in FirrInternal(gcd, allF))
                        yield return f1;

                    foreach (var f1 in FirrInternal(f / gcd, allF))
                        yield return f1;

                    break;
                }
            }
        }
        else
        {
            yield return f;
        }
    }

    public static IEnumerable<KPoly<K>> Firr<K>(KPoly<K> f, K a0)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var acc = a0.One;
        var allF = new List<K>() { a0.Zero };
        do
        {
            allF.Add(acc);
            acc *= a0;
        } while (!acc.Equals(a0.One));
        return FirrInternal(f, allF);
    }

    public static List<(KPoly<K> g, int q, int m)> FirrFsep<K>(KPoly<K> f, K a0)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        List<(KPoly<K> g, int q, int m)> all = new();
        foreach (var (g, q, m) in GianniTrager(f))
            all.AddRange(Firr(g, a0).Select(g0 => (g0, q, m)));

        return all;
    }

    public static void DisplayFactorization<K>(KPoly<K> f, K a0)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var p = f.P;
        var x = f.X;
        Console.WriteLine($"f = {f} mod ({p})");
        Console.WriteLine($"Disc(f) = {Ring.Discriminant(f)} mod ({p})");
        var firrm = FirrFsep(f, a0);
        List<(KPoly<K> g, int m)> all = new();
        foreach (var (g, q, m) in firrm)
        {
            var xq = x.Pow(q);
            all.Add((g.Substitute(xq), m));
        }

        string Display(KPoly<K> g, int m) => m == 1 ? $"({g})" : $"({g})^{m}";
        var prod = all.Aggregate(x.One, (acc, gm) => acc * gm.g.Pow(gm.m));
        var fact = all.OrderBy(e => e.g.Pow(e.m)).Select(e => Display(e.g, e.m)).Glue(" * ");
        Console.WriteLine($"Decomposition = {firrm.Glue(", ", "{{{0}}}")}");
        Console.WriteLine($"Fact(f) = {fact} mod ({p})");
        Console.WriteLine($"f = Fact(f) : {prod.Equals(f)}");

        Console.WriteLine();
    }

    public static Rational SquareNorm2(KMatrix<Rational> v)
    {
        if (v.M == 1)
            return (v * v.T)[0, 0];
        else if (v.N == 1)
            return (v.T * v)[0, 0];

        throw new ArgumentException();
    }

    static void SwapRows<T>(int i, int j, T[,] A)
    {
        var cols = A.GetLength(1);
        for (int k = 0; k < cols; k++)
            (A[i, k], A[j, k]) = (A[j, k], A[i, k]);
    }

    static KMatrix<Rational> RandMatrixRational(int n)
    {
        var m = new KMatrix<Rational>(Rational.KZero(), n, n);
        while (true)
        {
            var m0 = Ring.Matrix(n, Rational.KZero(), (n * n).Range().Select(i => IntExt.Rng.Next(n + 1)).ToArray());
            m = new(m0);
            if (!m.Det.IsZero())
                return m;
        }
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

    public static void SquareFreeFactorizationQ()
    {
        Monom.Display = MonomDisplay.Superscript;
        var x = FG.QPoly();

        {
            var f = (x + 1) * (x + 2) * (x + 3).Pow(2) * (x + 4).Pow(2) * (x + 5).Pow(3) * (x + 6).Pow(3);
            Console.WriteLine(f);
            Console.WriteLine("Musser Algo");
            var sff1 = MusserSFF(f).ToArray();
            Console.WriteLine(sff1.Glue("\n"));
            CheckSeparability(f, sff1);

            Console.WriteLine("Yun Algo");
            var sff2 = YunSFF(f).ToArray();
            Console.WriteLine(sff2.Glue("\n"));
            CheckSeparability(f, sff2);
        }

        {
            var f = x.Pow(3) * (x + 2).Pow(4) * (x.Pow(2) + 2 * x - 2).Pow(2) * (x.Pow(3) + 5);
            Console.WriteLine(f);
            Console.WriteLine("Musser Algo");
            var sff1 = MusserSFF(f).ToArray();
            Console.WriteLine(sff1.Glue("\n"));
            CheckSeparability(f, sff1);

            Console.WriteLine("Yun Algo");
            var sff2 = YunSFF(f).ToArray();
            Console.WriteLine(sff2.Glue("\n"));
            CheckSeparability(f, sff2);
        }
    }

    public static void SeparableFactorizationFp()
    {
        Monom.Display = MonomDisplay.Superscript;

        {
            var x = FG.ZPoly(3);
            var f = x.Pow(2) * (x + 1).Pow(3) * (x + 2).Pow(4);
            Console.WriteLine(f);
            var l = GianniTrager(f).ToArray();
            Console.WriteLine(l.Glue("\n"));
            CheckSeparability(f, l);
            CheckIrreductibility(l);
            Console.WriteLine(Ring.Discriminant(f));
        }

        {
            var x = FG.ZPoly(2);
            var f = x.Pow(2) * (x + 1).Pow(3) * (x.Pow(2) + 1).Pow(4);
            Console.WriteLine(f);
            var l = GianniTrager(f).ToArray();
            Console.WriteLine(l.Glue("\n"));
            CheckSeparability(f, l);
            CheckIrreductibility(l);
        }

        {
            var (x, t) = FG.FpT_Poly(3);
            // F = (X + 2T)7 (X3 + 2T)3 (X6 + T)
            var f = (x + 2 * t).Pow(7) * (x.Pow(3) + 2 * t).Pow(3) * (x.Pow(6) + t);
            Console.WriteLine(f);
            var l = GianniTrager(f).ToArray();
            Console.WriteLine(l.Glue("\n"));
            CheckSeparability(f, l);
        }

        {
            Monom.Display = MonomDisplay.StarPowFct;
            var x = FG.ZPoly(3);
            // x15 + 2x14 + 2x12 + x11 + 2x10 + 2x8 + x7 + 2x6 + 2x4
            var f = x.Pow(15) + 2 * x.Pow(14) + 2 * x.Pow(12) + x.Pow(11) + 2 * x.Pow(10) + 2 * x.Pow(8) + x.Pow(7) +
                    2 * x.Pow(6) + 2 * x.Pow(4);
            Console.WriteLine(f);
            var l = GianniTrager(f).ToArray();
            Console.WriteLine(l.Glue("\n"));
            CheckSeparability(f, l);
            CheckIrreductibility(l);
        }
    }

    public static void Separable2IrreductibleFp()
    {
        Monom.Display = MonomDisplay.StarCaret;

        for (int i = 0; i < 20; ++i)
        {
            var p = IntExt.Primes10000[IntExt.Rng.Next(10)]; // 2, 3, 5, 7, 11, 13, 17, 19, 23, 29
            var d = IntExt.Rng.Next((int)(Math.Log(50) / Math.Log(p))) + 1; // p^d < 30 => 4, 8, 9, 16, 25, 27, 32, 49
            var fq = new Fq(p.Pow(d), 'a');
            var gf = FG.Galois(p.Pow(d), 'a');
            var a0 = gf.GetGenerators().First();
            Console.WriteLine($"{fq} with {fq.F} = 0");
            var n = 2 + IntExt.Rng.Next(7);
            var f = RandPolySep(fq.One, p, n);
            Console.WriteLine($"f = {f} mod ({p})");

            var firr = Firr(f, a0).Order().ToArray();
            var prod = firr.Aggregate((a, b) => a * b);

            Console.WriteLine($"Disc(f) = {Ring.Discriminant(f)} mod ({p})");
            Console.WriteLine($"Fact(f) = {firr.Glue("*", "({0})")} mod ({p})");
            if (firr.Length > 1)
                Console.WriteLine($"f = Fact(f) : {prod.Equals(f)}");
            else
                Console.WriteLine("f is irreductible");

            Console.WriteLine();
        }
    }

    public static void FactorizationFp()
    {
        Monom.Display = MonomDisplay.StarCaret;

        for (int j = 0; j < 20; j++)
        {
            var p = IntExt.Primes10000[IntExt.Rng.Next(5)]; // 2, 3, 5, 7, 11
            var a0 = new Un(p).GetGenerators().First()[new ZnInt(p, 1)];
            var n = 2 + IntExt.Rng.Next(11);
            var degrees = IntExt.Partitions32[n].Where(l => l.All(i => i != 1)).OrderBy(i => IntExt.Rng.NextDouble()).First()
                .ToArray();
            var f = RandPoly(ZnInt.KZero(p), p, degrees);
            DisplayFactorization(f, a0);
        }
    }

    public static void FactorizationFp2()
    {
        Monom.Display = MonomDisplay.StarCaret;

        {
            var p = 3;
            var a0 = new Un(p).GetGenerators().First()[new ZnInt(p, 1)];
            var x = FG.ZPoly(p);
            var f = x.Pow(2) * (x + 1).Pow(3) * (x + 2).Pow(4);
            DisplayFactorization(f, a0);
        }

        {
            var p = 2;
            var a0 = new Un(p).GetGenerators().First()[new ZnInt(p, 1)];
            var x = FG.ZPoly(p);
            var f = x.Pow(2) * (x + 1).Pow(3) * (x.Pow(2) + 1).Pow(4);
            DisplayFactorization(f, a0);
        }

        {
            var p = 3;
            var a0 = new Un(p).GetGenerators().First()[new ZnInt(p, 1)];
            var x = FG.ZPoly(p);
            var f = x.Pow(15) + 2 * x.Pow(14) + 2 * x.Pow(12) + x.Pow(11) + 2 * x.Pow(10) + 2 * x.Pow(8) + x.Pow(7) +
                    2 * x.Pow(6) + 2 * x.Pow(4);
            DisplayFactorization(f, a0);
        }
    }
}
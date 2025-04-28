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
    static PolynomialFactorization()
    {
        Ring.DisplayPolynomial = MonomDisplay.StarCaret;
    }

    public static KPoly<K> RandPoly<K>(K scalar, int p, int n, bool monic = true) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var coefs = n.Range().Select(i => IntExt.Rng.Next(-p, p + 1) * scalar.One).ToList();
        if (monic)
            coefs.Add(scalar.One);
        else
            coefs.Add(IntExt.Rng.Next(-p, p + 1) * scalar.One);

        return new KPoly<K>('x', scalar, coefs.TrimSeq().ToArray());
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

    public static KPoly<K> RandPolyStrict<K>(K scalar, int p, int[] degrees)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        while (true)
        {
            var f = degrees.Select(n => RandPolySepStrict(scalar, p, n)).Aggregate((a, b) => a * b).Monic;
            if (f.Degree > 1)
                return f;
        }
    }

    public static KPoly<K> RandPolySep<K>(K scalar, int p, int n, bool monic = true)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        while (true)
        {
            var f = RandPoly(scalar, p, n, monic);
            if (f.Degree > 1 && !Ring.Discriminant(f).IsZero())
                return f;
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
            var f = RandPoly(ZnInt.ZnZero(p), p, n);
            Console.WriteLine($"{{0,{-7 * (n + 1)}}} is irreductible : {{1}}", f, IsIrreductibleFp(f));
        }
    }

    public static List<(KPoly<K> g, int q, int i)> YunSFFDetails<K>(KPoly<K> f)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var sep = IntFactorisation.YunSFF(f);
        Console.WriteLine($"Square Free Factors of F = {f}");
        Console.WriteLine(sep.Glue("\n"));
        Console.WriteLine();
        return sep;
    }

    static void CheckSeparability<K>(KPoly<K> f, (KPoly<K> g, int q, int m)[] fSep)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var x = f.X;
        var fSep2 = fSep.Where(e => e.g.Degree != 0).Select(e => (gxp: e.g.Substitute(x.Pow(e.q)), e.q, e.m)).ToArray();
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

        if (fSep.Where(e => e.g.Degree != 0).All(e => e.g.Degree >= 1 && !Ring.Discriminant(e.g).IsZero()))
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

    public static void DisplayFactorization<K>(KPoly<K> f, K a0)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var p = f.P;
        var x = f.X;
        Console.WriteLine($"f = {f} mod ({p})");
        Console.WriteLine($"Disc(f) = {Ring.Discriminant(f)} mod ({p})");
        var firrm = IntFactorisation.FirrFsep(f, a0);
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

    public static void SquareFreeFactorizationQ()
    {
        Ring.DisplayPolynomial = MonomDisplay.Superscript;
        var x = FG.QPoly();

        {
            var f = (x + 1) * (x + 2) * (x + 3).Pow(2) * (x + 4).Pow(2) * (x + 5).Pow(3) * (x + 6).Pow(3);
            Console.WriteLine(f);
            Console.WriteLine("Musser Algo");
            var sff1 = IntFactorisation.MusserSFF(f).ToArray();
            Console.WriteLine(sff1.Glue("\n"));
            CheckSeparability(f, sff1);

            Console.WriteLine("Yun Algo");
            var sff2 = IntFactorisation.YunSFF(f).ToArray();
            Console.WriteLine(sff2.Glue("\n"));
            CheckSeparability(f, sff2);
        }

        {
            var f = x.Pow(3) * (x + 2).Pow(4) * (x.Pow(2) + 2 * x - 2).Pow(2) * (x.Pow(3) + 5);
            Console.WriteLine(f);
            Console.WriteLine("Musser Algo");
            var sff1 = IntFactorisation.MusserSFF(f).ToArray();
            Console.WriteLine(sff1.Glue("\n"));
            CheckSeparability(f, sff1);

            Console.WriteLine("Yun Algo");
            var sff2 = IntFactorisation.YunSFF(f).ToArray();
            Console.WriteLine(sff2.Glue("\n"));
            CheckSeparability(f, sff2);
        }
    }

    public static void SeparableFactorizationFp()
    {
        Ring.DisplayPolynomial = MonomDisplay.Superscript;

        {
            var x = FG.ZPoly(3);
            var f = x.Pow(2) * (x + 1).Pow(3) * (x + 2).Pow(4);
            Console.WriteLine(f);
            var l = IntFactorisation.GianniTrager(f).ToArray();
            Console.WriteLine(l.Glue("\n"));
            CheckSeparability(f, l);
            CheckIrreductibility(l);
            Console.WriteLine(Ring.Discriminant(f));
        }

        {
            var x = FG.ZPoly(2);
            var f = x.Pow(2) * (x + 1).Pow(3) * (x.Pow(2) + 1).Pow(4);
            Console.WriteLine(f);
            var l = IntFactorisation.GianniTrager(f).ToArray();
            Console.WriteLine(l.Glue("\n"));
            CheckSeparability(f, l);
            CheckIrreductibility(l);
        }

        {
            var (x, t) = FG.FpT_Poly(3);
            // F = (X + 2T)7 (X3 + 2T)3 (X6 + T)
            var f = (x + 2 * t).Pow(7) * (x.Pow(3) + 2 * t).Pow(3) * (x.Pow(6) + t);
            Console.WriteLine(f);
            var l = IntFactorisation.GianniTrager(f).ToArray();
            Console.WriteLine(l.Glue("\n"));
            CheckSeparability(f, l);
        }

        {
            var x = FG.ZPoly(3);
            // x15 + 2x14 + 2x12 + x11 + 2x10 + 2x8 + x7 + 2x6 + 2x4
            var f = x.Pow(15) + 2 * x.Pow(14) + 2 * x.Pow(12) + x.Pow(11) + 2 * x.Pow(10) + 2 * x.Pow(8) + x.Pow(7) +
                    2 * x.Pow(6) + 2 * x.Pow(4);
            Console.WriteLine(f);
            var l = IntFactorisation.GianniTrager(f).ToArray();
            Console.WriteLine(l.Glue("\n"));
            CheckSeparability(f, l);
            CheckIrreductibility(l);
        }
    }

    public static void Separable2IrreductibleFp()
    {
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

            var firr = IntFactorisation.Firr(f, a0).Order().ToArray();
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
        for (int j = 0; j < 20; j++)
        {
            var p = IntExt.Primes10000[IntExt.Rng.Next(5)]; // 2, 3, 5, 7, 11
            var a0 = new Un(p).GetGenerators().First()[new ZnInt(p, 1)];
            var n = 2 + IntExt.Rng.Next(11);
            var degrees = IntExt.Partitions32[n].Where(l => l.All(i => i != 1)).OrderBy(i => IntExt.Rng.NextDouble()).First()
                .ToArray();
            var f = RandPoly(ZnInt.ZnZero(p), p, degrees);
            DisplayFactorization(f, a0);
        }
    }

    public static void FactorizationFp2()
    {
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

    public static void BenchFactorisationFiniteFields()
    {
        for (int i = 0; i < 20; ++i)
        {
            var p = IntExt.Primes10000[IntExt.Rng.Next(5)]; // 2, 3, 5, 7, 11
            var d = IntExt.Rng.Next((int)(Math.Log(50) / Math.Log(p))) + 1; // p^d < 50 => 4, 8, 9, 16, 25, 27, 32, 49
            (p, d) = (2, 1);
            var fq = new Fq(p.Pow(d), 'a');
            var gf = FG.Galois(p.Pow(d), 'a');
            var a0 = gf.GetGenerators().First();
            var n = 2 + IntExt.Rng.Next(6);
            var f0 = RandPolySep(fq.One, p, n);
            var f1 = RandPolySep(fq.One, p, n);
            var f2 = RandPolySep(fq.One, p, n);
            var f = f0 * f1 * f2;
            if (Ring.Discriminant(f).IsZero())
            {
                --i;
                continue;
            }

            Console.WriteLine($"Test[{i + 1}]");
            Console.WriteLine($"{fq} with {fq.F} = 0");
            Console.WriteLine($"f = {f} mod ({p})");
            Console.WriteLine($"Disc(f) = {Ring.Discriminant(f)} mod ({p})");

            var firr0 = IntFactorisation.Firr(f, a0).Order().ToArray();
            Console.WriteLine($"Fact1(f) = {firr0.Glue("*", "({0})")} mod ({p})");

            var firr1 = IntFactorisation.CantorZassenhausVShoup(f, a0, fq.Q).Order().ToArray();
            Console.WriteLine($"Fact2(f) = {firr1.Glue("*", "({0})")} mod ({p})");

            var firr2 = IntFactorisation.CantorZassenhausAECF(f, a0, fq.Q).Order().ToArray();
            Console.WriteLine($"Fact3(f) = {firr2.Glue("*", "({0})")} mod ({p})");

            var firr3 = IntFactorisation.BerlekampProbabilisticVShoup(f, a0, fq.Q).Order().ToArray();
            Console.WriteLine($"Fact4(f) = {firr3.Glue("*", "({0})")} mod ({p})");

            var firr4 = IntFactorisation.BerlekampProbabilisticAECF(f, a0, fq.Q).Order().ToArray();
            Console.WriteLine($"Fact5(f) = {firr4.Glue("*", "({0})")} mod ({p})");

            var check1 = firr0.Aggregate(f.One, (prod, fi) => fi * prod).Equals(f);
            var check2 = firr0.SequenceEqual(firr1);
            var check3 = firr0.SequenceEqual(firr2);
            var check4 = firr0.SequenceEqual(firr3);
            var check5 = firr0.SequenceEqual(firr4);
            Console.WriteLine($"Check1 : {check1}");
            Console.WriteLine($"Check2 : {check2}");
            Console.WriteLine($"Check3 : {check3}");
            Console.WriteLine($"Check4 : {check4}");
            Console.WriteLine($"Check5 : {check5}");

            var nb = 3;
            GlobalStopWatch.Bench(nb, "B1", () => IntFactorisation.Firr(f, a0).Order().ToArray());
            GlobalStopWatch.Bench(nb, "B2", () => IntFactorisation.BerlekampProbabilisticVShoup(f, a0, fq.Q).Order().ToArray());
            GlobalStopWatch.Bench(nb, "B3", () => IntFactorisation.BerlekampProbabilisticAECF(f, a0, fq.Q).Order().ToArray());
            GlobalStopWatch.Bench(nb, "CZ2", () => IntFactorisation.CantorZassenhausVShoup(f, a0, fq.Q).Order().ToArray());
            GlobalStopWatch.Bench(nb, "CZ3", () => IntFactorisation.CantorZassenhausAECF(f, a0, fq.Q).Order().ToArray());

            Console.WriteLine();
        }
    }
}
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Integers;

namespace FastGoat.Examples;

public static class AlgebraicFactorization
{
    static EPoly<K>[] Ebase<K>(EPoly<K> scalar) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var n = scalar.F.Degree;
        return n.Range().Select(i => scalar.X.Pow(i)).ToArray();
    }

    static KMatrix<K> Endo<K>(EPoly<K> b) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var bs = Ebase(b);
        var vs = new KMatrix<K>[bs.Length];
        for (int i = 0; i < bs.Length; i++)
        {
            var e = bs[i];
            var f = e * b;
            var m0 = new K[bs.Length, 1];
            for (int j = 0; j < bs.Length; j++)
                m0[j, 0] = f[j];

            vs[i] = new(m0);
        }

        return KMatrix<K>.MergeSameRows(vs);
    }

    static KMatrix<FracPoly<EPoly<K>>> MatrixEndo<K>(EPoly<K> b) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var x = FG.KFracPoly('x', b);
        var mb = Endo(b).ToEMatrix(b.F);
        var P = mb - x * mb.One;
        return P;
    }

    public static (EPoly<K> Tr, EPoly<K> Norm) TraceNorm<K>(EPoly<K> b)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var mb = Endo(b).ToEMatrix(b.F);
        var norm = mb.Det;
        var tr = mb.Trace;
        return (tr, norm);
    }

    public static void CharacPoly<K>(EPoly<K> b) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var P = MatrixEndo(b);

        Console.WriteLine($"With {b.F} = 0");
        Console.WriteLine($"X0 = {b}");
        var f1 = P.Det;
        var X = FG.KPoly('X', b[0]);
        var f2 = f1.Num.Coefs.Select((c, i) => c[0] * X.Pow(i)).Aggregate(X.Zero, (acc, a) => acc + a);
        var sep = PolynomialFactorization.YunSFF(f2);
        var sep2 = sep.Select(e => (e.i, e.g.Substitute(X.Pow(e.q))))
            .Select(e => e.i == 1 ? $"({e.Item2})" : $"({e.Item2})^{e.i}").Glue(" * ");

        Console.WriteLine($"Characteristic Polynomial of X0 is f(X) = {sep2}");
        Console.WriteLine("f(X0) = {0}", f1.Substitute(b.ToKPoly('x')));
        Console.WriteLine();
    }

    public static KPoly<K> Norm<K>(KPoly<EPoly<K>> A, char c = 'x') where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var n = A.Degree;
        var g0 = A[0].F;
        var x = FG.KFracPoly(c, A.KZero.Poly.KZero);
        var y = FG.KPoly('y', x.Zero);
        var ga = y.Zero;
        for (int i = 0; i <= n; i++)
        {
            var gi = A[i].Poly.Coefs;
            var gy = new KPoly<FracPoly<K>>(y.x, y.KZero, gi.Select(k => k * y.KOne).ToArray());
            ga = ga + gy * x.Pow(i);
        }

        var ug = new KPoly<FracPoly<K>>(y.x, y.KZero, g0.Coefs.Select(k => k * y.KOne).ToArray());
        var res = Ring.DeterminantByPivot(Ring.SylvesterMatrix(ug, ga));
        return res.Num;
    }

    public static void NormDetails<K>(KPoly<EPoly<K>> A, char c = 'X') where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var norm = Norm(A, c);

        Console.WriteLine($"With {A[0].F} = 0");
        Console.WriteLine($"P = {A}");
        var sep = PolynomialFactorization.YunSFF(norm);
        var sep2 = sep.Select(e => (e.i, e.g.Substitute(norm.X.Pow(e.q))))
            .Select(e => e.i == 1 ? $"({e.Item2})" : $"({e.Item2})^{e.i}").Glue(" * ");

        var pow = sep.Select(e => norm.X.Pow(e.q)).Glue(fmt: "({0})");
        Console.WriteLine($"Norm(P) = f = {sep2} = [{norm}]{pow}");

        var f2 = norm.Coefs.Select((e, i) => (A.KOne * e) * A.X.Pow(i)).Aggregate((a, b) => a + b);
        var (q, r) = f2.Div(A);
        Console.WriteLine($"f/P = {q} rem {r}");
        Console.WriteLine();
    }

    // Barry Trager, Algebraic Factoring
    public static (int s, KPoly<EPoly<K>> g, KPoly<K> r) SqfrNorm<K>(KPoly<EPoly<K>> f, char c = 'X')
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var a = f[0].X;
        var x = f.X;
        var g = f.Substitute(x);
        for (int s = 0; s < 50; s++)
        {
            var r = Norm(g, c);
            if (Ring.Gcd(r, r.Derivative).Degree == 0)
                return (s, g, r);

            g = g.Substitute(x - a);
        }

        throw new Exception();
    }

    // Barry Trager, Algebraic Factoring
    public static List<(KPoly<EPoly<Rational>> hi, int i, int q)> AlgebraicFactors(KPoly<EPoly<Rational>> f, char c = 'X')
    {
        Console.WriteLine("############# START Algebraic Factors #############");
        var (s, g, r) = SqfrNorm(f);
        var sepFactors = PolynomialFactorization.YunSFF(r);
        var L = new List<(KPoly<EPoly<Rational>> hi, int i, int q)>();
        var x = g.X;
        var a = g[0].X;
        var g0 = g.Substitute(x);

        var k = 0;
        foreach (var (h0, q, i) in sepFactors)
        {
            var hs = PolynomialFactorizationPart2.FirrQ(h0);
            foreach (var h1 in hs)
            {
                var h2 = h1.Substitute(x);
                var h3 = Ring.Gcd(g0, h2);
                g0 = g0 / h3;
                var h4 = h3.Substitute(x + s * a);
                L.Add((h4.Monic, q, i));
                Console.WriteLine(new { k, hi = h4.Monic, q, i });
            }

            ++k;
        }

        Console.WriteLine();
        Console.WriteLine($"With {a.F} = 0");

        var prod = L.Aggregate(g.One, (acc, h) => acc * h.hi);
        var seq = L.Select(h => h.hi).Glue(" * ", "({0})");
        Console.WriteLine($"{prod} = {seq}");
        Console.WriteLine($"Is equal {f.Substitute(g.X).Equals(prod)}");
        Console.WriteLine("#############  END  Algebraic Factors #############");
        Console.WriteLine();

        return L;
    }

    // Barry Trager, Algebraic Factoring
    public static (KPoly<K> r, KPoly<K> a, KPoly<K> b) PrimitiveElt<K>(KPoly<EPoly<K>> f)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var (s, g, r) = SqfrNorm(f, 'y');
        var y = FG.EPoly(r, 'y');
        var x = FG.KPoly('x', y);

        var g0 = x.Zero;
        for (int i = 0; i <= g.Degree; i++)
        {
            var c = g[i].Poly;
            var g1 = c.Coefs.Select((k, i0) => (k * x.KOne) * x.Pow(i0)).Aggregate(x.Zero, (acc, xi) => acc + xi);
            g0 += g1 * y.Pow(i);
        }

        var p0 = f[0].F.Coefs.Select((k, i0) => (k * x.KOne) * x.Pow(i0)).Aggregate(x.Zero, (acc, xi) => acc + xi);
        var gcd = Ring.Gcd(g0, p0);

        var a = (-gcd[0]) / gcd[1];
        var b = y - s * a;
        return (r, a.Poly, b.Poly);
    }

    static KPoly<EPoly<Rational>> Resolve(List<KPoly<Rational>> minPolys, KPoly<EPoly<Rational>> f, List<EPoly<Rational>> roots)
    {
        if (f.Degree == 0)
            return f;

        var i = minPolys.Count;
        var (r0, a0, b0) = PrimitiveElt(f);
        var L0 = PolynomialFactorizationPart2.FirrQ(r0);
        var r1 = L0[0];
        var (x1, y0) = FG.EPolyXC(r1, (char)(f[0].Poly.x + 1));
        var xa0 = a0.Substitute(y0)[0];
        var xb0 = b0.Substitute(y0)[0];
        var f0 = new KPoly<EPoly<Rational>>(x1.x, x1.KZero,
            f.Coefs.Select(k => new EPoly<Rational>(xa0.F, k.Substitute(xa0).Poly)).ToArray());
        var f1 = f0 / (x1 - xb0);
        if (f1.Degree == 0)
        {
            var xb = (-f[0]) / f[1];
            roots.Add(xb);
            Console.WriteLine("MinPolys");
            Console.WriteLine(minPolys.Glue("\n"));
            Console.WriteLine("Roots");
            Console.WriteLine(roots.Glue("\n"));
            Console.WriteLine();

            return f1;
        }


        minPolys.Add(r1);
        var tmp = roots.Select(s => s.Substitute(xa0)).ToList();
        roots.Clear();
        roots.AddRange(tmp);
        roots.Add(xb0);
        Console.WriteLine($"f{i} = {f} => f{i + 1} = {f1}");

        return f1;
    }

    public static void SplittingField(KPoly<Rational> f)
    {
        var (x, a) = FG.EPolyXC(f, 'a');
        var f0 = f.Substitute(x);
        var minPolys = new List<KPoly<Rational>>() { f };
        var roots = new List<EPoly<Rational>>() { a[0] };
        Console.WriteLine($"f0 = {f}");
        var tmp = f0 / (x - a);

        while (tmp.Degree > 0)
        {
            tmp = Resolve(minPolys, tmp, roots);
        }

        var mp = minPolys.Last();
        var mp0 = roots.Last().F;
        Console.WriteLine(new { mp, mp0 });
        var X = FG.KPoly('X', roots.Last());
        var f1 = roots.Select(r => X - r).Aggregate(X.One, (acc, xi) => acc * xi);
        Console.WriteLine($"prod = {f1}");
        Console.WriteLine($"f = {f}");
        Console.WriteLine();
    }

    public static void MinimalPolynomials()
    {
        {
            var x = FG.QPoly('a');
            var f = x.Pow(2) - 5;
            var a = FG.EPoly(f);
            var X = FG.KPoly('X', a);
            CharacPoly(a);
            NormDetails(X - a);
            CharacPoly((a + 3) / 2);
            NormDetails(X - (a + 3) / 2);
        }

        {
            var x = FG.QPoly('a');
            var f = x.Pow(2) - 3 * x + 1;
            var a = FG.EPoly(f);
            var X = FG.KPoly('X', a);
            CharacPoly(a);
            NormDetails(X - a);
            CharacPoly(2 * a - 3);
            NormDetails(X - (2 * a - 3));
        }

        {
            var x = FG.QPoly('a');
            var f = x.Pow(4) + x.Pow(3) + x.Pow(2) + x + 1;
            var a = FG.EPoly(f);
            var X = FG.KPoly('X', a);
            CharacPoly(a);
            NormDetails(X - a);
            CharacPoly(a + a.Inv());
            NormDetails(X - (a + a.Inv()));
        }

        {
            var a = FG.EQPoly('a', -2, 0, 0, 1);
            var X = FG.KPoly('X', a);
            NormDetails(X.Pow(3) - (a - 1));
        }
    }

    public static void AlgFactorization()
    {
        Monom.Display = MonomDisplay.Caret;
        {
            var i = FG.EQPoly('i', 1, 0, 1);
            var x = FG.KPoly('x', i);
            var A = x.Pow(4) - x.Pow(2) - 2;
            AlgebraicFactors(A);
            // x^4 + -x^2 + -2 = (x + -i) * (x + i) * (x^2 + -2)
        }

        {
            var a = FG.EQPoly('a', -1, -3, 0, 1);
            var x = FG.KPoly('x', a);
            var f = x.Pow(3) - 3 * x - 1;
            Console.WriteLine($"f={f}; A = f/(x-a) = {f.Div(x - a)}");
            var A = f / (x - a);
            AlgebraicFactors(f);
            // x^3 + -3·x + -1 = (x + a^2 + -2) * (x + -a) * (x + -a^2 + a + 2)
        }

        {
            var a = FG.EQPoly('a', 1, 1, 1, 1, 1);
            var x = FG.KPoly('x', a);
            var A = x.Pow(2) - 5;
            AlgebraicFactors(A);
            // x^2 + -5 = (x + 2·a^3 + 2·a^2 + 1) * (x + -2·a^3 + -2·a^2 + -1)
        }

        {
            var a = FG.EQPoly('a', -2, 0, 0, 1);
            var x = FG.KPoly('x', a);
            var A = x.Pow(3) + 3 * x.Pow(2) + 3 * x - 1;
            AlgebraicFactors(A);
            // x^3 + 3·x^2 + 3·x + -1 = (x + -a + 1) * (x^2 + (a + 2)·x + a^2 + a + 1)
        }

        {
            var i = FG.EQPoly('i', 1, 0, 1);
            var x = FG.KPoly('x', i);
            var A = x.Pow(4) + 25 * x.Pow(2) + 50 * x + 25;
            AlgebraicFactors(A);
            // x^4 + 25·x^2 + 50·x + 25 = (x^2 + -5·i·x + -5·i) * (x^2 + 5·i·x + 5·i)
        }

        {
            var a = FG.EQPoly('a', 2, 2, 1);
            var x = FG.KPoly('x', a);
            var A = x.Pow(4) + 4;
            AlgebraicFactors(A);
            // x^4 + 4 = (x + -a) * (x + a) * (x + a + 2) * (x + -a + -2)
        }
    }
    
    public static void PrimitiveEltExamples()
    {
        {
            var (a, b, y) = Ring.Polynomial("a", "b", "y", Rational.KZero());
            y.Indeterminates.SetOrder(MonomOrder.Lex);
            var (p1, p2, p3) = (a.Pow(2) - 2, b.Pow(2) - 3, a + b - y);
            var gb = Ring.ReducedGrobnerBasis(p1, p2, p3);
            Console.WriteLine(gb.Glue("\n"));
            Console.WriteLine();
        }

        {
            var (x, _) = FG.EPolyXC(Rational.KZero(), 'a', -2, 0, 1);
            var f = x.Pow(2) - 3;
            var (r, a, b) = PrimitiveElt(f);
            Console.WriteLine(new[] { $"r = {r}", $"a = {a}", $"b = {b}" }.Glue("\n"));
            Console.WriteLine();
        }
        
        {
            var (a, b, y) = Ring.Polynomial("a", "b", "y", Rational.KZero());
            y.Indeterminates.SetOrder(MonomOrder.Lex);
            var (p1, p2, p3) = (a.Pow(3) + a + 1, b.Pow(3) - b.Pow(2) + 4 * b - 3, a + b - y);
            var gb = Ring.ReducedGrobnerBasis(p1, p2, p3);
            Console.WriteLine(gb.Glue("\n"));
            Console.WriteLine();
        }

        {
            var (x, _) = FG.EPolyXC(Rational.KZero(), 'a', 1, 1, 0, 1);
            var f = x.Pow(3) - x.Pow(2) + 4 * x - 3;
            var (r, a, b) = PrimitiveElt(f);
            Console.WriteLine(new[] { $"r = {r}", $"a = {a}", $"b = {b}" }.Glue("\n"));
            Console.WriteLine();
        }

        {
            var (a, b, y) = Ring.Polynomial("a", "b", "y", Rational.KZero());
            y.Indeterminates.SetOrder(MonomOrder.Lex);
            var (p1, p2, p3) = (a.Pow(4) + a.Pow(3) + a.Pow(2) + a + 1, b.Pow(2) - 5, a + b - y);
            var gb = Ring.ReducedGrobnerBasis(p1, p2, p3);
            Console.WriteLine(gb.Glue("\n"));
            Console.WriteLine();
        }

        {
            var (x, _) = FG.EPolyXC(Rational.KZero(), 'a', 1, 1, 1, 1, 1);
            var f = x.Pow(2) - 5;
            var (r, a, b) = PrimitiveElt(f);
            Console.WriteLine(new[] { $"r = {r}", $"a = {a}", $"b = {b}" }.Glue("\n"));
            Console.WriteLine();
        }
    }

    public static void SplittingFieldCubicPolynomial()

    {
        var x = FG.QPoly('y');
        
        SplittingField(x.Pow(2) - 3 * x - 3);
        SplittingField(x.Pow(2) + x + 1);
        SplittingField(x.Pow(2) + 2 * x - 5);
        
        SplittingField(x.Pow(3) - 2);
        SplittingField(x.Pow(3) - 3);
        SplittingField(x.Pow(3) - 3 * x - 1);
        SplittingField(x.Pow(3) - 2 * x + 2);
        SplittingField(x.Pow(3) + 2 * x.Pow(2) - x - 1);
        SplittingField(x.Pow(3) + 4 * x.Pow(2) + 3 * x + 1);
        SplittingField(x.Pow(3) - x + 1);
        
        SplittingField(x.Pow(4) - 2);
        SplittingField(x.Pow(4) + 4 * x.Pow(2) + 2);
        SplittingField(x.Pow(4) - 4 * x.Pow(2) + 2);
        SplittingField(x.Pow(4) + x.Pow(3) + x.Pow(2) + x + 1);
        
        // SplittingField(x.Pow(4) + 3*x.Pow(3) - x.Pow(2) + x + 1); // wont work
    }
}

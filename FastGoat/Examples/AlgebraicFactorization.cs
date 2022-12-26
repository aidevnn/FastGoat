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

    static void CharacPoly<K>(EPoly<K> b) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var x = FG.KFracPoly('x', b);
        var mb = Endo(b).ToEMatrix(b.F);
        var P = mb - x * mb.One;

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

    static KPoly<K> Norm<K>(KPoly<EPoly<K>> A, char c = 'x') where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
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

    static void NormDetails<K>(KPoly<EPoly<K>> A, char c = 'X') where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var norm = Norm(A, c);

        Console.WriteLine($"With {A[0].F} = 0");
        Console.WriteLine($"P = {A}");
        var sep = PolynomialFactorization.YunSFF(norm);
        var sep2 = sep.Select(e => (e.i, e.g.Substitute(norm.X.Pow(e.q))))
            .Select(e => e.i == 1 ? $"({e.Item2})" : $"({e.Item2})^{e.i}").Glue(" * ");

        Console.WriteLine($"Norm(P) = f = {sep2} = [{norm}]");

        var f2 = norm.Coefs.Select((e, i) => (A.KOne * e) * A.X.Pow(i)).Aggregate((a, b) => a + b);
        var (q, r) = f2.Div(A);
        Console.WriteLine($"f/P = {q} rem {r}");
        Console.WriteLine();
    }

    // Barry Trager, Algebraic Factoring
    static (int s, KPoly<EPoly<K>> g, KPoly<K> r) SqfrNorm<K>(KPoly<EPoly<K>> f, char c = 'X')
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var a = f[0].X;
        var x = f.X;
        var g = f.Substitute(x);
        Console.WriteLine("###### START SqfrNorm");
        Console.WriteLine($"A = {g}");
        for (int s = 0; s < 50; s++)
        {
            var r = Norm(g, c);
            Console.WriteLine(new { s, g, r });
            if (Ring.Gcd(r, r.Derivative).Degree == 0)
            {
                Console.WriteLine("######  END  SqfrNorm");
                Console.WriteLine();
                return (s, g, r);
            }

            g = g.Substitute(x - a);
        }

        throw new Exception();
    }

    // Barry Trager, Algebraic Factoring
    static List<(KPoly<EPoly<Rational>> hi, int i, int q)> AlgebraicFactors(KPoly<EPoly<Rational>> f, char c = 'X')
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
            var hs = PolynomialFactorizationPart2.FirrZ(h0);
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

        var prod = L.Aggregate(g.One, (acc, h) => acc * h.hi);
        Console.WriteLine(new { prod });
        Console.WriteLine("#############  END  Algebraic Factors #############");
        Console.WriteLine();

        return L;
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
        {
            var i = FG.EQPoly('i', 1, 0, 1);
            var x = FG.KPoly('x', i);
            var A = x.Pow(4) - x.Pow(2) - 2;
            AlgebraicFactors(A);
        }

        {
            var a = FG.EQPoly('a', -1, -3, 0, 1);
            var x = FG.KPoly('x', a);
            var f = x.Pow(3) - 3 * x - 1;
            Console.WriteLine($"f={f}; A = f/(x-a) = {f.Div(x - a)}");
            var A = f / (x - a);
            AlgebraicFactors(A);
        }

        {
            var a = FG.EQPoly('a', 1, 1, 1, 1, 1);
            var x = FG.KPoly('x', a);
            var A = x.Pow(2) - 5;
            AlgebraicFactors(A);
        }
        
        {
            var a = FG.EQPoly('a', -2, 0, 0, 1);
            var x = FG.KPoly('x', a);
            var A = x.Pow(3) + 3 * x.Pow(2) + 3 * x - 1;
            AlgebraicFactors(A);
        }
    }
}
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Polynoms;

namespace FastGoat.Examples;

public static class AlgebraicFactorization
{
    public static EPoly<K>[] Ebase<K>(EPoly<K> scalar) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var n = scalar.F.Degree;
        return n.Range().Select(i => scalar.X.Pow(i)).ToArray();
    }

    public static KMatrix<K> Endo<K>(EPoly<K> b) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
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
        var y = FG.KPoly(A.x, x.Zero);
        var ga = y.Zero;
        for (int i = 0; i <= n; i++)
        {
            var gi = A[i].Poly.Coefs;
            var gy = new KPoly<FracPoly<K>>(y.x, y.KZero, gi.Select(k => k * y.KOne).ToArray());
            ga = ga + gy * x.Pow(i);
        }

        var ug = new KPoly<FracPoly<K>>(y.x, y.KZero, g0.Coefs.Select(k => k * y.KOne).ToArray());
        var res = Ring.FastResultant(ug, ga);
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
        // Console.WriteLine($"SqfrNorm({f})");
        for (int s = 0; s < 150; s++)
        {
            // Console.WriteLine($"s={s} Norm({g})");
            var r = Norm(g, c);
            if (Ring.Gcd(r, r.Derivative).Degree == 0)
                return (s, g, r);

            g = g.Substitute(x - a);
        }

        throw new Exception();
    }

    // Barry Trager, Algebraic Factoring
    public static List<KPoly<EPoly<Rational>>> AlgebraicFactors(KPoly<EPoly<Rational>> f, bool details = false)
    {
        var (s, g, r) = SqfrNorm(f);
        var L = new List<KPoly<EPoly<Rational>>>();
        var x = g.X;
        var a = g[0].X;
        var g0 = g.Substitute(x);
        
        var hs = PolynomialFactorizationPart2.FirrZ(r);
        if (hs.Length == 1)
        {
            L.Add(f.Substitute(x));
            return L;
        }
        
        foreach (var h1 in hs)
        {
            var h2 = h1.Substitute(x);
            var h3 = Ring.StableGcd(g0, h2);
            g0 = g0 / h3;
            var h4 = h3.Substitute(x + s * a);
            L.Add(h4.Monic);
        }

        if (details)
        {
            Console.WriteLine($"With {a.F} = 0");
            var prod = L.Aggregate(g.One, (acc, h) => acc * h);
            var seq = L.Glue(" * ", "({0})");
            Console.WriteLine($"{prod} = {seq}");
            Console.WriteLine($"Is equal {f.Substitute(g.X).Equals(prod)}");
            Console.WriteLine();
        }

        return L;
    }

    public static EPoly<K> AlphaPrimElt<K>(KPoly<EPoly<K>> A, KPoly<EPoly<K>> newX)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        if (A.IsZero())
            return newX[0].Zero;

        var y = newX[0].X;
        var p0 = A[0].X.F.Substitute(newX);
        var g0 = newX.Zero;
        for (int i = 0; i <= A.Degree; i++)
        {
            var g1 = A[i].Poly.Substitute(newX);
            g0 += g1 * y.Pow(i);
        }

        var gcd = Ring.StableGcd(p0, g0);
        if (gcd.Degree == 0)
            throw new ArgumentException($"A={A} MinPoly={p0} A0={g0}");
        
        return (-gcd[0]) / gcd[1];
    }

    // Barry Trager, Algebraic Factoring
    public static (KPoly<K> r, KPoly<K> a, KPoly<K> b) PrimitiveElt<K>(KPoly<EPoly<K>> f)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var (s, g, r) = SqfrNorm(f, 'y');
        var y = FG.EPoly(r, 'y');
        var x = FG.KPoly('x', y);
        var a = AlphaPrimElt(g, x);
        var b = y - s * a;
        return (r, a.Poly, b.Poly);
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
            AlgebraicFactors(A, true);
            // x^4 + -x^2 + -2 = (x + -i) * (x + i) * (x^2 + -2)
        }

        {
            var a = FG.EQPoly('a', -1, -3, 0, 1);
            var x = FG.KPoly('x', a);
            var f = x.Pow(3) - 3 * x - 1;
            Console.WriteLine($"f={f}; A = f/(x-a) = {f.Div(x - a)}");
            var A = f / (x - a);
            AlgebraicFactors(f, true);
            // x^3 + -3·x + -1 = (x + a^2 + -2) * (x + -a) * (x + -a^2 + a + 2)
        }

        {
            var a = FG.EQPoly('a', 1, 1, 1, 1, 1);
            var x = FG.KPoly('x', a);
            var A = x.Pow(2) - 5;
            AlgebraicFactors(A, true);
            // x^2 + -5 = (x + 2·a^3 + 2·a^2 + 1) * (x + -2·a^3 + -2·a^2 + -1)
        }

        {
            var a = FG.EQPoly('a', -5, 0, 1);
            var x = FG.KPoly('x', a);
            var A = x.Pow(4) + x.Pow(3) + x.Pow(2) + x + 1;
            AlgebraicFactors(A, true);
            // x^4 + x^3 + x^2 + x + 1 = (x^2 + (1/2·a + 1/2)·x + 1) * (x^2 + (-1/2·a + 1/2)·x + 1)
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
            AlgebraicFactors(A, true);
            // x^4 + 25·x^2 + 50·x + 25 = (x^2 + -5·i·x + -5·i) * (x^2 + 5·i·x + 5·i)
        }

        {
            var a = FG.EQPoly('a', 2, 2, 1);
            var x = FG.KPoly('x', a);
            var A = x.Pow(4) + 4;
            AlgebraicFactors(A, true);
            // x^4 + 4 = (x + -a) * (x + a) * (x + a + 2) * (x + -a + -2)
        }

        {
            var a = FG.EQPoly('a', 3, 0, 1);
            var x = FG.KPoly('x', a);
            var A = x.Pow(6) - 1;
            AlgebraicFactors(A, true);
            // x^6 + -1 = (x + 1/2·a + 1/2) * (x + 1/2·a + -1/2) * (x + -1/2·a + 1/2) * (x + -1/2·a + -1/2) * (x + 1) * (x + -1)
        }
    }

    public static (KPoly<K> minPoly, KPoly<K> a, KPoly<K> b) PrimEltGb<K>(KPoly<EPoly<K>> f)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var d = ((char)(f[0].F.x + 1)).ToString();
        var (X, Y, T) = Ring.Polynomial("X", "Y", d, f.KZero.KZero, MonomOrder.Lex);
        var (x, y, t) = Y.Indeterminates;
        var P = f.Coefs.Select((k, i) => k.Poly.ToPolynomial(Y) * X.Pow(i)).Aggregate(X.Zero, (acc, xi) => acc + xi);
        var M = f[0].F.ToPolynomial(Y);
        for (int s = 1; s < 50; s++)
        {
            var PM = X + s * Y - T;
            var gb = Ring.ReducedGrobnerBasis(P, M, PM);
            if (gb.Length != 3)
                continue;

            var (ft, fy, fx) = (gb[0], gb[1], gb[2]);
            if (fx.DegreeOf(y) == 0 && fx.DegreeOf(x) == 1)
            {
                var r = FG.KPoly('y', f.KZero.KZero);
                var minPoly = ft.ToKPoly(t).Substitute(r);
                var rx = (X - fx).ToKPoly(t).Substitute(r);
                var ry = (Y - fy).ToKPoly(t).Substitute(r);
                return (minPoly, ry, rx);
            }
        }

        throw new Exception();
    }

    // Barry Trager, Algebraic Factoring
    public static List<EPoly<Rational>> SplittingField(KPoly<Rational> P,bool details = false)
    {
        var (X, y) = FG.EPolyXc(P, 'a');
        var P0 = P.Substitute(X);
        var roots = new List<EPoly<Rational>>();
        var polys = new List<KPoly<EPoly<Rational>>>() { P0 };
        var minPoly = P;
        var b = y;
        EPoly<Rational> new_y;
        var idx = 0;

        if (details)
            Console.WriteLine($"Polynomial P = {P0} in Q(a)[X]");

        while (true)
        {
            var BPoly = X.Zero;
            polys[idx] = polys[idx] / (X - b);
            roots.Add(b);
            var (quo,rem) = polys[idx].Div(X + b);
            if (rem.IsZero())
            {
                polys[idx] = quo;
                roots.Add(-b);
                if (quo.Degree == 0)
                    polys.Clear();
            }
            var (k, new_s) = (0, 0);
            var newFactors = new List<KPoly<EPoly<Rational>>>();

            foreach (var pi in polys)
            {
                var (s, g, R) = SqfrNorm(pi);
                var L = R.Degree < 13 ? PolynomialFactorizationPart2.FirrQ(R) : new[] { R };
                foreach (var qj in L.Order())
                {
                    var f = Ring.StableGcd(g, qj.Substitute(X));
                    if (qj.Degree > minPoly.Degree)
                    {
                        minPoly = qj;
                        idx = k;
                        new_s = s;
                        BPoly = f;
                    }

                    g /= f;
                    f = f.Substitute(X + s * y);
                    if (f.Degree == 1)
                    {
                        var r0 = -f[0] / f[1];
                        roots.Add(r0);
                    }
                    else
                    {
                        newFactors.Add(f);
                        ++k;
                    }
                }
            }

            (X, new_y) = FG.EPolyXc(minPoly, 'a');
            var a = AlphaPrimElt(BPoly, X);
            b = new_y - new_s * a;

            if (newFactors.Count == 0)
            {
                if (details)
                {
                    Console.WriteLine($"With {new_y.F} = 0");
                    roots.Println("Roots");
                    Console.WriteLine("Verif [Prod(X - ri)] = {0}", roots.Select(r => X - r).Aggregate((Xi, Xj) => Xi * Xj));
                    Console.WriteLine();
                }

                return roots;
            }

            roots = roots.Select(r => r.Substitute(a)).ToList();
            polys = newFactors.Select(f => f.SubstituteP0b(X, a)).ToList();
            y = new_y;
        }
    }

    // Barry Trager, Algebraic Factoring
    public static List<EPoly<ZnInt>> SplittingFieldFp(KPoly<ZnInt> P,bool details = false)
    {
        var prime = Un.FirstGen(P.KOne.P);
        var (X, y) = FG.EPolyXc(P, 'a');
        var P0 = P.Substitute(X);
        var roots = new List<EPoly<ZnInt>>();
        var polys = new List<KPoly<EPoly<ZnInt>>>() { P0 };
        var minPoly = P;
        var b = y;
        EPoly<ZnInt> new_y;
        var idx = 0;

        if (details)
            Console.WriteLine($"Polynomial P = {P0} in Z(a)[X] mod {prime.P}");

        while (true)
        {
            var BPoly = X.Zero;
            polys[idx] = polys[idx] / (X - b);
            roots.Add(b);
            var (k, new_s) = (0, 0);
            var newFactors = new List<KPoly<EPoly<ZnInt>>>();

            foreach (var pi in polys)
            {
                var (s, g, R) = SqfrNorm(pi);
                var L = PolynomialFactorization.Firr(R, prime);
                foreach (var qj in L.Order())
                {
                    var f = Ring.StableGcd(g, qj.Substitute(X));
                    if (qj.Degree > minPoly.Degree)
                    {
                        minPoly = qj;
                        idx = k;
                        new_s = s;
                        BPoly = f;
                    }

                    g /= f;
                    f = f.Substitute(X + s * y);
                    if (f.Degree == 1)
                    {
                        var r0 = -f[0] / f[1];
                        roots.Add(r0);
                    }
                    else
                    {
                        newFactors.Add(f);
                        ++k;
                    }
                }
            }

            (X, new_y) = FG.EPolyXc(minPoly, 'a');
            var a = AlphaPrimElt(BPoly, X);
            b = new_y - new_s * a;

            if (newFactors.Count == 0)
            {
                if (details)
                {
                    Console.WriteLine($"With {new_y.F} = 0");
                    roots.Println("Roots");
                    Console.WriteLine("Verif [Prod(X - ri)] = {0}", roots.Select(r => X - r).Aggregate((Xi, Xj) => Xi * Xj));
                    Console.WriteLine();
                }

                return roots;
            }

            roots = roots.Select(r => r.Substitute(a)).ToList();
            polys = newFactors.Select(f => f.SubstituteP0b(X, a)).ToList();
            y = new_y;
        }
    }

    public static List<EPoly<ZnInt>> SplittingFieldFp2(KPoly<ZnInt> P,bool details = false)
    {
        var (X, a) = FG.EPolyXc(P, 'a');
        var gf = new GFp("Gf", a);
        var g = Group.Generate(gf, a);
        
        var P0 = P.Substitute(X);
        var roots = new List<EPoly<ZnInt>>();
        foreach (var e in g)
        {
            if (P0.Div(X - e).rem.IsZero())
                roots.Add(e);
        }
    
        if (details)
        {
            Console.WriteLine($"Polynomial P = {P0} in Z(a)[X] mod {gf.P}");
            Console.WriteLine("Splitting field multiplicative group");
            DisplayGroup.Head(g);
            roots.Println("Roots");
            Console.WriteLine("Verif [Prod(X - ri)] = {0}", roots.Select(r => X - r).Aggregate((Xi, Xj) => Xi * Xj));
            Console.WriteLine();
        }

        return roots;
    }
    
    // Barry Trager, Algebraic Factoring
    public static void Resolvent<K>(KPoly<K> P) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var (X, y) = FG.EPolyXc(P, 'a');
        var Poly = P.Substitute(X);
        var minPoly = P;
        var b = y;
        var roots = new List<EPoly<K>>();
        while (true)
        {
            Poly = Poly / (X - b);
            roots.Add(b);
            var (quo,rem) = Poly.Div(X + b);
            if (rem.IsZero())
            {
                Poly = quo;
                roots.Add(-b);
            }
            if (Poly.Degree <= 1)
            {
                if (Poly.Degree == 1)
                    roots.Add(-Poly[0] / Poly[1]);
                
                Console.WriteLine($"A Resolvant for [P({P.x}) = {P}] is [R({minPoly.x}) = {minPoly}]");
                roots.Println("Roots");
                Console.WriteLine("Prod = {0}", roots.Select(r => X - r).Aggregate((ri, rj) => ri * rj));
                
                return;
            }

            var (s, g, R) = SqfrNorm(Poly);
            (X, y) = FG.EPolyXc(R, 'a');
            var a = AlphaPrimElt(g, X);
            b = y - s * a;
            Poly = Poly.SubstituteP0b(X, a);
            minPoly = R;
            roots = roots.Select(r => r.Substitute(a)).ToList();
        }
    }

    public static void PrimitiveEltExamples()
    {
        {
            var (b, _) = FG.EPolyXC(Rational.KZero(), 'a', 'b', -2, 0, 1);
            var f = b.Pow(2) - 3;
            var (r0, a0, b0) = PrimitiveElt(f);
            Console.WriteLine($"With {b[0].F} = 0 and {f} = 0");
            Console.WriteLine($"Trager Primitive element y with Q(y) = Q(a,b)");
            Console.WriteLine(new[] { "SqfrNorm", $"{r0} = 0", $"a = {a0}", $"b = {b0}" }.Glue("\n"));
            var (r1, a1, b1) = PrimEltGb(f);
            Console.WriteLine(new[] { "Grobner Basis", $"{r1} = 0", $"a = {a1}", $"b = {b1}" }.Glue("\n"));
            Console.WriteLine();
        }

        {
            var (b, _) = FG.EPolyXC(Rational.KZero(), 'a', 'b', 1, 1, 0, 1);
            var f = b.Pow(3) - b.Pow(2) + 4 * b - 3;
            var (r0, a0, b0) = PrimitiveElt(f);
            Console.WriteLine($"With {b[0].F} = 0 and {f} = 0");
            Console.WriteLine($"Primitive element y with Q(y) = Q(a,b)");
            Console.WriteLine(new[] { "Trager SqfrNorm", $"{r0} = 0", $"a = {a0}", $"b = {b0}" }.Glue("\n"));
            var (r1, a1, b1) = PrimEltGb(f);
            Console.WriteLine(new[] { "Grobner Basis", $"{r1} = 0", $"a = {a1}", $"b = {b1}" }.Glue("\n"));
            Console.WriteLine();
        }

        {
            var (b, _) = FG.EPolyXC(Rational.KZero(), 'a', 'b', 1, 1, 1, 1, 1);
            var f = b.Pow(2) - 5;
            var (r0, a0, b0) = PrimitiveElt(f);
            Console.WriteLine($"With {b[0].F} = 0 and {f} = 0");
            Console.WriteLine($"Primitive element y with Q(y) = Q(a,b)");
            Console.WriteLine(new[] { "Trager SqfrNorm", $"{r0} = 0", $"a = {a0}", $"b = {b0}" }.Glue("\n"));
            var (r1, a1, b1) = PrimEltGb(f);
            Console.WriteLine(new[] { "Grobner Basis", $"{r1} = 0", $"a = {a1}", $"b = {b1}" }.Glue("\n"));
            Console.WriteLine();
        }
    }

    public static void SplittingFieldQuarticPolynomial()
    {
        var x = FG.QPoly();

        SplittingField(x.Pow(2) - 3, true);
        SplittingField(x.Pow(2) - 3 * x - 3, true);
        SplittingField(x.Pow(2) + x + 1, true);
        SplittingField(x.Pow(2) + 2 * x - 5, true);
        
        SplittingField(x.Pow(3) - 2, true);
        SplittingField(x.Pow(3) - 3, true);
        SplittingField(x.Pow(3) - 3 * x - 1, true);
        SplittingField(x.Pow(3) - 2 * x + 2, true);
        SplittingField(x.Pow(3) + 2 * x.Pow(2) - x - 1, true);
        SplittingField(x.Pow(3) + 4 * x.Pow(2) + 3 * x + 1, true);
        SplittingField(x.Pow(3) - x + 1, true);
        
        SplittingField(x.Pow(4) + 4 * x.Pow(2) + 2, true);
        SplittingField(x.Pow(4) - 4 * x.Pow(2) + 2, true);
        SplittingField(x.Pow(4) + x.Pow(3) + x.Pow(2) + x + 1, true);
        SplittingField(x.Pow(4) - 2, true);
        SplittingField(x.Pow(4) + 2, true);
        SplittingField(x.Pow(4) + 5, true);
        SplittingField(x.Pow(4) + 3 * x.Pow(2) + 3, true);

        // SplittingField(x.Pow(4) + x + 1, true); // Time:70s
        // SplittingField(x.Pow(4) + 3*x.Pow(3) - x.Pow(2) + x + 1, true); // Time:113s
    }
}
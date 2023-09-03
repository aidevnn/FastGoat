using System.Numerics;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Padic;

namespace FastGoat.UserGroup.Polynoms;

public static partial class IntFactorisation
{
    public static EPoly<K>[] GetBase<K>(EPoly<K> a) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var bs = Array.Empty<EPoly<K>>();
        var n = a.F.Degree;
        for (int i = 0; i <= n - 1; ++i)
        {
            var ai = a.Pow(i);
            var bs0 = bs.Append(ai).ToArray();
            var mat0 = KMatrix<K>.MergeSameRows(bs0.Select(e => e.Poly.ToVMatrix(n - 1)).ToArray());
            if (mat0.NullSpace().nullity == 0)
                bs = bs0.ToArray();
        }

        return bs;
    }

    public static KPoly<K> Rewrite<K>(EPoly<K>[] bs, EPoly<K> b, char x = 'a') where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        if (bs.Any(a => !a.F.Equals(b.F)))
            throw new("Elements must belong to the same field");

        var n = b.F.Degree;
        var mat = KMatrix<K>.MergeSameRows(bs.Append(b).Select(e => e.Poly.ToVMatrix(n)).ToArray());
        var ns = mat.NullSpace();
        if (ns.nullity != 0)
        {
            return Ring.ReducedRowsEchelonForm(mat).A0.Cols.Last().ToKPoly(x);
        }

        return b.Zero.Poly.SubstituteChar('a');
    }

    public static KPoly<K> Rewrite<K>(EPoly<K> a, EPoly<K> b, char x = 'a') where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        return Rewrite(GetBase(a), b);
    }

    public static (EPoly<K>[], KPoly<K>) GetBaseAndMinPolynomial<K>(EPoly<K> a, char x = 'a')
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var bs = GetBase(a);
        var n = bs.Length;
        var P = Rewrite(bs, a.Pow(n), x);
        return (bs, P.X.Pow(n) - P);
    }

    public static (EPoly<Rational>[], KPoly<Rational>) ExtDegree(EPoly<Rational> a)
    {
        var (bs, minPol) = GetBaseAndMinPolynomial(a);
        var c = a.Poly.x == 'a' ? 'y' : 'a';
        Console.WriteLine($"[Q({c})/Q] = {bs.Length} with {c}={a} and {minPol} = 0");
        return (bs, minPol);
    }

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

    public static KPoly<K> CharacPoly<K>(EPoly<K> b, bool details = true) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var P = MatrixEndo(b);
        Console.WriteLine($"With {b.F} = 0");
        Console.WriteLine($"X0 = {b}");
        var f1 = P.Det;

        var X = FG.KPoly('X', b[0]);
        var f2 = f1.Num.Coefs.Select((c, i) => c[0] * X.Pow(i)).Aggregate(X.Zero, (acc, a) => acc + a);
        var sep = YunSFF(f2);

        if (details)
        {
            var sep2 = sep.Select(e => (e.i, e.g.Substitute(X.Pow(e.q))))
                .Select(e => e.i == 1 ? $"({e.Item2})" : $"({e.Item2})^{e.i}").Glue(" * ");

            Console.WriteLine($"Characteristic Polynomial of X0 is f(X) = {sep2}");
            Console.WriteLine("f(X0) = {0}", f1.Substitute(b.ToKPoly('x')));
            Console.WriteLine();
        }

        var (g, q, _) = sep[0];
        return g.Substitute(X.Pow(q));
    }

    public static KPoly<K> CharacPoly2<K>(EPoly<K> b, char c = 'x') where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var C = FG.KPoly(c, b);
        var norm = Norm(C - b, c);
        var sep = YunSFF(norm);
        return sep.Select(e => (e.i, e.g.Substitute(norm.X.Pow(e.q)))).First().Item2.SubstituteChar(c);
    }

    public static KPoly<K> MinPolynomial<K>(KPoly<EPoly<K>> A, char c = 'x') where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var norm = Norm(A, c);
        var sep = YunSFF(norm);
        if (sep.Count > 1)
            throw new();

        return sep.Select(e => (e.i, e.g.Substitute(norm.X.Pow(e.q)))).First().Item2.SubstituteChar(c);
    }

    public static KPoly<K> Norm<K>(KPoly<EPoly<K>> A, char c = 'x') where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var n = A.Degree;
        var g0 = A[0].F;
        var x = FG.KPoly(c, A.KZero.Poly.KZero);
        var y = FG.KPoly(A.x, x.Zero);
        var ga = y.Zero;
        for (int i = 0; i <= n; i++)
        {
            var xi = x.Pow(i);
            var gi = A[i].Poly.Coefs;
            var gy = new KPoly<KPoly<K>>(y.x, y.KZero, gi.Select(k => k * xi).ToArray());
            ga = ga + gy;
        }

        var yKOne = y.KOne;
        var ug = new KPoly<KPoly<K>>(y.x, y.KZero, g0.Coefs.Select(k => k * yKOne).ToArray());
        var res = Ring.FastResultant(ug, ga);
        return res;
    }

    public static KPoly<Rational> NormRationals(KPoly<EPoly<Rational>> A, char c = 'a')
    {
        var n = A.Degree;
        var g0 = A[0].F;
        var x = FG.KPoly(c, A.KZero.Poly.KZero);
        var y = FG.KPoly(A.x, x.Zero);
        var ga = y.Zero;
        for (int i = 0; i <= n; i++)
        {
            var xi = x.Pow(i);
            var gi = A[i].Poly.Coefs;
            var gy = new KPoly<KPoly<Rational>>(y.x, y.KZero, gi.Select(k => k * xi).ToArray());
            ga = ga + gy;
        }

        var yKOne = y.KOne;
        var ug = new KPoly<KPoly<Rational>>(y.x, y.KZero, g0.Coefs.Select(k => k * yKOne).ToArray());
        var res = Ring.FastResultant(ug, ga).SubstituteChar(y.x);
        return res;
    }

    public static void NormDetails<K>(KPoly<EPoly<K>> A, char c = 'X') where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var norm = Norm(A, c);

        Console.WriteLine($"With {A[0].F} = 0");
        Console.WriteLine($"P = {A}");
        var sep = YunSFF(norm);
        var sep2 = sep.Select(e => (e.i, e.g.Substitute(norm.X.Pow(e.q))))
            .Select(e => e.i == 1 ? $"({e.Item2})" : $"({e.Item2})^{e.i}").Glue(" * ");

        var pow = sep.Select(e => norm.X.Pow(e.q)).Glue(fmt: "({0})");
        Console.WriteLine($"Norm(P) = f = {sep2} = [{norm}]{pow}");
        //
        // var f2 = norm.Coefs.Select((e, i) => (A.KOne * e) * A.X.Pow(i)).Aggregate((a, b) => a + b);
        // var (q, r) = f2.Div(A);
        // Console.WriteLine($"f/P = {q} rem {r}");
        Console.WriteLine();
    }

    // Barry Trager, Algebraic Factoring
    public static (int s, KPoly<EPoly<K>> g, KPoly<K> r) SqfrNorm<K>(KPoly<EPoly<K>> f)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var a = f[0].X;
        var x = f.X;
        // Console.WriteLine($"SqfrNorm({f})");
        var coefs = 150.Range(-74).OrderBy(i => Int32.Abs(i)).ThenAscending().ToArray();
        foreach (var s in coefs.Where(e => e != 0))
        {
            try
            {
                var g = f.Substitute(x - a * s);
                // Console.WriteLine($"s={s} Norm({g})");
                // Console.WriteLine($"s={s})");
                // Console.WriteLine($"Norm({g}");
                var r = Norm(g);
                // if (r.Coefs.Any(c0 => c0 is Rational c1 && !c1.IsInteger()))
                //     continue;

                // Console.WriteLine($" = {r}");
                if (Ring.FastGCD(r, r.Derivative).Degree == 0) // dilemma between Ring.Gcd and Ring.FastGCD
                {
                    return (s, g.Monic, r);
                }
            }
            catch (Exception)
            {
                // Console.WriteLine(e);
            }
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

        var (nf, c0) = ConstCoef(r, monic: true);
        var hs = FirrZ2(nf, details).Select(f0 => f0.Substitute(f0.X / c0)).ToArray();
        if (hs.Length == 1)
        {
            L.Add(f.Substitute(x));
            if (details)
                Console.WriteLine($"f = {L[0]} is irreductible");
        }
        else
        {
            foreach (var h1 in hs)
            {
                var h2 = h1.Substitute(x);
                var h3 = Ring.FastGCD(g0, h2);
                g0 = (g0 / h3).Monic;
                var h4 = h3.Substitute(x + s * a);
                L.Add(h4.Monic);
            }
        }

        if (details)
        {
            Console.WriteLine($"f = {f} with f({a.F.x}) = 0");
            Console.WriteLine($"Square free norm : Norm(f({x.x} - {s}*{a.F.x}) = {r}");
            Console.WriteLine($"         = {hs.Glue(" * ", "({0})")}");
            Console.WriteLine();
            var prod = L.Aggregate(g.One, (acc, h) => acc * h);
            var seq = L.Glue(" * ", "({0})");
            Console.WriteLine($"{prod} = {seq}");
            Console.WriteLine($"Are equals {f.Substitute(g.X).Equals(prod)}");
            Console.WriteLine();
        }

        return L;
    }

    public static List<KPoly<EPoly<Rational>>> AlgebraicFactors(KPoly<Rational> f, bool details = false)
    {
        var (X, _) = FG.EPolyXc(f.Monic, 'y');
        return AlgebraicFactors(f.Substitute(X), details);
    }

    public static List<EPoly<Rational>> AlgebraicRoots(KPoly<EPoly<Rational>> f, bool details = false)
    {
        return AlgebraicFactors(f, details).Where(p => p.Degree == 1).Select(p => -p[0] / p[1]).Order().ToList();
    }

    public static List<EPoly<Rational>> AlgebraicRoots(KPoly<Rational> f, bool details = false)
    {
        var (X, _) = FG.EPolyXc(f.Monic, 'y');
        return AlgebraicRoots(f.Substitute(X), details);
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

        var gcd = Ring.FastGCD(p0, g0);
        if (gcd.Degree != 1)
            throw new ArgumentException($"A={A} MinPoly={p0} A0={g0} gcd={gcd}");

        return (-gcd[0]) / gcd[1];
    }

    // Barry Trager, Algebraic Factoring
    public static (KPoly<K> r, KPoly<K> a, KPoly<K> b) PrimitiveElt<K>(KPoly<EPoly<K>> f)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var (s, g, r) = SqfrNorm(f);
        var y = FG.EPoly(r, 'y');
        var x = FG.KPoly('x', y);
        var a = AlphaPrimElt(g, x);
        var b = y - s * a;
        return (r, a.Poly, b.Poly);
    }

    public static (KPoly<Rational> F, KPoly<EPoly<Rational>> X, EPoly<Rational> a, EPoly<Rational> b)[]
        PrimitiveElt(KPoly<Rational> f, KPoly<Rational> g)
    {
        var (X0, _) = FG.EPolyXc(f, 'y');
        var (r, a0, b0) = PrimitiveElt(g.Substitute(X0));
        return FirrZ2(r).Select(r0 =>
        {
            var (X1, y1) = FG.EPolyXc(r0, 'y');
            return (y1.F, X1, a0.Substitute(y1), b0.Substitute(y1));
        }).ToArray();
    }

    public static (KPoly<K> minPoly, KPoly<K> a, KPoly<K> b) PrimEltGb<K>(KPoly<EPoly<K>> f)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var d = ((char)(f[0].F.x + 1)).ToString();
        var (X, Y, T) = Ring.Polynomial(f.KZero.KZero, MonomOrder.Lex, "X", "Y", d).Deconstruct();
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
    public static List<EPoly<Rational>> SplittingField(KPoly<Rational> P, bool details = false)
    {
        GlobalStopWatch.Restart();
        if (FirrZ2(P).Length > 1)
            throw new($"{P} isnt an irreductible polynomial");

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
            var (quo, rem) = polys[idx].Div(X + b);
            if (rem.IsZero())
            {
                polys[idx] = quo;
                roots.Add(-b);
                if (quo.Degree == 0)
                {
                    polys.RemoveAt(idx);
                    if (idx >= polys.Count)
                        idx--;
                }
            }

            if (polys.Count > 0 && polys[idx].Degree == 1)
            {
                var pl = polys[idx];
                roots.Add(-pl[0] / pl[1]);
                polys.RemoveAt(idx);
                if (idx >= polys.Count)
                    idx--;
            }

            var (k, new_s) = (0, 0);
            var newFactors = new List<KPoly<EPoly<Rational>>>();

            foreach (var pi in polys)
            {
                var (s, g, R) = SqfrNorm(pi);
                var L = FirrZ2(R.Monic, details);
                foreach (var qj in L.OrderBy(e => e.Degree).ThenBy(e => e.NormInf()))
                {
                    var f = Ring.FastGCD(g, qj.Substitute(X));
                    if (qj.Degree > minPoly.Degree)
                    {
                        minPoly = qj;
                        idx = k;
                        new_s = s;
                        BPoly = f;
                    }

                    g = (g / f).Monic;
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

            (minPoly, var c0) = ConstCoef(minPoly, monic: true);
            BPoly = BPoly.Substitute(BPoly.X * (c0 * BPoly.KOne));
            (X, new_y) = FG.EPolyXc(minPoly, 'a');
            var a = AlphaPrimElt(BPoly, X);
            b = new_y * c0 - new_s * a;

            if (newFactors.Count == 0)
            {
                if (details)
                {
                    Console.WriteLine($"With {new_y.F} = 0");
                    roots.Println("Roots");
                    Console.WriteLine("Verif [Prod(X - ri)] = {0}", roots.Select(r => X - r).Aggregate((Xi, Xj) => Xi * Xj));
                    GlobalStopWatch.Show();
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
    public static List<EPoly<ZnInt>> SplittingFieldFp(KPoly<ZnInt> P, bool details = false)
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
                var L = Firr(R, prime);
                foreach (var qj in L.Order())
                {
                    var f = Ring.FastGCD(g, qj.Substitute(X));
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

    public static List<EPoly<ZnInt>> SplittingFieldFp2(KPoly<ZnInt> P, bool details = false)
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
        GlobalStopWatch.Restart();
        GlobalStopWatch.AddLap();
        var (X, y) = FG.EPolyXc(P, 'a');
        var Poly = P.Substitute(X);
        var minPoly = P;
        var b = y;
        var roots = new List<EPoly<K>>();
        while (true)
        {
            Poly = Poly / (X - b);
            roots.Add(b);
            var (quo, rem) = Poly.Div(X + b);
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
                GlobalStopWatch.Show("End");

                return;
            }

            var (s, g, R) = SqfrNorm(Poly);
            (X, y) = FG.EPolyXc(R, 'a');
            var a = AlphaPrimElt(g, X);
            b = y - s * a;
            Poly = Poly.SubstituteP0b(X, a);
            minPoly = R;
            roots = roots.Select(r => r.Substitute(a)).ToList();
            GlobalStopWatch.Show("Loop");
        }
    }
}
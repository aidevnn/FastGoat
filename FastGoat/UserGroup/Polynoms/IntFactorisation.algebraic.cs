using System.Numerics;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Floats;
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

    public static KPoly<K> Rewrite<K>(EPoly<K>[] bs, EPoly<K> b, char x = 'a')
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
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

    public static KPoly<K> Rewrite<K>(EPoly<K> a, EPoly<K> b, char x = 'a')
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
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

    public static KPoly<K> CheckQuadraticMinPolynomial<K>(EPoly<K> a, char x = 'a')
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var n = a.F.Degree;
        var mat = KMatrix<K>.MergeSameRows(3.Range().Select(k => (a.Pow(k)).Poly.ToVMatrix(n)).ToArray());
        var ns = mat.NullSpace();
        if (ns.nullity != 0)
        {
            var P = Ring.ReducedRowsEchelonForm(mat).A0.Cols.Last().ToKPoly(x);
            return P.X.Pow(2) - P;
        }

        return a.Zero.Poly.SubstituteChar(x);
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

    public static KPoly<K> CharacPoly<K>(EPoly<K> b) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var P = MatrixEndo(b);
        Console.WriteLine($"With {b.F} = 0");
        Console.WriteLine($"X0 = {b}");
        var f1 = P.Det;

        var X = FG.KPoly('X', b[0]);
        var f2 = f1.Num.Coefs.Select((c, i) => c[0] * X.Pow(i)).Aggregate(X.Zero, (acc, a) => acc + a);
        var sep = YunSFF(f2);

        if (Logger.Level != LogLevel.Off)
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

    public static KPoly<K> MinPolynomial<K>(KPoly<EPoly<K>> A, char c = 'x')
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var norm = Norm(A, c);
        var sep = YunSFF(norm).Where(e => e.g.Degree != 0).ToList();
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

    // Barry Trager, Algebraic Factoring
    public static (int s, KPoly<EPoly<K>> g, KPoly<K> r) SqfrNorm<K>(KPoly<EPoly<K>> f, bool integerOnly = false)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var a = f[0].X;
        var x = f.X;
        // Console.WriteLine($"SqfrNorm({f})");
        var coefs = 150.Range(-74).OrderBy(i => Int32.Abs(i)).ThenAscending().ToArray();
        foreach (var s in coefs)
        {
            try
            {
                var g = f.Substitute(x - a * s);
                var r = Norm(g);
                if (integerOnly && r is KPoly<Rational> r0 && r0.Coefs.Any(c => !c.Denom.IsOne))
                    continue;

                if (Ring.FastGCD(r, r.Derivative).Degree == 0) // dilemma between Ring.Gcd and Ring.FastGCD
                    return (s, g.Monic, r);
            }
            catch (Exception)
            {
                // Console.WriteLine(e);
            }
        }

        throw new Exception();
    }

    // Barry Trager, Algebraic Factoring
    public static List<KPoly<EPoly<Rational>>> AlgebraicFactors(KPoly<EPoly<Rational>> f)
    {
        var (s, g, r) = SqfrNorm(f, integerOnly: f.Degree > 4);
        var L = new List<KPoly<EPoly<Rational>>>();
        var x = g.X;
        var a = g[0].X;
        var g0 = g.Substitute(x);

        var (nf, c0) = ConstCoef(r, monic: true);
        var hs = FirrZ2(nf).Select(f0 => f0.Substitute(f0.X / c0)).ToArray();
        if (hs.Count(e => e.Degree > 0) == 1)
        {
            L.Add(f.Substitute(x));
            if (Logger.Level != LogLevel.Off)
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

        if (Logger.Level != LogLevel.Off)
        {
            Console.WriteLine($"f = {f} with f({a.F.x}) = 0");
            Console.WriteLine($"Square free norm : Norm(f({x - s * a})) = {r}");
            Console.WriteLine($"         = {hs.Glue(" * ", "({0})")}");
            Console.WriteLine();
            var prod = L.Aggregate(g.One, (acc, h) => acc * h);
            var seq = L.Glue(" * ", "({0})");
            Console.WriteLine($"{prod} = {seq}");
            Console.WriteLine($"Are equals {f.Substitute(g.X).Monic.Equals(prod.Monic)}");
            Console.WriteLine();
        }

        return L;
    }

    public static List<KPoly<EPoly<Rational>>> AlgebraicFactors(KPoly<Rational> f)
    {
        var (X, _) = FG.EPolyXc(f.Monic, 'y');
        return AlgebraicFactors(f.Substitute(X));
    }

    public static List<EPoly<Rational>> AlgebraicRoots(KPoly<EPoly<Rational>> f)
    {
        return AlgebraicFactors(f).Where(p => p.Degree == 1).Select(p => -p[0] / p[1]).Order().ToList();
    }

    public static List<EPoly<Rational>> AlgebraicRoots(KPoly<Rational> f)
    {
        var (X, _) = FG.EPolyXc(f.Monic, 'y');
        return AlgebraicRoots(f.Substitute(X));
    }

    public static (Rational disc, KPoly<EPoly<Rational>>[] roots) FactorsQuadratic(KPoly<Rational> P, char a = 'a')
    {
        if (P.Degree != 2 || !P.LC.Equals(Rational.KOne()))
            throw new("P must be monic and quadratic");

        var (P1, c) = ConstCoef(P, monic: true);
        P1 = P1.Monic;
        var D = P1[1].Pow(2) - 4 * P1[0];
        var (numD, _) = new Rational(D.Sign * IntExt.LcmBigInt(D.Num * D.Sign, D.Denom)).Decomp();
        var numD0 = numD.Where(e => e.Item2 % 2 != 0).ToArray();
        if (numD0.Length == 0)
            throw new("P must be irreductible");

        var D0 = numD0.Select(e => new Rational(e.Item1)).Aggregate((a0, a1) => a0 * a1);
        var (X, y) = FG.EPolyXc(P.X.Pow(2) - D0, a, P.x);
        var P2 = P1.Substitute(X);
        var res = AlgebraicFactors(P2).Select(p => p.Substitute(p.X / (p.KOne * c)).Monic).ToArray();
        if (Logger.Level != LogLevel.Off)
        {
            res.Println($"    D = {D} and {y} = Sqrt({D0})");
            Console.WriteLine();
        }

        return (D0, res);
    }

    public static Polynomial<Rational, Xi> PrettyPrintCnf(Cnf cf)
    {
        var e = cf.Simplify().E;
        if (e.Degree == 0)
            return e.Poly.ToPolynomial(Ring.Polynomial(Rational.KZero()));

        var minPol = CheckQuadraticMinPolynomial(e, 'x');
        if (minPol.Degree != 2)
            return cf.E.Poly.ToPolynomial(Ring.Polynomial(Rational.KZero(), MonomOrder.RevLex,
                $"{Cnf.RootsOfUnit}{cf.N}")[0]);

        var (disc, roots) = FactorsQuadratic(minPol);
        var (x0, x1) = (-roots[0][0] / roots[0][1], -roots[1][0] / roots[1][1]);
        var cplx0 = new Cplx(Complex.FromPolarCoordinates(cf.Module, cf.Phase));
        var sqrtD = new Cplx(Complex.Sqrt((double)disc));

        var adisc = disc * disc.Sign;
        var ind = adisc.Equals(disc.One)
            ? (disc.Sign == 1 ? "a" : "I")
            : (disc.Sign == 1 ? $"\u221a{disc}" : $"I\u221a{-disc}");
        if (x0.Poly.Substitute(sqrtD).Equals(cplx0))
            return x0.Poly.ToPolynomial(Ring.Polynomial(Rational.KZero(), MonomOrder.RevLex, ind)[0]);
        else
            return x1.Poly.ToPolynomial(Ring.Polynomial(Rational.KZero(), MonomOrder.RevLex, ind)[0]);
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
        var (x, y, t) = Y.Indeterminates.Deconstruct();
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
    public static List<EPoly<Rational>> SplittingField(KPoly<Rational> P)
    {
        GlobalStopWatch.Restart();
        if (FirrZ2(P).Count(e => e.Degree > 0) > 1)
        {
            throw new($"{P} isnt an irreductible polynomial");
        }

        var (X, y) = FG.EPolyXc(P, 'a');
        var P0 = P.Substitute(X);
        var roots = new List<EPoly<Rational>>();
        var polys = new List<KPoly<EPoly<Rational>>>() { P0 };
        var minPoly = P;
        var b = y;
        EPoly<Rational> new_y;
        var idx = 0;

        if (Logger.Level != LogLevel.Off)
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
                var (s, g, R) = SqfrNorm(pi, integerOnly: true);
                var L = FirrZ2(R.Monic);
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
                if (Logger.Level != LogLevel.Off)
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
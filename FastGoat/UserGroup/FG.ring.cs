using System.Numerics;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Characters;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Matrix;
using FastGoat.UserGroup.Padic;
using FastGoat.UserGroup.Polynoms;

namespace FastGoat.UserGroup;

public static partial class FG
{
    public static Rational Comb(int k, int n)
    {
        var c = Rational.KOne();
        if (k > n || k < 0)
            return c;

        var num = c;
        var denom = c;
        k = 2 * k > n ? n - k : k;
        for (int i = 0; i < k; i++)
        {
            num *= (n - i);
            denom *= (i + 1);
        }

        return num / denom;
    }
    
    public static ((int, int)[]num, (int, int)[] denom) Decomp(this Rational r)
    {
        var num = IntExt.PrimesDec(r.Num).Select(e => (e.Key, e.Value)).ToArray();
        var denom = IntExt.PrimesDec(r.Denom).Select(e => (e.Key, e.Value)).ToArray();
        return (num, denom);
    }

    public static KPoly<ZnInt> ZPoly(int p, char x = 'x') => new KPoly<ZnInt>(x, ZnInt.ZnZero(p)).X;
    public static KPoly<ZnBInt> ZbPoly(int p, char x = 'x') => new KPoly<ZnBInt>(x, ZnBInt.ZnZero(p)).X;
    public static KPoly<Rational> QPoly(char x = 'x') => new(x);
    public static KPoly<Cplx> CplxPoly(char x = 'x') => new(x);
    public static KPoly<BigCplx> BCplxPoly(int o = 40, char x = 'x') => new KPoly<BigCplx>(x, BigCplx.BcZero(o)).X;

    public static Cplx NSolve(KPoly<Cplx> P, double epsilon = 1e-14, int maxLoop = 200)
    {
        var dP = P.Derivative;
        var ai = Cplx.CZero;
        var aj = new Cplx(Double.Pi + Complex.ImaginaryOne);
        var i = 0;
        while ((ai - aj).NormInf > epsilon && i++ < maxLoop)
        {
            ai = aj;
            aj = ai - (P.Substitute(ai) / dP.Substitute(ai)); // Newton iteration
        }

        return aj;
    }

    public static Cplx[] NRoots(KPoly<Cplx> P, double epsilon = 1e-14, int maxLoop = 200)
    {
        var P0 = P;
        var roots = new List<Cplx>();
        while (P0.Degree > 0)
        {
            var a0 = NSolve(P0, epsilon);
            roots.Add(a0);
            P0 /= P0.X - a0;
        }

        return roots.ToArray();
    }

    public static BigCplx NSolve(KPoly<BigCplx> P, int maxLoop = 200)
    {
        var o = P.KZero.O;
        var aj = new BigCplx(BigReal.Pi(o), BigReal.E(o));
        return NSolve(P, aj, maxLoop);
    }

    public static BigCplx NSolve(KPoly<BigCplx> P, BigCplx a0, int maxLoop = 200)
    {
        var dP = P.Derivative;
        var i = 0;
        var ai = a0 * 1000;
        var aj = a0;
        while (!(ai - aj).IsZero() && i++ < maxLoop)
        {
            ai = aj;
            aj = ai - (P.Substitute(ai) / dP.Substitute(ai)); // Newton iteration
        }

        return aj;
    }

    public static BigCplx[] NRoots(KPoly<BigCplx> P, int maxLoop = 200)
    {
        var P0 = P;
        var roots = new List<BigCplx>();
        while (P0.Degree > 0)
        {
            var a0 = NSolve(P0, maxLoop);
            var a1 = NSolve(P, a0, maxLoop);
            roots.Add(a1);
            P0 /= P0.X - a1;
        }

        return roots.ToArray();
    }

    public static KPoly<Cplx> ToCPoly(this KPoly<Rational> P) => new(P.x, Cplx.CZero, P.Coefs.Select(c => c * Cplx.COne).ToArray());

    public static KPoly<BigReal> ToBrPoly(this KPoly<Rational> P, int O = 40) =>
        new(P.x, BigReal.BrZero(O), P.Coefs.Select(c => BigReal.BrOne(O) * BigReal.FromRational(c, O)).ToArray());

    public static KPoly<BigCplx> ToBcPoly(this KPoly<Rational> P, int O = 40) =>
        new(P.x, BigCplx.BcZero(O), P.Coefs.Select(c => BigCplx.BcOne(O) * BigCplx.FromRational(c, O)).ToArray());

    public static KPoly<BigCplx> ToRoundBcPoly(this KPoly<BigCplx> P, int d = 40) =>
        new(P.x, BigCplx.Round(P.KZero, d), P.Coefs.Select(c => BigCplx.Round(c, d)).ToArray());

    public static KPoly<BigCplx> ToBcPoly(this KPoly<BigCplx> P, int O = 40) =>
        new(P.x, BigCplx.BcZero(O), P.Coefs.Select(c => c.ToBigCplx(O)).ToArray());

    public static KPoly<Rational> ToIntPoly(this KPoly<BigCplx> P, int err = 4)
    {
        var one = P.KOne;
        if (err < 4 || one.O < 2 * err)
            throw new($"Digits err {err} must be >=4 and <=O/2");

        Rational Cv(BigCplx r)
        {
            if (!r.ImaginaryPart.ToBigReal(one.O - err).IsZero())
                throw new("Only real numbers allowed");

            if (r.ToBigCplx(one.O - err).IsZero())
                return Rational.KZero();

            return r.RealPart.ToBigReal(one.O - err).ToRational.RoundEven;
        }

        return new(P.x, Rational.KZero(), P.Coefs.Select(Cv).ToArray());
    }

    public static KPoly<Rational> ToIntPoly(this KPoly<BigReal> P)
    {
        var one = P.KOne;

        Rational Cv(BigReal r)
        {
            if (r.IsZero())
                return Rational.KZero();

            return BigReal.Round(r, one.O - 4).ToRational;
        }

        return new(P.x, Rational.KZero(), P.Coefs.Select(Cv).ToArray());
    }

    public static KPoly<Cplx> ToCPoly(this KPoly<EPoly<Rational>> P, Cplx e) =>
        new(P.x, Cplx.CZero, P.Coefs.Select(c => c.Poly.Substitute(e)).ToArray());

    public static KPoly<Rational> ToAbsKPoly(this KPoly<Rational> P) =>
        new(P.x, P.KZero, P.Coefs.Select(c => Rational.Absolute(c)).ToArray());

    public static Rational NormB(this KPoly<Rational> P, int b)
    {
        if (b < 1)
            throw new();

        return P.Coefs.Aggregate(P.KZero, (acc, c) => acc + Rational.Absolute(c).Pow(b));
    }

    public static Rational NormInf(this KPoly<Rational> P)
    {
        return P.Coefs.Max(c => Rational.Absolute(c));
    }

    public static int NbCoeffs<K>(this KPoly<K> f) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        return f.Coefs.Count(c => !c.IsZero());
    }

    public static KPoly<K> KPoly<K>(char x, K scalar) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        return new KPoly<K>(x, scalar).X;
    }

    public static FracPoly<K> KFracPoly<K>(char x, K scalar) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        return new FracPoly<K>(x, scalar).X;
    }

    public static EPoly<K> EPoly<K>(KPoly<K> f) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        return new(f);
    }

    public static EPoly<K> EPoly<K>(KPoly<K> f, char x) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        return new(new KPoly<K>(x, f.KZero, f.Coefs));
    }

    public static EPoly<Rational> CyclotomicEPoly(int n) => FG.EPoly(FG.CyclotomicPolynomial(n), 'a');

    public static KPoly<K> KPoly<K>(char x, params K[] coefs)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var coefs0 = coefs.Reverse().SkipWhile(i => i.IsZero()).Reverse().ToArray();
        if (coefs0.Length < 2)
            throw new GroupException(GroupExceptionType.GroupDef);

        return new KPoly<K>(x, coefs[0].Zero, coefs0);
    }

    private static KPoly<K> KPoly<K>(char x, K scalar, dynamic[] coefs)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var coefs0 = new List<K>();
        var cs = coefs[0] is Array a ? a : coefs;
        foreach (var c in cs)
        {
            if (c is int c0)
                coefs0.Add(c0 * scalar.One);
            else if (c is K c1)
                coefs0.Add(c1);
            else
                throw new ArgumentException();
        }

        var coefs1 = coefs0.Reverse<K>().SkipWhile(i => i.IsZero()).Reverse().ToArray();
        if (coefs1.Length < 2)
            throw new GroupException(GroupExceptionType.GroupDef);

        return new KPoly<K>(x, scalar.Zero, coefs1);
    }

    public static KPoly<K> KPoly<K>(K scalar, char x, params dynamic[] coefs)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        return KPoly(x, scalar, coefs);
    }

    public static EPoly<K> EPoly<K>(K scalar, char x, params dynamic[] coefs)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        return new EPoly<K>(KPoly(x, scalar, coefs));
    }

    public static (KPoly<EPoly<K>> X, KPoly<EPoly<K>> c) EPolyXC<K>(KPoly<K> f, char c, char x = 'X')
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var c0 = EPoly(f, c);
        var x0 = KPoly(x, c0);
        return (x0, c0 * x0.One);
    }

    public static (KPoly<EPoly<K>> X, EPoly<K> c) EPolyXc<K>(KPoly<K> f, char c, char x = 'X')
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var c0 = EPoly(f, c);
        var x0 = KPoly(x, c0);
        return (x0, c0);
    }

    public static (KPoly<EPoly<K>> X, KPoly<EPoly<K>> c) EPolyXC<K>(KPoly<K> f)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        return EPolyXC(f, f.x);
    }

    public static (KPoly<EPoly<K>> X, KPoly<EPoly<K>> c) EPolyXC<K>(char x, params K[] coefs)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        return EPolyXC(KPoly(x, coefs), x);
    }

    public static (KPoly<EPoly<K>> X, KPoly<EPoly<K>> c) EPolyXC<K>(K scalar, char x, params int[] coefs)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        return EPolyXC(KPoly(scalar, x, coefs), x);
    }

    public static (KPoly<EPoly<K>> X, KPoly<EPoly<K>> c) EPolyXC<K>(K scalar, char a, char b, params int[] coefs)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        return EPolyXC(KPoly(scalar, a, coefs), a, b);
    }

    public static EPoly<K> EPoly<K>(K scalar, char x, params int[] coefs)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        return EPoly(KPoly(scalar, x, coefs));
    }

    public static EPoly<Rational> EQPoly(char x, params int[] coefs) => EPoly(Rational.KZero(), x, coefs);

    public static FracPoly<Rational> QFracPoly(char x = 'x') => KFracPoly(x, Rational.KZero());

    public static FracPoly<ZnInt> ZFracPoly(int p, char x = 'x') => KFracPoly(x, ZnInt.ZnZero(p));

    public static (KPoly<FracPoly<ZnInt>> x, FracPoly<ZnInt> t) FpT_Poly(int p, (char x, char t) xt)
    {
        var t = ZPoly(p, xt.t);
        var t0 = new FracPoly<ZnInt>(t);
        return (KPoly(xt.x, t0), t0);
    }

    public static (KPoly<FracPoly<Rational>> X, FracPoly<Rational> T) QT_Poly(char x, char t)
    {
        var t0 = QPoly(t);
        var t1 = new FracPoly<Rational>(t0);
        return (KPoly(x, t1), t1);
    }

    public static (KPoly<FracPoly<ZnInt>> x, FracPoly<ZnInt> t) FpT_Poly(int p) => FpT_Poly(p, ('x', 't'));

    public static (KPoly<EPoly<ZnInt>> x, EPoly<ZnInt> a) FqX_Poly(int q, (char x, char a) xa)
    {
        var a = FqX(q, xa.a);
        return (KPoly(xa.x, a), a);
    }

    public static (KPoly<EPoly<ZnInt>> x, EPoly<ZnInt> a) FqX_Poly(int q) => FqX_Poly(q, ('x', 'a'));

    public static GL GLnp(int n, int p)
    {
        return new GL(n, p);
    }

    public static GLn<K> GLnK<K>(int n, K scalar) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        return new GLn<K>(n, scalar);
    }

    public static GLn<K> GLnK<K>(string name, int n, K scalar) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        return new GLn<K>(name, n, scalar);
    }

    public static EPolynomial<Rational>[] NumberFieldQ(Polynomial<Rational, Xi>[] basis)
    {
        var e0 = basis[0];
        var polyBasis = new PolynomialBasis<Rational, Xi>(e0.Indeterminates, basis);
        return e0.Indeterminates.Select(xi => new Polynomial<Rational, Xi>(new Monom<Xi>(e0.Indeterminates, xi, 1), e0.KOne))
            .Select(xi => new EPolynomial<Rational>(xi, polyBasis)).ToArray();
    }

    public static EPolynomial<Rational> NumberFieldQ(Polynomial<Rational, Xi> e0) => NumberFieldQ(new[] { e0 })[0];

    public static (EPolynomial<Rational>, EPolynomial<Rational>) NumberFieldQ(Polynomial<Rational, Xi> e0, Polynomial<Rational, Xi> e1)
    {
        var nbf = NumberFieldQ(new[] { e0, e1 });
        var x0 = e0.ExtractIndeterminate;
        var i0 = e0.Indeterminates.ToList().FindIndex(xi => xi.Equals(x0));
        var x1 = e1.ExtractIndeterminate;
        var i1 = e0.Indeterminates.ToList().FindIndex(xi => xi.Equals(x1));
        return (nbf[i0], nbf[i1]);
    }

    public static EPolynomial<Rational> NumberFieldQ(KPoly<Rational> e, string x)
    {
        var (_, t) = Ring.Polynomial(x, "_t_", Rational.KZero());
        var a = e.ToPolynomial(t.Indeterminates, t.Indeterminates[0]);
        return NumberFieldQ(a);
    }

    public static (EPolynomial<Rational> x, EPolynomial<Rational> p0) NumberFieldQ(KPoly<Rational> e, string x, string p0)
    {
        var all = new[] { x, p0, "_t_" };
        var xis = Ring.Polynomial(Rational.KZero(), MonomOrder.Lex, all);
        var x0 = xis[0];
        var a = e.ToPolynomial(x0.Indeterminates, x0.Indeterminates[0]);
        var nf = NumberFieldQ(a);
        return (nf, new EPolynomial<Rational>(xis[1], nf.Basis));
    }

    public static (EPolynomial<Rational>, EPolynomial<Rational>) NumberFieldQ((KPoly<Rational>, string) e0,
        (KPoly<Rational>, string) e1)
    {
        var (_, _, t) = Ring.Polynomial(e0.Item2, e1.Item2, "_t_", Rational.KZero());
        var a = e0.Item1.ToPolynomial(t.Indeterminates, t.Indeterminates[0]);
        var b = e1.Item1.ToPolynomial(t.Indeterminates, t.Indeterminates[1]);
        var nbf = NumberFieldQ(new[] { a, b });
        return (nbf[0], nbf[1]);
    }

    public static (EPolynomial<Rational> e0, EPolynomial<Rational> e1, EPolynomial<Rational> p0) NumberFieldQ(
        (KPoly<Rational>, string) e0,
        (KPoly<Rational>, string) e1, string p0, params string[] others)
    {
        var all = new[] { e0.Item2, e1.Item2, p0, "_t_" }.ToArray();
        var xis = Ring.Polynomial(Rational.KZero(), MonomOrder.Lex, all);
        var x0 = xis[0];
        var a = e0.Item1.ToPolynomial(x0.Indeterminates, x0.Indeterminates[0]);
        var b = e1.Item1.ToPolynomial(x0.Indeterminates, x0.Indeterminates[1]);
        var nbf = NumberFieldQ(new[] { a, b });
        return (nbf[0], nbf[1], new EPolynomial<Rational>(xis[2], nbf[0].Basis));
    }

    // Algebre Tome 1, page 415
    public static KPoly<Rational>[] GaloisGroupPolynomialsList(int n)
    {
        var x = QPoly();
        if (n == 1)
            return new[] { x };
        if (n == 2)
            return new[] { x.Pow(2) + x + 1 };
        if (n == 3)
            return new[]
            {
                x.Pow(3) + x.Pow(2) - 2 * x - 1,
                x.Pow(3) + 2
            };
        if (n == 4)
            return new[]
            {
                x.Pow(4) + x.Pow(3) + x.Pow(2) + x + 1,
                x.Pow(4) + 1, x.Pow(4) - 2,
                x.Pow(4) + 8 * x + 12,
                x.Pow(4) + x + 1
            };
        if (n == 5)
            return new[]
            {
                x.Pow(5) + x.Pow(4) - 4 * x.Pow(3) - 3 * x.Pow(2) + 3 * x + 1,
                x.Pow(5) - 5 * x + 12,
                x.Pow(5) + 2,
                x.Pow(5) + 20 * x + 16,
                x.Pow(5) - x + 1
            };
        if (n == 6)
            return new[]
            {
                x.Pow(6) + x.Pow(5) + x.Pow(4) + x.Pow(3) + x.Pow(2) + x + 1, x.Pow(6) + 108, x.Pow(6) + 2,
                x.Pow(6) - 3 * x.Pow(2) - 1, x.Pow(6) + 3 * x.Pow(3) + 3,
                x.Pow(6) - 3 * x.Pow(2) + 1, x.Pow(6) - 4 * x.Pow(2) - 1,
                x.Pow(6) - 3 * x.Pow(5) + 6 * x.Pow(4) - 7 * x.Pow(3) + 2 * x.Pow(2) + x - 4,
                x.Pow(6) + 2 * x.Pow(3) - 2,
                x.Pow(6) + 6 * x.Pow(4) + 2 * x.Pow(3) + 9 * x.Pow(2) + 6 * x - 4,
                x.Pow(6) + 2 * x.Pow(2) + 2,
                x.Pow(6) - 2 * x.Pow(5) - 5 * x.Pow(2) - 2 * x - 1,
                x.Pow(6) + 2 * x.Pow(4) + 2 * x.Pow(3) + x.Pow(2) + 2 * x + 2,
                x.Pow(6) - x.Pow(5) - 10 * x.Pow(4) + 30 * x.Pow(3) - 31 * x.Pow(2) + 7 * x + 9,
                x.Pow(6) + 24 * x - 20,
                x.Pow(6) + x + 1
            };
        if (n == 7)
            return new[]
            {
                x.Pow(7) + x.Pow(6) - 12 * x.Pow(5) - 7 * x.Pow(4) + 28 * x.Pow(3) + 14 * x.Pow(2) - 9 * x + 1,
                x.Pow(7) + 7 * x.Pow(3) + 7 * x.Pow(2) + 7 * x - 1,
                x.Pow(7) - 14 * x.Pow(5) + 56 * x.Pow(3) - 56 * x + 22,
                x.Pow(7) + 2,
                x.Pow(7) - 7 * x.Pow(3) + 14 * x.Pow(2) - 7 * x + 1, x.Pow(7) + 7 * x.Pow(4) + 14 * x + 3, x.Pow(7) + x + 1
            };

        throw new();
    }
}
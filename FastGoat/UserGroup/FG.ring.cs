using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Matrix;
using FastGoat.UserGroup.Padic;

namespace FastGoat.UserGroup;

public static partial class FG
{
    public static KPoly<ZnInt> ZPoly(int p, char x = 'x') => new KPoly<ZnInt>(x, ZnInt.KZero(p)).X;
    public static KPoly<ZnBInt> ZbPoly(int p, char x = 'x') => new KPoly<ZnBInt>(x, ZnBInt.KZero(p)).X;
    public static KPoly<Rational> QPoly(char x = 'x') => new KPoly<Rational>(x);

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

    public static (KPoly<EPoly<K>> X,EPoly<K> c) EPolyXc<K>(KPoly<K> f, char c, char x = 'X')
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

    public static FracPoly<ZnInt> ZFracPoly(int p, char x = 'x') => KFracPoly(x, ZnInt.KZero(p));

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
    
    public static (EPolynomial<Rational> e0, EPolynomial<Rational> e1, EPolynomial<Rational> p0) NumberFieldQ((KPoly<Rational>, string) e0,
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
}
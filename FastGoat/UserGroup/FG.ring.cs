using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Integers;

namespace FastGoat.UserGroup;

public static partial class FG
{
    public static KPoly<ZnInt> ZPoly(int p, char x = 'x') => new KPoly<ZnInt>(x, ZnInt.KZero(p)).X;
    public static KPoly<Rational> QPoly(char x = 'x') => new KPoly<Rational>(x);
    public static KPoly<Padic> PadicPoly(int p, int o, char x = 'x') => new KPoly<Padic>(x, new Padic(p, o)).X;

    public static KPoly<K> KPoly<K>(char x, K scalar) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        return new KPoly<K>(x, scalar).X;
    }

    public static EPoly<K> EPoly<K>(KPoly<K> f) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        return new(f);
    }

    public static Fraction<KPoly<K>> FracPoly<K>(char x, K scalar) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var x0 = KPoly(x, scalar);
        var t0 = new Fraction<KPoly<K>>(x0);
        return t0;
    }

    public static Fraction<KPoly<Rational>> QFracPoly(char x = 'x') => FracPoly(x, Rational.KZero());

    public static Fraction<KPoly<ZnInt>> ZFracPoly(int p, char x = 'x') => FracPoly(x, ZnInt.KZero(p));

    public static (KPoly<Fraction<KPoly<ZnInt>>> x, Fraction<KPoly<ZnInt>> t) FpT_Poly(int p, (char x, char t) xt)
    {
        var t = ZPoly(p, xt.t);
        var t0 = new Fraction<KPoly<ZnInt>>(t);
        return (KPoly(xt.x, t0), t0);
    }

    public static (KPoly<Fraction<KPoly<ZnInt>>> x, Fraction<KPoly<ZnInt>>t) FpT_Poly(int p) => FpT_Poly(p, ('x', 't'));

    public static (KPoly<EPoly<ZnInt>> x, EPoly<ZnInt> a) FqX_Poly(int q, (char x, char a) xa)
    {
        var a = FqX(q, xa.a);
        return (KPoly(xa.x, a), a);
    }

    public static (KPoly<EPoly<ZnInt>> x, EPoly<ZnInt> a) FqX_Poly(int q) => FqX_Poly(q, ('x', 'a'));
}
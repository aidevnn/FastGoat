using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Integers;

namespace FastGoat.UserGroup;

public static partial class FG
{
    public static KPoly<ZnInt> ZPoly(int p, char x = 'x') => new KPoly<ZnInt>(x, ZnInt.KZero(p)).X;
    public static KPoly<Rational> QPoly(char x = 'x') => new KPoly<Rational>(x);

    public static KPoly<K> KPoly<K>(char x, K scalar) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        return new KPoly<K>(x, scalar).X;
    }

    public static (KPoly<Fraction<KPoly<ZnInt>>> x, Fraction<KPoly<ZnInt>>t) FpT_Poly(int p, (char x, char t) xt)
    {
        var t = ZPoly(p, xt.t);
        var t0 = new Fraction<KPoly<ZnInt>>(t);
        return (KPoly(xt.x, t0), t0);
    }

    public static (KPoly<Fraction<KPoly<ZnInt>>> x, Fraction<KPoly<ZnInt>>t) FpT_Poly(int p) => FpT_Poly(p, ('x', 't'));
}
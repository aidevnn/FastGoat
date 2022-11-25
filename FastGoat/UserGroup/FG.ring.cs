using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Integers;

namespace FastGoat.UserGroup;

public static partial class FG
{
    public static KPoly<ZnInt> ZPoly(int p, char x = 'x') => new KPoly<ZnInt>(x, ZnInt.KZero(p)).X;
    public static KPoly<Rational> QPoly(char x = 'x') => new KPoly<Rational>(x);
}
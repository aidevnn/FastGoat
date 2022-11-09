using FastGoat.Commons;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Integers;

namespace FastGoat.Structures;

public static class Ring
{
    public static IEnumerable<T> TrimSeq<T>(this IEnumerable<T> seq) where T : struct, IElt<T>, IRingElt<T>
    {
        var stack = new Stack<T>(seq);
        while (stack.Count != 0 && stack.Peek().IsZero())
            stack.Pop();

        return stack.Reverse();
    }

    public static T Gcd<T>(T a, T b) where T : struct, IElt<T>, IRingElt<T>
    {
        if (b.IsZero())
            return a.CompareTo(a.Opp()) == -1 ? a.Opp() : a;

        return Gcd(b, a.Div(b).rem);
    }

    public static (T x, T y) Bezout<T>(T a, T b) where T : struct, IElt<T>, IRingElt<T>
    {
        if (b.IsZero())
            return a.CompareTo(a.Opp()) == -1 ? (a.One.Opp(), a.Zero) : (a.One, a.Zero);

        var (q, r) = a.Div(b);
        var (x0, y0) = Bezout(b, r);
        return (y0, x0.Add(q.Mul(y0).Opp()));
    }

    public static KPoly<ZnInt> ZPolynomial(int p = 0, char x = 'x')
    {
        if (p != 0 && !IntExt.Primes10000.Contains(p))
            throw new GroupException(GroupExceptionType.GroupDef);

        var kZero = new ZnInt(p, 0);
        return new KPoly<ZnInt>(x, kZero, new[] { kZero, kZero.One });
    }

    public static KPoly<Rational> QPolynomial(char x = 'x')
    {
        var kZero = new Rational(0, 1);
        return new KPoly<Rational>(x, kZero, new[] { kZero, kZero.One });
    }

    public static EPoly<T> EPoly<T>(KPoly<T> f, char x)
        where T : struct, IElt<T>, IRingElt<T>, IFieldElt<T>
    {
        var f0 = new KPoly<T>(x, f.KZero, f.Coefs);
        return new EPoly<T>(f0);
    }

    public static (KPoly<EPoly<T>>, EPoly<T>) ExtPolynomial<T>(KPoly<T> f, char x)
        where T : struct, IElt<T>, IRingElt<T>, IFieldElt<T>
    {
        var f0 = new KPoly<T>(x, f.KZero, f.Coefs);
        var fx = new EPoly<T>(f0);
        var kx = new KPoly<EPoly<T>>(f.x, fx.Zero, new[] { fx.Zero, fx.One });

        return (kx, fx.X);
    }

    public static Polynomial<K, Xi>[] Polynomial<K>(K zero, char x0, params char[] xi)
        where K : struct, IFieldElt<K>, IElt<K>, IRingElt<K>
    {
        return new Polynomial<K, Xi>(xi.Prepend(x0).Select(c => new Xi(c)), zero).Xi();
    }

    public static Polynomial<K, Xi> Polynomial<K>(char x1, K zero)
        where K : struct, IFieldElt<K>, IElt<K>, IRingElt<K>
    {
        return Polynomial(zero, x1)[0];
    }

    public static Polynomial<K, Xi> Polynomial<K>(K zero)
        where K : struct, IFieldElt<K>, IElt<K>, IRingElt<K>
    {
        return Polynomial(zero, 'X')[0];
    }

    public static (Polynomial<K, Xi> x1, Polynomial<K, Xi> x2) Polynomial<K>(char x1, char x2, K zero)
        where K : struct, IFieldElt<K>, IElt<K>, IRingElt<K>
    {
        var xi = Polynomial(zero, x1, x2);
        return (xi[0], xi[1]);
    }

    public static (Polynomial<K, Xi> x1, Polynomial<K, Xi> x2, Polynomial<K, Xi> x3) 
        Polynomial<K>(char x1, char x2, char x3, K zero)
        where K : struct, IFieldElt<K>, IElt<K>, IRingElt<K>
    {
        var xi = Polynomial(zero, x1, x2, x3);
        return (xi[0], xi[1], xi[2]);
    }

}
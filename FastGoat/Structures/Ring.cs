using FastGoat.Commons;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Integers;

namespace FastGoat.Structures;

public static partial class Ring
{
    public static IEnumerable<T> TrimSeq<T>(this IEnumerable<T> seq) where T : IElt<T>, IRingElt<T>
    {
        var stack = new Stack<T>(seq);
        while (stack.Count != 0 && stack.Peek().IsZero())
            stack.Pop();

        return stack.Reverse();
    }

    public static T Gcd<T>(T a, T b) where T : IElt<T>, IRingElt<T>
    {
        if (b.IsZero())
            return a.CompareTo(a.Opp()) == -1 ? a.Opp() : a;

        return Gcd(b, a.Div(b).rem);
    }

    public static (T x, T y) Bezout<T>(T a, T b) where T : IElt<T>, IRingElt<T>
    {
        if (b.IsZero())
            return a.CompareTo(a.Opp()) == -1 ? (a.One.Opp(), a.Zero) : (a.One, a.Zero);

        var (q, r) = a.Div(b);
        var (x0, y0) = Bezout(b, r);
        return (y0, x0.Add(q.Mul(y0).Opp()));
    }
    
    public static Monom<Xi> Xi(char c) => new Monom<Xi>(new Xi(c));

    public static Monom<T> Monom<T>(T e) where T : struct, IElt<T>
    {
        return new Monom<T>(e);
    }

    public static Polynomial<K, Xi>[] Polynomial<K>(K zero, char[] xi)
        where K : struct, IFieldElt<K>, IElt<K>, IRingElt<K>
    {
        return xi.Select(c => Polynomial(c, zero)).ToArray();
    }

    public static Polynomial<K, Xi>[] Polynomial<K>(K zero, char x0, params char[] xi)
        where K : struct, IFieldElt<K>, IElt<K>, IRingElt<K>
    {
        return Polynomial(zero, xi.Prepend(x0).ToArray());
    }

    public static Polynomial<K, Xi> PolynomialZero<K>(K k0)
        where K : struct, IFieldElt<K>, IElt<K>, IRingElt<K>
    {
        return new Polynomial<K, Xi>(k0.Zero);
    }

    public static Polynomial<K, Xi> Polynomial<K>(char x1, K zero)
        where K : struct, IFieldElt<K>, IElt<K>, IRingElt<K>
    {
        return new Polynomial<K, Xi>(new Xi(x1), zero.One);
    }

    public static Polynomial<K, Xi> Polynomial<K>(K zero)
        where K : struct, IFieldElt<K>, IElt<K>, IRingElt<K>
    {
        return Polynomial('X', zero);
    }

    public static (Polynomial<K, Xi> x1, Polynomial<K, Xi> x2) Polynomial<K>(char x1, char x2, K zero)
        where K : struct, IFieldElt<K>, IElt<K>, IRingElt<K>
    {
        return (Polynomial(x1, zero), Polynomial(x2, zero));
    }

    public static (Polynomial<K, Xi> x1, Polynomial<K, Xi> x2, Polynomial<K, Xi> x3)
        Polynomial<K>(char x1, char x2, char x3, K zero)
        where K : struct, IFieldElt<K>, IElt<K>, IRingElt<K>
    {
        return (Polynomial(x1, zero), Polynomial(x2, zero), Polynomial(x3, zero));
    }

    public static (int d, Polynomial<ZnInt, Xi> poly)[] ConwayPolys(int p)
    {
        var all = PolynomExt.GetAll(p);
        var x = Polynomial(ZnInt.KZero(p));
        return all.Select(e => (e.d, e.coefs.Select((a, i) => a * x.Pow(i)).Aggregate((a, b) => a.Add(b)))).ToArray();
    }

    public static KPoly<Rational> QX(char c = 'X')
    {
        return new KPoly<Rational>(c);
    }
    public static KPoly<ZnInt> ZX(char c = 'X', int p = 0)
    {
        return new KPoly<ZnInt>(c, ZnInt.KZero(p));
    }
}
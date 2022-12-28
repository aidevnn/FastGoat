using FastGoat.Commons;
using FastGoat.Structures.VecSpace;

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

    public static T Gcd<T>(T[] arr)where T : IElt<T>, IRingElt<T>
    {
        if (arr.Length == 0)
            throw new ArgumentException();
        
        if (arr.Length == 1)
            return arr[0];

        return Gcd(arr.First(), Gcd(arr.Skip(1).ToArray()));
    }

    public static (T x, T y) Bezout<T>(T a, T b) where T : IElt<T>, IRingElt<T>
    {
        if (b.IsZero())
            return a.CompareTo(a.Opp()) == -1 ? (a.One.Opp(), a.Zero) : (a.One, a.Zero);

        var (q, r) = a.Div(b);
        var (x0, y0) = Bezout(b, r);
        return (y0, x0.Add(q.Mul(y0).Opp()));
    }
    
    public static Monom<Xi> Xi(char c, int n = 1) => new Monom<Xi>(new Xi(c), n);

    public static Indeterminates<Xi> Indeterminates(params string[] xs) => new(xs.Select(s => new Xi(s)).ToArray());

    public static Indeterminates<Xi> Indeterminates(params char[] xs) =>
        Indeterminates(xs.Select(c => c.ToString()).ToArray());
    public static Polynomial<K, Xi> Polynomial<K>(K zero, Indeterminates<Xi> indeterminates)
        where K : struct, IFieldElt<K>, IElt<K>, IRingElt<K>
    {
        return new(indeterminates, zero);
    }
    public static Polynomial<K, Xi>[] Polynomial<K>(K zero, string[] xi)
        where K : struct, IFieldElt<K>, IElt<K>, IRingElt<K>
    {
        var indeterminates = Indeterminates(xi);
        return xi.Select(c => new Polynomial<K,Xi>(new Monom<Xi>(indeterminates, new(c), 1), zero.One)).ToArray();
    }

    public static Polynomial<K, Xi>[] Polynomial<K>(K zero, char[] xi)
        where K : struct, IFieldElt<K>, IElt<K>, IRingElt<K>
    {
        return Polynomial(zero, xi.Select(c => $"{c}").ToArray());
    }

    public static Polynomial<K, Xi>[] Polynomial<K>(K zero, string x0, params string[] xi)
        where K : struct, IFieldElt<K>, IElt<K>, IRingElt<K>
    {
        return Polynomial(zero, xi.Prepend(x0).ToArray());
    }

    public static Polynomial<K, Xi> Polynomial<K>(K zero)
        where K : struct, IFieldElt<K>, IElt<K>, IRingElt<K>
    {
        return Polynomial(zero, "X")[0];
    }

    public static (Polynomial<K, Xi> x1, Polynomial<K, Xi> x2) Polynomial<K>(string x1, string x2, K zero)
        where K : struct, IFieldElt<K>, IElt<K>, IRingElt<K>
    {
        var polys = Polynomial(zero, x1, x2);
        return (polys[0], polys[1]);
    }

    public static (Polynomial<K, Xi> x1, Polynomial<K, Xi> x2, Polynomial<K, Xi> x3)
        Polynomial<K>(string x1, string x2, string x3, K zero)
        where K : struct, IFieldElt<K>, IElt<K>, IRingElt<K>
    {
        var polys = Polynomial(zero, x1, x2, x3);
        return (polys[0], polys[1], polys[2]);
    }
    
    public static (Polynomial<K, Xi> x1, Polynomial<K, Xi> x2, Polynomial<K, Xi> x3, Polynomial<K, Xi> x4)
        Polynomial<K>(string x1, string x2, string x3, string x4, K zero)
        where K : struct, IFieldElt<K>, IElt<K>, IRingElt<K>
    {
        var polys = Polynomial(zero, x1, x2, x3, x4);
        return (polys[0], polys[1], polys[2], polys[3]);
    }
}
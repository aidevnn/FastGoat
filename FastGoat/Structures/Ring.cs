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

        if (a.CompareTo(b) == -1)
        {
            var a0 = b.Div(a).rem;
            if (a0.Equals(b))
            {
                Console.WriteLine(new { a, b, a0 });
                return a.One;
            }

            return Gcd(a, a0);
        }
        else
        {
            return Gcd(b, a.Div(b).rem);
        }
    }

    public static T Gcd<T>(T[] arr) where T : IElt<T>, IRingElt<T>
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

    public static Polynomial<K, Xi> Polynomial<K>(K scalar, Indeterminates<Xi> indeterminates)
        where K : struct, IFieldElt<K>, IElt<K>, IRingElt<K>
    {
        return new(indeterminates, scalar);
    }

    public static Polynomial<K, Xi>[] Polynomial<K>(K scalar, string[] xi)
        where K : struct, IFieldElt<K>, IElt<K>, IRingElt<K>
    {
        var indeterminates = Indeterminates(xi);
        return xi.Select(c => new Polynomial<K, Xi>(new Monom<Xi>(indeterminates, new(c), 1), scalar.One)).ToArray();
    }

    public static Polynomial<K, Xi>[] Polynomial<K>(K scalar, char[] xi)
        where K : struct, IFieldElt<K>, IElt<K>, IRingElt<K>
    {
        return Polynomial(scalar, xi.Select(c => $"{c}").ToArray());
    }

    public static Polynomial<K, Xi>[] Polynomial<K>(K scalar, string x0, params string[] xi)
        where K : struct, IFieldElt<K>, IElt<K>, IRingElt<K>
    {
        return Polynomial(scalar, xi.Prepend(x0).ToArray());
    }

    public static Polynomial<K, Xi> Polynomial<K>(string x, K scalar)
        where K : struct, IFieldElt<K>, IElt<K>, IRingElt<K>
    {
        return Polynomial(scalar, x)[0];
    }

    public static Polynomial<K, Xi> Polynomial<K>(K scalar)
        where K : struct, IFieldElt<K>, IElt<K>, IRingElt<K>
    {
        return Polynomial("X", scalar);
    }

    public static (Polynomial<K, Xi> x1, Polynomial<K, Xi> x2) Polynomial<K>(string x1, string x2, K scalar)
        where K : struct, IFieldElt<K>, IElt<K>, IRingElt<K>
    {
        var polys = Polynomial(scalar, x1, x2);
        return (polys[0], polys[1]);
    }

    public static (Polynomial<K, Xi> x1, Polynomial<K, Xi> x2, Polynomial<K, Xi> x3)
        Polynomial<K>(string x1, string x2, string x3, K scalar)
        where K : struct, IFieldElt<K>, IElt<K>, IRingElt<K>
    {
        var polys = Polynomial(scalar, x1, x2, x3);
        return (polys[0], polys[1], polys[2]);
    }

    public static (Polynomial<K, Xi> x1, Polynomial<K, Xi> x2, Polynomial<K, Xi> x3, Polynomial<K, Xi> x4)
        Polynomial<K>(string x1, string x2, string x3, string x4, K scalar)
        where K : struct, IFieldElt<K>, IElt<K>, IRingElt<K>
    {
        var polys = Polynomial(scalar, x1, x2, x3, x4);
        return (polys[0], polys[1], polys[2], polys[3]);
    }

    public static (Polynomial<K, Xi> x1, Polynomial<K, Xi> x2, Polynomial<K, Xi> x3, Polynomial<K, Xi> x4, Polynomial<K, Xi> x5)
        Polynomial<K>(string x1, string x2, string x3, string x4, string x5, K scalar)
        where K : struct, IFieldElt<K>, IElt<K>, IRingElt<K>
    {
        var polys = Polynomial(scalar, x1, x2, x3, x4, x5);
        return (polys[0], polys[1], polys[2], polys[3], polys[4]);
    }
    
    public static KPoly<K> ToKPoly<K>(this Polynomial<K, Xi> f, Xi x)
        where K : struct, IFieldElt<K>, IElt<K>, IRingElt<K>
    {
        var d = f.DegreeOf(x);
        var coefs = (d + 1).Range().Select(i => f[new(f.Indeterminates, x, i)]).ToArray();
        return new KPoly<K>(x.xi[0], f.KZero, coefs);
    }

    public static KPoly<K> ToKPoly<K>(this Polynomial<K, Xi> f, Polynomial<K, Xi> x)
        where K : struct, IFieldElt<K>, IElt<K>, IRingElt<K>
    {
        return ToKPoly(f, x.ExtractIndeterminate);
    }

    public static Polynomial<K, Xi> ToPolynomial<K>(this KPoly<K> f, Indeterminates<Xi> indeterminates, Xi xi)
        where K : struct, IFieldElt<K>, IElt<K>, IRingElt<K>
    {
        var mnm = new Monom<Xi>(indeterminates, xi, 1);
        return f.Coefs.Select((k, i) => new Polynomial<K, Xi>(mnm.Pow(i), k)).Aggregate((a, b) => a + b);
    }
}
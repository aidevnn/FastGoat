using System.Numerics;
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
            if(a.IsZero())
                return b.CompareTo(b.Opp()) == -1 ? b.Opp() : b;
            
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

    public static Polynomial<K, Xi>[] Polynomial<K>(K scalar, string[] xi)
        where K : struct, IFieldElt<K>, IElt<K>, IRingElt<K>
    {
        var indeterminates = Indeterminates(xi);
        return xi.Select(c => new Polynomial<K, Xi>(new Monom<Xi>(indeterminates, new(c), 1), scalar.One)).ToArray();
    }

    public static Polynomial<K, Xi>[] Polynomial<K>(K scalar, Indeterminates<Xi> indeterminates)
        where K : struct, IFieldElt<K>, IElt<K>, IRingElt<K>
    {
        return indeterminates.Content.Select(c => new Polynomial<K, Xi>(new Monom<Xi>(indeterminates, c, 1), scalar.One)).ToArray();
    }

    public static Polynomial<K, Xi>[] Polynomial<K>(K scalar, MonomOrder order, string[] xi)
        where K : struct, IFieldElt<K>, IElt<K>, IRingElt<K>
    {
        var indeterminates = new Indeterminates<Xi>(order, xi.Select(x => new Xi(x)).ToArray());
        return xi.Select(c => new Polynomial<K, Xi>(new Monom<Xi>(indeterminates, new(c), 1), scalar.One)).ToArray();
    }

    public static (Polynomial<K, Xi> X,Polynomial<K, Xi>[] Xis) Polynomial<K>(K scalar, MonomOrder order, (int n, string h) e, string x0)
        where K : struct, IFieldElt<K>, IElt<K>, IRingElt<K>
    {
        var digits = $"{e.n - 1}".Length;
        var zeros = Enumerable.Repeat(0, digits).Glue();
        var fmt = $"{{0}}{{1,{digits}:{zeros}}}";
        var xi = e.n.Range().Select(i => string.Format(fmt, e.h, i)).Append(x0).ToArray();
        var indeterminates = new Indeterminates<Xi>(order, xi.Select(x => new Xi(x)).ToArray());
        var xis = xi.Select(c => new Polynomial<K, Xi>(new Monom<Xi>(indeterminates, new(c), 1), scalar.One)).ToArray();
        var X0 = xis.Last();
        var Xis = xis.SkipLast(1).ToArray();
        return (X0, Xis);
    }

    public static Polynomial<K, Xi>[] Polynomial<K>(K scalar, MonomOrder order, (int n, string h) e)
        where K : struct, IFieldElt<K>, IElt<K>, IRingElt<K>
    {
        var digits = $"{e.n - 1}".Length;
        var zeros = Enumerable.Repeat(0, digits).Glue();
        var fmt = $"{{0}}{{1,{digits}:{zeros}}}";
        var xi = e.n.Range().Select(i => string.Format(fmt, e.h, i)).ToArray();
        var indeterminates = new Indeterminates<Xi>(order, xi.Select(x => new Xi(x)).ToArray());
        return xi.Select(c => new Polynomial<K, Xi>(new Monom<Xi>(indeterminates, new(c), 1), scalar.One)).ToArray();
    }

    public static EPolynomial<K>[] EPolynomial<K>(K scalar, MonomOrder order, (int n, string h) e)
        where K : struct, IFieldElt<K>, IElt<K>, IRingElt<K>
    {
        var xis = Polynomial(scalar, order, e);
        var basis = new PolynomialBasis<K, Xi>(xis[0].Indeterminates);
        return xis.Select(xi => new EPolynomial<K>(xi, basis)).ToArray();
    }

    public static EPolynomial<K>[] EPolynomial<K>(K scalar, MonomOrder order, string[] xis)
        where K : struct, IFieldElt<K>, IElt<K>, IRingElt<K>
    {
        var xis0 = Polynomial(scalar, order, xis);
        var basis = new PolynomialBasis<K, Xi>(xis0[0].Indeterminates);
        return xis0.Select(xi => new EPolynomial<K>(xi, basis)).ToArray();
    }

    public static (EPolynomial<K> X,EPolynomial<K>[] xis) EPolynomial<K>(K scalar, MonomOrder order, (int n, string h) e, string x0)
        where K : struct, IFieldElt<K>, IElt<K>, IRingElt<K>
    {
        var (x, xis) = Polynomial(scalar, order, e, x0);
        var basis = new PolynomialBasis<K, Xi>(x.Indeterminates);
        return (new EPolynomial<K>(x, basis), xis.Select(xi => new EPolynomial<K>(xi, basis)).ToArray());
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

    public static (Polynomial<K, Xi> x1, Polynomial<K, Xi> x2)
        Polynomial<K>(string x1, string x2, K scalar, MonomOrder order = MonomOrder.GrLex)
        where K : struct, IFieldElt<K>, IElt<K>, IRingElt<K>
    {
        var polys = Polynomial(scalar, order, new[] { x1, x2 });
        return (polys[0], polys[1]);
    }

    public static (Polynomial<K, Xi> x1, Polynomial<K, Xi> x2, Polynomial<K, Xi> x3)
        Polynomial<K>(string x1, string x2, string x3, K scalar, MonomOrder order = MonomOrder.GrLex)
        where K : struct, IFieldElt<K>, IElt<K>, IRingElt<K>
    {
        var polys = Polynomial(scalar, order, new[] { x1, x2, x3 });
        return (polys[0], polys[1], polys[2]);
    }

    public static (Polynomial<K, Xi> x1, Polynomial<K, Xi> x2, Polynomial<K, Xi> x3, Polynomial<K, Xi> x4) Polynomial<K>(string x1,
        string x2, string x3, string x4, K scalar, MonomOrder order = MonomOrder.GrLex)
        where K : struct, IFieldElt<K>, IElt<K>, IRingElt<K>
    {
        var polys = Polynomial(scalar, order, new[] { x1, x2, x3, x4 });
        return (polys[0], polys[1], polys[2], polys[3]);
    }

    public static (Polynomial<K, Xi> x1, Polynomial<K, Xi> x2, Polynomial<K, Xi> x3, Polynomial<K, Xi> x4, Polynomial<K, Xi> x5)
        Polynomial<K>(string x1, string x2, string x3, string x4, string x5, K scalar, MonomOrder order = MonomOrder.GrLex)
        where K : struct, IFieldElt<K>, IElt<K>, IRingElt<K>
    {
        var polys = Polynomial(scalar, order, new[] { x1, x2, x3, x4, x5 });
        return (polys[0], polys[1], polys[2], polys[3], polys[4]);
    }

    public static (Polynomial<K, Xi> x1, Polynomial<K, Xi> x2, Polynomial<K, Xi> x3,
        Polynomial<K, Xi> x4, Polynomial<K, Xi> x5, Polynomial<K, Xi> x6)
        Polynomial<K>(string x1, string x2, string x3, string x4, string x5, string x6, K scalar, MonomOrder order = MonomOrder.GrLex)
        where K : struct, IFieldElt<K>, IElt<K>, IRingElt<K>
    {
        var polys = Polynomial(scalar, order, new[] { x1, x2, x3, x4, x5, x6 });
        return (polys[0], polys[1], polys[2], polys[3], polys[4], polys[5]);
    }
    
    
    public static (EPolynomial<K> x1, EPolynomial<K> x2)
        EPolynomial<K>(string x1, string x2, K scalar, MonomOrder order = MonomOrder.GrLex)
        where K : struct, IFieldElt<K>, IElt<K>, IRingElt<K>
    {
        var polys = EPolynomial(scalar, order, new[] { x1, x2 });
        return (polys[0], polys[1]);
    }

    public static (EPolynomial<K> x1, EPolynomial<K> x2, EPolynomial<K> x3)
        EPolynomial<K>(string x1, string x2, string x3, K scalar, MonomOrder order = MonomOrder.GrLex)
        where K : struct, IFieldElt<K>, IElt<K>, IRingElt<K>
    {
        var polys = EPolynomial(scalar, order, new[] { x1, x2, x3 });
        return (polys[0], polys[1], polys[2]);
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

    public static Polynomial<K, Xi> ToPolynomial<K>(this KPoly<K> f, Polynomial<K, Xi> x)
        where K : struct, IFieldElt<K>, IElt<K>, IRingElt<K>
    {
        return f.Coefs.Select((k, i) => k * x.Pow(i)).Aggregate(x.Zero, (a, b) => a + b);
    }

    public static Polynomial<EPoly<K>, Xi> ToPolynomial<K>(this KPoly<K> f, Polynomial<EPoly<K>, Xi> x)
        where K : struct, IFieldElt<K>, IElt<K>, IRingElt<K>
    {
        return f.Coefs.Select((k, i) => k * x.KOne * x.Pow(i)).Aggregate(x.Zero, (a, b) => a + b);
    }

    public static KPoly<EPoly<K>> SubstituteP0b<K>(this KPoly<EPoly<K>> P, KPoly<EPoly<K>> s, EPoly<K> a)
        where K : struct, IFieldElt<K>, IElt<K>, IRingElt<K>
    {
        return P.Coefs.Select((c, i) => c.Poly.Substitute(a) * s.Pow(i)).Aggregate(s.Zero, (acc, c) => acc + c);
    }

    public static Complex Substitute(this KPoly<Rational> f, Complex c) => f.Coefs.Select((k, i) => k * Complex.Pow(c, i))
        .Aggregate(Complex.Zero, (sum, ci) => sum + ci);
}
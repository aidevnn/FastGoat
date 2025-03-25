using System.Numerics;
using FastGoat.Commons;
using FastGoat.Structures.VecSpace;

namespace FastGoat.Structures;

public static partial class Ring
{
    public static bool IsOne<T>(this T t) where T : IElt<T>, IRingElt<T>
    {
        return (t - t.One).IsZero();
    }

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

    public static T Lcm<T>(T a, T b) where T : IElt<T>, IRingElt<T>
    {
        return (a * b) / Gcd(a, b);
    }

    public static T Gcd<T>(T[] arr) where T : IElt<T>, IRingElt<T>
    {
        if (arr.Length == 0)
            throw new ArgumentException();

        if (arr.Length == 1)
            return arr[0];

        return Gcd(arr.First(), Gcd(arr.Skip(1).ToArray()));
    }

    public static T Lcm<T>(T[] arr) where T : IElt<T>, IRingElt<T>
    {
        if (arr.Length == 0)
            throw new ArgumentException();

        if (arr.Length == 1)
            return arr[0];

        return Lcm(arr.First(), Lcm(arr.Skip(1).ToArray()));
    }

    public static (T x, T y) Bezout<T>(T a, T b) where T : IElt<T>, IRingElt<T>
    {
        if (b.IsZero())
            return a.CompareTo(a.Opp()) == -1 ? (a.One.Opp(), a.Zero) : (a.One, a.Zero);

        var (q, r) = a.Div(b);
        var (x0, y0) = Bezout(b, r);
        return (y0, x0.Add(q.Mul(y0).Opp()));
    }

    public static K FastPow<K>(this K a, int k) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        if (k == 0)
            return a.One;

        if (k < 0)
            return a.Inv().FastPow(-k);
        
        var (r, a0, e0) = (a.One, a, k);
        while (e0 > 0)
        {
            if (e0 % 2 == 1)
                r *= a0;

            e0 >>= 1;
            a0 *= a0;
        }

        return r;
    }

    public static Monom<Xi> Xi(string c, int n = 1) => new Monom<Xi>(new Xi(c), n);

    public static Indeterminates<Xi> Indeterminates(params string[] xs) => new(xs.Select(s => new Xi(s)).ToArray());

    public static Indeterminates<Xi> Indeterminates(params char[] xs) =>
        Indeterminates(xs.Select(c => c.ToString()).ToArray());

    public static Polynomial<K, Xi> ToPolynomial<K>(this Monom<Xi> xi, Polynomial<K, Xi> f)
        where K : struct, IFieldElt<K>, IElt<K>, IRingElt<K>
    {
        if (!xi.Indeterminates.Equals(f.Indeterminates))
            throw new();
        return new Polynomial<K, Xi>(xi, f.KOne) * f;
    }

    public static Polynomial<K, Xi> ToPolynomial<K>(this Xi xi, Polynomial<K, Xi> f)
        where K : struct, IFieldElt<K>, IElt<K>, IRingElt<K>
    {
        return new Polynomial<K, Xi>(new Monom<Xi>(f.Indeterminates, xi), f.KOne) * f;
    }

    public static Polynomial<K, Xi>[] Polynomial<K>(K scalar, params string[] xi)
        where K : struct, IFieldElt<K>, IElt<K>, IRingElt<K>
    {
        if (xi.Length == 0)
            throw new();

        var indeterminates = Indeterminates(xi);
        return xi.Select(c => new Polynomial<K, Xi>(new Monom<Xi>(indeterminates, new(c), 1), scalar.One)).ToArray();
    }

    public static Polynomial<K, Xi>[] Polynomial<K>(K scalar, MonomOrder order, params string[] xi)
        where K : struct, IFieldElt<K>, IElt<K>, IRingElt<K>
    {
        if (xi.Length == 0)
            throw new();

        var indeterminates = new Indeterminates<Xi>(order, xi.Select(x => new Xi(x)).ToArray());
        return xi.Select(c => new Polynomial<K, Xi>(new Monom<Xi>(indeterminates, new(c), 1), scalar.One)).ToArray();
    }

    public static Polynomial<K, Xi>[] Polynomial<K>(K scalar, char[] xi)
        where K : struct, IFieldElt<K>, IElt<K>, IRingElt<K>
    {
        return Polynomial(scalar, xi.Select(c => $"{c}").ToArray());
    }

    public static Polynomial<K, Xi> Polynomial<K>(K scalar)
        where K : struct, IFieldElt<K>, IElt<K>, IRingElt<K>
    {
        return Polynomial(scalar, new[] { "X" })[0];
    }

    public static (Polynomial<K, Xi> X, Polynomial<K, Xi>[] Xis) Polynomial<K>(K scalar, MonomOrder order, (int n, string h) e,
        string x0)
        where K : struct, IFieldElt<K>, IElt<K>, IRingElt<K>
    {
        var digits = $"{e.n - 1}".Length;
        var zeros = Enumerable.Repeat(0, digits).Glue();
        var fmt = $"{{0}}{{1,{digits}:{zeros}}}";
        var lex = (order & MonomOrder.Lex) == MonomOrder.Lex;
        var xi = e.n.Range().Select(i => string.Format(fmt, e.h, i)).Append(x0).ToArray();
        xi = lex ? xi : xi.OrderDescending().ToArray();
        var indeterminates = new Indeterminates<Xi>(order, xi.Select(x => new Xi(x)).ToArray());
        var xis = xi.Select(c => new Polynomial<K, Xi>(new Monom<Xi>(indeterminates, new(c), 1), scalar.One)).ToArray();
        var X0 = xis.First(e0 => !e0[new Monom<Xi>(indeterminates, new Xi(x0))].IsZero());
        var Xis = xis.Where(e0 => e0[new Monom<Xi>(indeterminates, new Xi(x0))].IsZero()).ToArray();
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
        var his = Polynomial(scalar, order, e);
        return his.Select(hi => new EPolynomial<K>(hi, hi.One)).ToArray();
    }

    public static EPolynomial<K>[] EPolynomial<K>(K scalar, (int n, string h) e)
        where K : struct, IFieldElt<K>, IElt<K>, IRingElt<K>
    {
        var his = Polynomial(scalar, MonomOrder.GrLex, e);
        return his.Select(hi => new EPolynomial<K>(hi, hi.One)).ToArray();
    }

    public static EPolynomial<K>[] PolynomialModI<K>(string[] vs, Polynomial<K, Xi> b0, params Polynomial<K, Xi>[] bs)
        where K : struct, IFieldElt<K>, IElt<K>, IRingElt<K>
    {
        var xis = vs.Select(expr => new Xi(expr)).ToArray();
        b0.Indeterminates.ExtendAppend(xis);
        var basis = new PolynomialBasis<K, Xi>(bs.Prepend(b0).ToArray());
        return xis.Select(xi => new EPolynomial<K>(b0.X(xi), basis)).ToArray();
    }

    public static EPolynomial<K>[] PolynomialModI<K>(KPoly<K> P, string x0, params string[] xis)
        where K : struct, IFieldElt<K>, IElt<K>, IRingElt<K>
    {
        var ind = Indeterminates(P.x.ToString());
        return PolynomialModI(xis.Prepend(x0).ToArray(), P.ToPolynomial(ind, ind.First()));
    }

    public static EPolynomial<K>[] PolynomialModI<K>((int n, string h) e, Polynomial<K, Xi> b0, params Polynomial<K, Xi>[] bs)
        where K : struct, IFieldElt<K>, IElt<K>, IRingElt<K>
    {
        var digits = $"{e.n - 1}".Length;
        var zeros = Enumerable.Repeat(0, digits).Glue();
        var fmt = $"{{0}}{{1,{digits}:{zeros}}}";
        var xis = e.n.Range().Select(i => string.Format(fmt, e.h, i)).ToArray();
        return PolynomialModI(xis, b0, bs);
    }

    public static EPolynomial<K>[] PolynomialModI<K>((int n, string h) e, Polynomial<K, Xi>[] bs)
        where K : struct, IFieldElt<K>, IElt<K>, IRingElt<K>
    {
        return PolynomialModI(e, bs[0], bs.Skip(1).ToArray());
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
        var g0 = s.One;
        var acc = s.Zero;
        for (int i = 0; i <= P.Degree; i++)
        {
            acc += P.Coefs[i].Poly.Substitute(a) * g0;
            g0 *= s;
        }

        return acc;
        return P.Coefs.Select((c, i) => c.Poly.Substitute(a) * s.Pow(i)).Aggregate(s.Zero, (acc, c) => acc + c);
    }

    public static KPoly<K> SubstituteChar<K>(this KPoly<K> f, char c)
        where K : struct, IFieldElt<K>, IElt<K>, IRingElt<K>
    {
        return new(c, f.KZero, f.Coefs);
    }
    
    public static KPoly<K> NewtonInverse<K>(KPoly<K> F, int N) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        if (F.Degree >= N)
            throw new($"F={F} and N={N}");

        if (F[0].IsZero())
            throw new($"F={F} is non invertible");

        if (N == 1)
            return F[0].Inv() * F.One;

        var mid = (N / 2) + (N % 2);
        var F0 = F.Div(F.X.Pow(mid)).rem;
        var G = NewtonInverse(F0, mid);
        return (G + (1 - F * G) * G).Div(F.X.Pow(N)).rem;
    }

    public static Vec<T> ToVec<T>(this IEnumerable<T> ts) where T : IElt<T>, IRingElt<T> => new(ts.ToArray());
    public static Vec<K> ToVec<K>(this KPoly<K> poly, int size) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        return new(poly.CoefsExtended(size));
    }
    public static Vec<K> ToVec<K>(this KPoly<K> poly) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
        => new(poly.Coefs.ToArray());
    public static Vec<Vec<K>> Transpose<K>(this Vec<Vec<K>> vec) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        if (vec.Select(e => e.Length).Distinct().Count() != 1)
            throw new("Cannot transpose");

        var n = vec[0].Length;
        return n.SeqLazy().Select(j => vec.Select(v => v[j]).ToVec()).ToVec();
    }
    public static MonomDisplay DisplayPolynomial { get; set; } = MonomDisplay.Default;
}
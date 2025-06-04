using System.Numerics;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Padic;
using GFelt = FastGoat.Structures.VecSpace.EPoly<FastGoat.UserGroup.Integers.ZnInt>;

namespace FastGoat.UserGroup.Polynoms;

public static partial class IntFactorisation
{
    static Polynomial<K, Xi>[] AllMonoms<K>(Polynomial<K, Xi>[] variables, int n)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var (x, y) = variables.Deconstruct();
        var sum = (n + 1).Range().Grid2D().Where(e => e.t1 + e.t2 <= n)
            .Select(e => x.Pow(e.t1) * y.Pow(e.t2)).Distinct().ToVec().Sum();
        return sum.Coefs.Select(e => new Polynomial<K, Xi>(e.Key, e.Value.One)).Order()
            .ToArray();
    }

    static Polynomial<K, Xi> RandPol<K>(Polynomial<K, Xi>[] monoms, int maxMonoms, int amplitude)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var m0 = monoms[0];
        var (y, x) = m0.Indeterminates.Deconstruct();
        while (true)
        {
            var f = IntExt.Rng.GetItems(monoms, IntExt.Rng.Next(2, maxMonoms + 1))
                .Select(e => IntExt.Rng.Next(-amplitude / 2, amplitude / 2 + 1) * e)
                .Aggregate(m0.Zero, (acc, e) => acc + e);
            var fy = f.ToKPoly(y);
            if (fy.Degree == 0 || f.NbIndeterminates != 2)
                continue;

            return f;
        }
    }

    public static IEnumerable<Polynomial<ZnInt, Xi>> RandZnIntPolynomials(int p, int maxMonomDegree = 4,
        int nbMonoms = 4, int nbFacts = 3, int nbRandPolys = 30)
    {
        var (Y, X) = Ring.Polynomial(ZnInt.ZnZero(p), MonomOrder.Lex, "Y", "X").Deconstruct();
        var monoms = AllMonoms([X, Y], maxMonomDegree);
        var pols = (10 * nbRandPolys).SeqLazy().Select(_ => RandPol(monoms, nbMonoms, p)).Distinct()
            .Where(IsFxyResultantZero).Take(nbRandPolys).ToArray();

        return (100 * nbRandPolys).SeqLazy().Select(_ => nbFacts.SeqLazy()
            .Select(_ => pols[IntExt.Rng.Next(pols.Length)])
            .Aggregate(X.One, (acc, e) => e * acc).Monic()
        );
    }

    public static IEnumerable<Polynomial<Rational, Xi>> RandRationalPolynomials(int amplitude, int maxMonomDegree = 4,
        int nbMonoms = 4, int nbFacts = 3, int nbRandPolys = 30)
    {
        var (Y, X) = Ring.Polynomial(Rational.KZero(), MonomOrder.Lex, "Y", "X").Deconstruct();
        var monoms = AllMonoms([X, Y], maxMonomDegree);
        var pols = (10 * nbRandPolys).SeqLazy().Select(_ => RandPol(monoms, nbMonoms, amplitude)).Distinct()
            .Where(IsFxyResultantZero).Take(nbRandPolys).ToArray();

        return (100 * nbRandPolys).SeqLazy().Select(_ => nbFacts.SeqLazy()
            .Select(_ => pols[IntExt.Rng.Next(pols.Length)]).Aggregate(X.One, (acc, e) => e * acc)
        ).Select(f => f.Primitive());
    }

    public static Polynomial<K, T> Swap<K, T>(Polynomial<K, T> F, T xi, T xj)
        where T : struct, IElt<T>
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var Fswap = new Polynomial<K, T>(F.Indeterminates, F.KZero,
            new(F.Coefs.ToDictionary(e => Monom<T>.Swap(e.Key, xi, xj), e => e.Value)));
        if (!xi.Equals(xj))
        {
            Console.WriteLine(new { Finit = F });
            Console.WriteLine(new { Fswap });
        }

        return Fswap;
    }

    public static KPoly<GFelt> ToGF(this KPoly<ZnInt> f, GFelt a)
    {
        return f.Select(c => c * a.One).ToKPoly();
    }

    public static KPoly<GFelt> ToGF(this KPoly<ZnInt> f) => f.ToGF(FG.FqX(f.P, 'a'));

    public static Polynomial<GFelt, T> ToGF<T>(this Polynomial<ZnInt, T> f, GFelt a) where T : struct, IElt<T>
    {
        var coefs = f.Coefs.ToDictionary(e => e.Key, e => e.Value * a.One);
        return new(f.Indeterminates, a.Zero, new(coefs));
    }

    public static Polynomial<GFelt, T> ToGF<T>(this Polynomial<ZnInt, T> f) where T : struct, IElt<T>
    {
        return f.ToGF(FG.FqX(f.P, 'a'));
    }

    public static Polynomial<GFelt, T> ToGF<T>(this Polynomial<GFelt, T> f, GFelt a) where T : struct, IElt<T>
    {
        var a0 = f.KZero;
        if (a0.P != a.P || a0.F.Degree != 1)
            throw new($"a0.P={a0.P} a0.F={a0.F}");
        var coefs = f.Coefs.ToDictionary(e => e.Key, e => e.Value[0] * a.One);
        return new(f.Indeterminates, a.Zero, new(coefs));
    }

    static KPoly<EPoly<K>> ToEPoly<K, T>(this Polynomial<K, T> f, Polynomial<K, T> P1, Polynomial<K, T> P2)
        where K : struct, IFieldElt<K>, IElt<K>, IRingElt<K>
        where T : struct, IElt<T>
    {
        var (x2, x1) = P1.Indeterminates.Deconstruct();
        var (Y, e2) = FG.EPolyXc(P1.ToKPoly(x1), $"{x1}"[0], $"{x2}"[0]);
        var decC2 = Ring.Decompose(P2, x2).Item1;
        var P2a = decC2.Select(e => e.Key.ToKPoly(x2).Substitute(Y) * (Y.One * e.Value.ToKPoly(x1).Substitute(e2)))
            .Aggregate((ai, aj) => ai + aj);

        if (P2.Equals(f))
            return P2a;

        var (decC, decY) = Ring.Decompose(f.Div(P2).rem, x2);
        var acc = Y.Zero;
        foreach (var (i, polY) in decY)
        {
            var polX = decC[polY].ToKPoly(x1).Substitute(e2);
            acc += Y.Pow(i) * polX;
        }

        return acc.Div(P2a).rem;
    }

    public static KPoly<KPoly<K>> ToKPoly<K, T>(this Polynomial<K, T> f, T x, T y)
        where K : struct, IFieldElt<K>, IElt<K>, IRingElt<K>
        where T : struct, IElt<T>
    {
        var X = FG.KPoly($"{x}"[0], f.KOne);
        var Y = FG.KPoly($"{y}"[0], X);
        var (decC, decY) = Ring.Decompose(f, y);
        return decY.Select(e => Y.Pow(e.Key) * (decC[e.Value].ToKPoly(x) * X.One * Y.One))
            .Aggregate((ai, aj) => ai + aj);
    }

    public static KPoly<FracPoly<K>> ToFracPoly<K, T>(this Polynomial<K, T> f, Polynomial<K, T> I, T x, T y)
        where K : struct, IFieldElt<K>, IElt<K>, IRingElt<K>
        where T : struct, IElt<T>
    {
        var X = FG.KFracPoly($"{x}"[0], f.KOne);
        var Y = FG.KPoly($"{y}"[0], X);
        var (decC, decY) = Ring.Decompose(f, y);
        return decY.Select(e => Y.Pow(e.Key) * (decC[e.Value].Div(I).rem.ToKPoly(x) * X.One * Y.One))
            .Aggregate((ai, aj) => ai + aj);
    }

    public static KPoly<FracPoly<K>> ToFracPoly<K, T>(this Polynomial<K, T> f, T x, T y)
        where K : struct, IFieldElt<K>, IElt<K>, IRingElt<K>
        where T : struct, IElt<T>
    {
        var X = FG.KFracPoly($"{x}"[0], f.KOne);
        var Y = FG.KPoly($"{y}"[0], X);
        var (decC, decY) = Ring.Decompose(f, y);
        return decY.Select(e => Y.Pow(e.Key) * (decC[e.Value].ToKPoly(x) * X.One * Y.One))
            .Aggregate((ai, aj) => ai + aj);
    }

    public static SPoly<K> ToSPoly<K, T>(this Polynomial<K, T> P, T xi, int ord)
        where T : struct, IElt<T>
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        return new(ord, P.ToKPoly(xi));
    }

    public static KPoly<SPoly<K>> ToSPoly<K, T>(this Polynomial<K, T> f, int ord, T x, T y)
        where K : struct, IFieldElt<K>, IElt<K>, IRingElt<K>
        where T : struct, IElt<T>
    {
        var X = new SPoly<K>(ord, FG.KPoly($"{x}"[0], f.KOne));
        var Y = FG.KPoly($"{y}"[0], X);
        var (decC, decY) = Ring.Decompose(f, y);
        return decY.Select(e => Y.Pow(e.Key) * (decC[e.Value].ToKPoly(x) * X.One * Y.One))
            .Aggregate((ai, aj) => ai + aj);
    }

    public static Polynomial<K, T> ToPolynomial<K, T>(this KPoly<EPoly<K>> f, Indeterminates<T> ind, T x, T y)
        where K : struct, IFieldElt<K>, IElt<K>, IRingElt<K>
        where T : struct, IElt<T>
    {
        var X = new Polynomial<K, T>(new Monom<T>(ind, x), f.KOne.KOne);
        var Y = new Polynomial<K, T>(new Monom<T>(ind, y), f.KOne.KOne);
        return f.Select((c, i) => Y.Pow(i) * c.Poly.Substitute(X)).Aggregate(X.Zero, (acc, xi) => acc + xi);
    }

    public static Polynomial<K, T> ToPolynomial<K, T>(this KPoly<SPoly<K>> f, Indeterminates<T> ind, T x, T y)
        where K : struct, IFieldElt<K>, IElt<K>, IRingElt<K>
        where T : struct, IElt<T>
    {
        var X = new Polynomial<K, T>(new Monom<T>(ind, x), f.KOne.KOne);
        var Y = new Polynomial<K, T>(new Monom<T>(ind, y), f.KOne.KOne);
        return f.Select((c, i) => Y.Pow(i) * c.Poly.Substitute(X)).Aggregate(X.Zero, (acc, xi) => acc + xi);
    }

    public static KPoly<FracPoly<K>> ConvertToJaggedArrayPoly<K, T>(Polynomial<K, T> F, T x, T y)
        where T : struct, IElt<T>
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var Y = FG.KPoly($"{y}"[0], FG.KFracPoly($"{x}"[0], F.KOne));
        var X = Y.KOne.X * Y.One;
        var Fcoefs = Ring.Decompose(F, y).Item1;
        return Fcoefs.Select(e => e.Key.ToKPoly(y).Substitute(Y) * e.Value.ToKPoly(x).Substitute(X))
            .Aggregate((acc, e) => e + acc);
    }

    public static Polynomial<K, T> ConvertToDictionaryPoly<K, T>(KPoly<FracPoly<K>> F, Polynomial<K, T> X,
        Polynomial<K, T> Y)
        where T : struct, IElt<T>
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var lcm = F.Coefs.Select(f => f.Denom).Aggregate(Ring.Lcm);
        var F0 = new FracPoly<K>(lcm) * F;
        if (F0.Coefs.Any(c => !(c.Denom - 1).IsZero()))
            throw new();
        return F0.Coefs.Select((c0, i) => c0.Num.Substitute(X) * Y.Pow(i)).Aggregate((acc, e) => e + acc);
    }

    static Polynomial<Rational, Xi> ZPoly2QPoly(Polynomial<ZnBInt, Xi> f, ZnBInt c)
    {
        var coefs1 = f.Coefs.ToDictionary(e => e.Key, e => e.Value * c)
            .ToDictionary(
                e => e.Key,
                e => new Rational(e.Value.K * 2 <= e.Value.Mod ? e.Value.K : e.Value.K - e.Value.Mod))
            .Where(e => !e.Value.IsZero())
            .ToDictionary(e => e.Key, e => e.Value);
        return new Polynomial<Rational, Xi>(f.Indeterminates, Rational.KZero(), new(coefs1));
    }

    public static Polynomial<Rational, Xi> Primitive(this Polynomial<Rational, Xi> f)
    {
        if (f.IsZero())
            return f;

        var arrGcd = f.Coefs.Values.Where(e => !e.IsZero()).Select(e => e.Absolute.Num).Distinct().Order().ToArray();
        var arrLcm = f.Coefs.Values.Select(e => e.Absolute.Denom).Distinct().Order().ToArray();
        return f * new Rational(f.LeadingDetails.lc.Sign * IntExt.LcmBigInt(arrLcm), IntExt.GcdBigInt(arrGcd));
    }

    public static Polynomial<ZnInt, Xi> ToZnInt(this Polynomial<Rational, Xi> f, int p)
    {
        var z = new ZnInt(p, 0);
        var coefs = f.Coefs.Select(e => (e.Key, e.Value.ToZnInt(p)))
            .Where(e => !e.Item2.IsZero()).ToDictionary(e => e.Key, e => e.Item2);
        return new(f.Indeterminates, z, new(coefs));
    }

    public static (bool Lt, K gk, Polynomial<K, T> F) RewritingPolynomialLeadingTerm<K, T>(Polynomial<K, T> F,
        IEnumerable<K> seq)
        where T : struct, IElt<T>
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var ((y, Y), (x, X)) = F.AllIndeterminatesAndVariables.Deconstruct();
        if (F.NbIndeterminates == 1 || F.CoefMax(y).Degree == 0)
            return (true, F.KZero, F);

        foreach (var gk in seq)
        {
            var subs1 = new List<(T, Polynomial<K, T>)>()
            {
                (x, X + gk * Y),
                (y, Y)
            };

            var F1 = F.Substitute(subs1);
            if (F1.CoefMax(y).Degree != 0)
                continue;

            return (true, gk, F1);
        }

        if (Logger.Level != LogLevel.Off)
            Console.WriteLine($"Warning LT problem F={F}");

        return (false, F.KZero, F);
    }

    public static (bool Lt, K gk, Polynomial<K, T> F) RewritingPolynomialLeadingTerm<K, T>(Polynomial<K, T> F, K g,
        int q)
        where T : struct, IElt<T>
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        return RewritingPolynomialLeadingTerm(F, (q - 1).SeqLazy(1).Select(k => g.Pow(k)));
    }

    public static (bool Lt, Rational gk, Polynomial<Rational, Xi> F)
        RewritingPolynomialLeadingTerm(Polynomial<Rational, Xi> F)
    {
        var ((x, X), (y, Y)) = F.IndeterminatesAndVariables.Deconstruct();
        var degY = F.DegreeOf(y);
        var seq = (2 * degY + 1).Range(-degY).Select(i => new Rational(i))
            .OrderBy(e => e.Absolute).ThenDescending().ToArray();
        return RewritingPolynomialLeadingTerm(F, seq);
    }

    public static bool IsFxyResultantZero<K, T>(Polynomial<K, T> F)
        where T : struct, IElt<T>
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var (y, x) = F.Indeterminates.Deconstruct();
        if (F.NbIndeterminates < 2)
            return false;

        var F0y = F.Substitute(F.Zero, x).ToKPoly(y);
        var disc = Ring.Discriminant(F0y);
        return F0y.Degree > 0 && !disc.IsZero();
    }

    public static (bool resNotNul, K gk, Polynomial<K, T> F) RewritingPolynomialResultantZero<K, T>(Polynomial<K, T> F,
        IEnumerable<K> seq)
        where T : struct, IElt<T>
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var ((y, Y), (x, X)) = F.AllIndeterminatesAndVariables.Deconstruct();
        if (F.NbIndeterminates == 1)
            return (true, F.KZero, F);

        foreach (var gk in seq)
        {
            var F1 = F.Substitute(X + gk, x);
            if (!IsFxyResultantZero(F1))
                continue;

            return (true, gk, F1);
        }

        if (Logger.Level != LogLevel.Off)
            Console.WriteLine($"Warning Non separable polynomial");

        return (false, F.KZero, F);
    }

    public static (bool resNotNul, K gk, Polynomial<K, T> F)
        RewritingPolynomialResultantZero<K, T>(Polynomial<K, T> F, K g, int q)
        where T : struct, IElt<T>
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        return RewritingPolynomialResultantZero(F, (q - 1).SeqLazy(1).Select(k => g.Pow(k)).Prepend(g.Zero));
    }

    public static (bool resNotNul, Rational gk, Polynomial<Rational, Xi> F)
        RewritingPolynomialResultantZero(Polynomial<Rational, Xi> F)
    {
        var ((x, X), (y, Y)) = F.IndeterminatesAndVariables.Deconstruct();
        var degY = F.DegreeOf(y);
        var seq = (2 * degY + 1).Range(-degY).Select(i => new Rational(i))
            .OrderBy(e => e.Absolute).ThenDescending().ToArray();
        return RewritingPolynomialResultantZero(F, seq);
    }

    public static (Modulus po, KPoly<ZnBInt>[] firr, KPoly<ZnBInt>[] firr1) FirrRational(KPoly<Rational> F0y,
        Rational lc, int p,
        int o)
    {
        var k = NumberTheory.PrimitiveRootMod(p);
        var a0 = ZnBInt.ZnZero(p) + k;
        var po = a0.Details;
        var f0y = QPoly2ZnInt(F0y, po);
        var firr = CantorZassenhausAECF(f0y.Monic, a0, p).ToArray();

        var all = new List<KPoly<ZnBInt>>(firr);
        while (po.O < o && all.Count > 1)
        {
            po *= 2;
            var tmp = new List<KPoly<ZnBInt>>();
            var fa = QPoly2ZnInt(F0y, po) * lc.ToZnBInt(po).Inv();
            foreach (var g in all)
            {
                var gi = ZPoly2ZnInt(g, po);
                var y = FG.EPoly(gi);
                var dgi = gi.Derivative.Substitute(y);
                var fi = fa.Substitute(y);
                var dfi = fa.Derivative.Substitute(y);
                var ri = (dgi * fi * dfi.Inv()).Poly;
                tmp.Add(gi + ri);
            }

            all.Clear();
            all = tmp.ToList();
        }

        return (po, firr, all.ToArray());
    }

    static Polynomial<K, Xi>[] HenselLiftingStep<K>(Polynomial<K, Xi> F, Polynomial<K, Xi>[] fi,
        Polynomial<K, Xi> I) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var xis = F.ExtractAllIndeterminates;
        if (xis.Length != 2)
            throw new();

        var (x, y) = xis.Deconstruct();
        var P0 = F;
        var P1 = F.D(y);
        var tmp = new List<Polynomial<K, Xi>>();
        foreach (var f in fi)
        {
            var ord = I.Degree;
            try
            {
                var F01a = P0.Div(f).rem.Div(I).rem.ToSPoly(ord, x, y);
                var F11a = P1.Div(f).rem.Div(I).rem.ToSPoly(ord, x, y);
                var f0a = f.Div(I).rem.ToSPoly(ord, x, y);
                var df0a = f0a.Derivative;
                var F11ai = Ring.FastInv(F11a, f0a);
                if (!(F11a * F11ai).Div(f0a).rem.IsOne())
                    throw new("Invalid inverse #1. F11a * F11ai must equal 1");

                var R1 = (df0a * F01a * F11ai).Div(f0a).rem;
                var fr = (f0a + R1).Monic.ToPolynomial(f.Indeterminates, x, y).Div(I).rem;
                tmp.Add(fr);
            }
            catch (Exception e)
            {
                if (Logger.Level != LogLevel.Off)
                    Console.WriteLine($"    {e.Message}");

                var f0a = f.Div(I).rem.ToFracPoly(x, y);
                var F01a = P0.Div(f).rem.Div(I).rem.ToFracPoly(x, y);
                var F11a = P1.Div(f).rem.Div(I).rem.ToFracPoly(x, y);
                var df0a = f0a.Derivative;
                var F11ai = Ring.FastInv(F11a, f0a);
                if (!(F11a * F11ai).Div(f0a).rem.IsOne())
                    throw new("Invalid inverse #2. F11a * F11ai must equal 1");

                var R1 = (df0a * F01a * F11ai).Div(f0a).rem
                    .Select(c => new FracPoly<K>(c.Num * Ring.NewtonInverse(c.Denom, c.Denom.Degree + 1)))
                    .ToKPoly($"{y}"[0]);
                var fr = ConvertToDictionaryPoly((f0a + R1).Monic, F.X(x), F.X(y)).Div(I).rem;
                tmp.Add(fr);
            }
        }

        return tmp.Order().ToArray();
    }

    public static Polynomial<K, Xi>[] HenselLifting<K>(Polynomial<K, Xi> F, Polynomial<K, Xi>[] firr, int o)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var (y, x) = F.Indeterminates.Deconstruct();
        var all = firr.ToArray();

        var I = F.X(x);
        while (I.Degree <= o)
        {
            I *= I;
            all = HenselLiftingStep(F, all, I);
        }

        return all;
    }

    public static Polynomial<K, Xi>[] HenselLifting<K>(Polynomial<K, Xi> F, Polynomial<K, Xi>[] firr)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var (_, x) = F.Indeterminates.Deconstruct();
        var o = F.DegreeOf(x) + 1;
        return HenselLifting(F, firr, o);
    }

    public static Polynomial<K, Xi>[] Recombinaison<K>(Polynomial<K, Xi> F, Polynomial<K, Xi>[] fi)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var (y, x) = F.Indeterminates.Deconstruct();
        var o = F.DegreeOf(x) + 1;
        var xo = F.X(x).Pow(o);

        var facts = new List<Polynomial<K, Xi>>();
        var rem = new HashSet<Polynomial<K, Xi>>(fi);
        var nbCombs = 1;
        while (rem.Count != 0)
        {
            var sz = rem.Count;
            foreach (var comb in rem.AllCombinationsFromM(nbCombs).ToArray())
            {
                var fact = comb.Aggregate(F.One, (acc, li) => acc * li).Div(xo).rem;
                if (F.Div(fact).rem.IsZero())
                {
                    facts.Add(fact);
                    rem.ExceptWith(comb);
                    nbCombs = comb.Length;
                    break;
                }
            }

            if (rem.Count == sz)
                break;
        }

        var prod = facts.Aggregate(F.One, (acc, fj) => acc * fj);
        if (prod.Degree < F.Degree)
        {
            var (quo1, rem1) = F.Div(prod);
            if (rem1.IsZero())
                facts.Add(quo1);
            else
                throw new($"Recombinaison F={F} Prod={prod} (quo,rem)={(quo1, rem1)}");
        }

        return facts.Order().ToArray();
    }

    public static Polynomial<Rational, Xi>[] Recombinaison(Polynomial<Rational, Xi> F, Polynomial<ZnBInt, Xi>[] fi,
        ZnBInt c)
    {
        var factsQ = new List<Polynomial<Rational, Xi>>();
        var rem = new HashSet<Polynomial<ZnBInt, Xi>>(fi);
        var nbCombs = 1;

        while (rem.Count != 0)
        {
            var sz = rem.Count;
            foreach (var comb in rem.AllCombinationsFromM(nbCombs).ToArray())
            {
                var fact = comb.Aggregate((acc, e) => e * acc);
                var factQ1 = ZPoly2QPoly(fact, c.One);
                if (F.Div(factQ1).rem.IsZero())
                {
                    factsQ.Add(factQ1.Primitive());
                    rem.ExceptWith(comb);
                    nbCombs = comb.Length;
                    break;
                }

                var factQc = ZPoly2QPoly(fact, c);
                if (F.Div(factQc).rem.IsZero())
                {
                    factsQ.Add(factQc.Primitive());
                    rem.ExceptWith(comb);
                    nbCombs = comb.Length;
                    break;
                }
            }

            if (rem.Count == sz)
                break;
        }

        return factsQ.Order().ToArray();
    }

    public static Polynomial<Rational, Xi>[] FactorsFxyStep(Polynomial<Rational, Xi> F)
    {
        if (F.Coefs.Any(e => !e.Value.IsInteger()))
            throw new($"F = {F}");

        var (y, x) = F.Indeterminates.Deconstruct();
        var F0y = F.Substitute(F.Zero, x).ToKPoly(y);
        var lc0 = F.LeadingDetails.lc;
        var pMin = 2 + F.DegreeOf(y) * (2 * F.DegreeOf(x) - 1);
        if (Logger.Level != LogLevel.Off)
        {
            Console.WriteLine($"FactorsFxyStep : {F}");
            Console.WriteLine($"F(0,{y,2})        : {F0y.ToPolynomial(F.Indeterminates, y)}");
        }

        foreach (var (p, o) in PSigma(F0y, 600).Where(e => e.p >= pMin))
        {
            try
            {
                var (po, _, firr1) = FirrRational(F0y, lc0, p, o);
                if (Logger.Level != LogLevel.Off)
                    firr1.Select(f => f.ToPolynomial(F.Indeterminates, y))
                        .Println($"Factors F(0,{y}) in Z/({p}^{o})Z[{x},{y}]");

                var firrPAdic = firr1.Select(f => f.Monic.ToPolynomial(F.Indeterminates, y)).ToArray();
                var coefs = F.Coefs.ToDictionary(e => e.Key, e => e.Value.ToZnBInt(po))
                    .Where(e => !e.Value.IsZero())
                    .ToDictionary(e => e.Key, e => e.Value);
                var Fp = new Polynomial<ZnBInt, Xi>(F.Indeterminates, ZnBInt.ZnZero(p, po.O), new(coefs));
                var lc1 = Fp.LeadingDetails.lc;
                Fp *= lc1.Inv();

                var lifts = HenselLifting(Fp, firrPAdic);
                var factsFp = Recombinaison(Fp, lifts);

                var factsQ = Recombinaison(F, factsFp, lc1);
                if (Logger.Level != LogLevel.Off)
                {
                    lifts.Println("HenselLifting");
                    factsFp.Println($"Recombinaison in Z/({p}^{o})Z[{x},{y}]");
                }

                if (factsQ.Length == 0)
                    throw new("Message:Rational Recombinaison"); // TODO

                if (Logger.Level != LogLevel.Off)
                    factsQ.Println($"Recombinaison in Q[{x},{y}]");

                return factsQ.Select(f => f.Primitive()).Order().ToArray();
            }
            catch (Exception e)
            {
                if (Logger.Level != LogLevel.Off)
                    Console.WriteLine($"########### P = {p,-5} and O = {o} wont work. {e.Message}");
            }
        }

        return [F];
    }

    public static Polynomial<Rational, Xi>[] FactorsFxy(Polynomial<Rational, Xi> F, bool rewrite = false)
    {
        var (i0, i1, F1, F2) = (F.KZero, F.KZero, F, F);
        var ((y, Y), (x, X)) = F2.AllIndeterminatesAndVariables.Deconstruct();
        if (rewrite)
        {
            var degY = F.DegreeOf(y);
            var seq = (2 * degY + 1).Range(-degY).Select(i => new Rational(i))
                .OrderBy(i => i.Absolute).ThenDescending().ToArray();
            (_, i0, F1) = RewritingPolynomialLeadingTerm(F, seq);
            F1 = F1.Primitive();
            (_, i1, F2) = RewritingPolynomialResultantZero(F1, seq);
            F2 = F2.Primitive();
        }

        if (Logger.Level != LogLevel.Off)
        {
            Console.WriteLine($"FactorsFxy F   = {F}");
            var Sub = X.Substitute((X + i0 * Y, x), (X + i1, x));
            if (!X.Equals(Sub))
                Console.WriteLine($"Substitute {x} <-- {Sub}");

            if (!F.Equals(F2))
                Console.WriteLine($"Rewrited       = {F2}");
        }

        var facts = new List<Polynomial<Rational, Xi>>();
        var Frem = F2;
        while (Frem.Degree > 0)
        {
            var factsQ = FactorsFxyStep(Frem);
            facts.AddRange(factsQ);
            var prod = factsQ.Aggregate(F2.One, (acc, f) => f * acc);
            Frem = Primitive(Frem.Div(prod).quo);
        }

        if (rewrite)
        {
            if (i1 != 0)
                facts = facts.Select(f => Primitive(f.Substitute(X - i1, X))).Order().ToList();
            if (i0 != 0)
                facts = facts.Select(f => Primitive(f.Substitute(X - i0 * Y, X))).Order().ToList();
        }

        var lt = F / facts.Aggregate(F.One, (acc, f) => f * acc);
        if (!(lt - 1).IsZero())
            facts.Add(lt);

        var factsFinal = facts.Order().ToArray();
        var check = F.Equals(factsFinal.Aggregate((acc, e) => e * acc));
        if (!check)
            throw new("Wrong factorisation");

        return factsFinal;
    }

    public static (Polynomial<K, T> F, Polynomial<K, T> ctx, Polynomial<K, T> cty) CT<K, T>(Polynomial<K, T> F)
        where T : struct, IElt<T>
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var (y, x) = F.Indeterminates.Deconstruct();
        var decY = Ring.Decompose(F, y).Item1;
        var gcdY = decY.Keys.Order().FirstOrDefault(e => !decY[e].IsZero(), F.One);
        var seqX = decY.Values.Where(e => !e.IsZero()).Select(e => e.ToKPoly(x)).ToArray();
        var gcdX = seqX.Length == 0 ? F.One : Ring.FastGCD(seqX).ToPolynomial(F.Indeterminates, x).Monic();
        var F1 = gcdX * gcdY;
        var (F2, rem) = F.Div(F1);
        if (rem.IsZero())
            return (F2, gcdX, gcdY);

        Console.WriteLine(new { F });
        decY.Println();
        Console.WriteLine(new { gcdX });
        Console.WriteLine(new { gcdY });
        Console.WriteLine(new { F1 });
        Console.WriteLine(new { F2 });
        Console.WriteLine(new { rem });
        throw new();
    }

    public static (Polynomial<K, T> g, int m) Deflate<K, T>(Polynomial<K, T> F)
        where T : struct, IElt<T>
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var p = F.P;
        if (p == 0)
            return (F, 1);

        if (F.Coefs.Keys.All(mn => mn.ToTuples().All(e => e.Item2 % p == 0)))
        {
            var coefs = F.Coefs.Keys.Where(mn => mn.Degree > 0).ToDictionary(
                mn => mn,
                mn => mn.ToTuples().Where(e => e.Item2 != 0)
                    .Select(e => new Monom<T>(F.Indeterminates, e.Item1, e.Item2 / p)).Aggregate((mi, mj) => mi.Mul(mj))
            );

            var coefs1 = coefs.ToDictionary(e => e.Value, e => F.Coefs[e.Key]);
            var F1 = new Polynomial<K, T>(F.Indeterminates, F.KZero, new(coefs1)) + F.ConstTerm;
            return (F1, p);
        }

        return (F, 1);
    }

    static (Polynomial<K, T> g, int m)[] SFFStep<K, T>(Polynomial<K, T> F, T x, T y)
        where T : struct, IElt<T>
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var (X, Y) = (F.X(x), F.X(y));
        var sff = MusserSFF(ConvertToJaggedArrayPoly(F, x, y)).Where(e => e.g.Degree > 0)
            .Select(e => (g: ConvertToDictionaryPoly(e.g, X, Y), i: e.i * e.q)).OrderBy(e => e.i)
            .ToList();

        var prod = sff.Aggregate(F.One, (acc, sj) => acc * sj.g.Pow(sj.i));
        var (quo, rem) = F.Div(prod);
        if (!rem.IsZero())
            throw new();
        if (quo.Degree > 0)
        {
            if (quo.P > 0)
            {
                var dquo = Deflate(quo);
                sff.Add((dquo.g, dquo.m));
            }
            else
                sff.Add((quo, 1));
        }

        return sff.OrderBy(e => e.Item1).ToArray();
    }

    static (Polynomial<K, T> g, int m)[] SFFrec<K, T>(Polynomial<K, T> F)
        where T : struct, IElt<T>
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var (y, x) = F.Indeterminates.Deconstruct();
        var sffXY = SFFStep(F, x, y);
        var sffYX = SFFStep(F, y, x);
        if (sffYX.Sum(e => e.m) > sffXY.Sum(e => e.m))
            return sffYX;

        return sffXY;
    }

    public static (Polynomial<K, T> g, int m)[] SFF<K, T>(Polynomial<K, T> F, bool rec = false)
        where T : struct, IElt<T>
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        if (!rec)
        {
            var (y, x) = F.Indeterminates.Deconstruct();
            return SFFStep(F, x, y);
        }

        var seq = new List<(Polynomial<K, T> g, int m)>() { (F, 1) };
        int sz = 0;
        while (seq.Count != sz)
        {
            sz = seq.Count;
            seq = seq.SelectMany(e => SFFrec(e.g).Select(f => (f.g, f.m * e.m))).ToList();
        }

        return seq.ToArray();
    }

    static Polynomial<K, Xi>[] FactorStepFxyFq<K>(Polynomial<K, Xi> F, K g, BigInteger q)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var (Y, X) = F.AllVariables.Deconstruct();
        var (y, x) = F.Indeterminates.Deconstruct();

        if (F.DegreeOf(y) == 0)
        {
            var Fx = F.ToKPoly(x);
            return FirrFsepCantorZassenhausAECF(Fx, g, q)
                .SelectMany(e => Enumerable.Repeat(e.g, e.m * e.q))
                .Select(e => e.ToPolynomial(F.Indeterminates, x))
                .ToArray();
        }

        var Fy = F.Substitute(F.Zero, X).ToKPoly(Y);
        var firr = FirrFsepCantorZassenhausAECF(Fy, g, q)
            .Select(f => (f.g, m: f.q * f.m))
            .SelectMany(f => Enumerable.Repeat(f.g.ToPolynomial(F.Indeterminates, Y.ExtractIndeterminate), f.m))
            .ToArray();

        if (Logger.Level != LogLevel.Off)
            firr.Println($"Fy = {Fy}");

        if (firr.Length == 1)
            return [F];

        var lifts = HenselLifting(F, firr);
        if (Logger.Level != LogLevel.Off)
            lifts.Println("Lifts Fy");

        var facts = Recombinaison(F, lifts);
        if (Logger.Level != LogLevel.Off)
            facts.Println("Recombs F");

        return facts;
    }

    static Polynomial<ZnInt, Xi>[] FactorStepFxyFp(Polynomial<ZnInt, Xi> F)
    {
        if (Logger.Level != LogLevel.Off)
            Console.WriteLine($"FactorStepFxyFp : {F}");

        if (F.NbIndeterminates < 2)
            return [F];

        var g0 = NumberTheory.PrimitiveRootMod(F.P) * F.KOne;
        var p = F.P;
        var (lt0, i10, P10) = RewritingPolynomialLeadingTerm(F, g0, p);
        var (resNotNull0, i20, F20) = RewritingPolynomialResultantZero(P10, g0, p);
        if (lt0 && resNotNull0)
        {
            var (y, x) = F20.Indeterminates.Deconstruct();
            var (Y, X) = F20.AllVariables.Deconstruct();
            var facts = FactorStepFxyFq(F20, g0, p);
            var facts3 = facts.Select(f => f.Substitute((X - i20, x), (X - i10 * Y, x)).Monic()).ToList();
            if (Logger.Level != LogLevel.Off)
                facts3.Println($"StepZnInt F = {F}");

            return facts3.ToArray();
        }
        else
        {
            if (Logger.Level != LogLevel.Off)
                Console.WriteLine(new { q = p, lt = lt0, resNotNull = resNotNull0 });

            var g1 = FG.FqX(p, 'a');
            var (lt, i1, P1) = (false, g1.Zero, F.ToGF(g1));
            var (resNotNull, i2, F2) = (false, g1.Zero, F.ToGF(g1));
            var q = p;
            while (!lt || !resNotNull)
            {
                q *= F.P;
                g1 = FG.FqX(q, 'a');
                (lt, i1, P1) = RewritingPolynomialLeadingTerm(F.ToGF(g1), g1, q);
                (resNotNull, i2, F2) = RewritingPolynomialResultantZero(P1, g1, q);

                if (Logger.Level != LogLevel.Off)
                    Console.WriteLine(new { q, lt, resNotNull, F2 });

                if (q > 500)
                    throw new();
            }

            var (y, x) = F2.Indeterminates.Deconstruct();
            var (Y, X) = F2.AllVariables.Deconstruct();
            var facts = FactorStepFxyFq(F2, g1, q);

            var facts3 = facts.Select(f => f.Substitute((X - i2, x), (X - i1 * Y, x)).Monic()).ToList();
            var factFq = facts3.Where(e => e.Coefs.Any(c => c.Value.Degree > 0))
                .Aggregate(F2.One, (acc, aj) => acc * aj);
            if (factFq.Degree > 0)
            {
                if (Logger.Level != LogLevel.Off)
                    facts3.Println($"Problems mod={F2.KOne.F}");

                facts3 = facts3.Where(e => e.Coefs.All(c => c.Value.Degree == 0)).Append(factFq).ToList();
            }

            var o = F.KOne;
            var facts4 = facts3.Select(f =>
                    new Polynomial<ZnInt, Xi>(F.Indeterminates, F.KZero,
                        new(f.Coefs.ToDictionary(e => e.Key, e => e.Value[0] * o))))
                .ToArray();

            if (Logger.Level != LogLevel.Off)
                facts4.Println($"StepGFelt F = {F}");

            return facts4;
        }
    }

    public static (Polynomial<ZnInt, Xi> g, int m)[] FactorFxyFp(Polynomial<ZnInt, Xi> F)
    {
        var (F1, ctx, cty) = CT(F);
        var (y, x) = F.Indeterminates.Deconstruct();
        var sff = SFF(F1, rec: true).ToList();
        if (ctx.Degree > 0)
        {
            var Fx = ctx.ToKPoly(x);
            var g = Fx.KOne * NumberTheory.PrimitiveRootMod(Fx.P);
            sff.AddRange(FirrFsepCantorZassenhausAECF(Fx, g, Fx.P)
                .Select(e => (e.g.ToPolynomial(ctx.Indeterminates, x), m: e.m * e.q)));
        }

        if (cty.Degree > 0)
            sff.Add((F.X(y), cty.Degree));

        if (Logger.Level != LogLevel.Off)
        {
            Console.WriteLine($"FactorFxyFp : {F}");
            sff.Println("SFF");
        }

        return sff.SelectMany(f => FactorStepFxyFp(f.g).Select(g => (g.Monic(), f.m)))
            .Select(e => (e.m, d: Deflate(e.Item1)))
            .Select(e => (e.d.g, e.m * e.d.m))
            .GroupBy(e => e.g).Select(e => (g: e.Key, m: e.Sum(f => f.Item2)))
            .ToArray();
    }
}
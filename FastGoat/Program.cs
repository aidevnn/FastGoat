using System.CodeDom;
using System.Collections;
using System.Diagnostics;
using System.Numerics;
using FastGoat.Commons;
using FastGoat.Examples;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.Naming;
using FastGoat.Structures.Subgroups;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Characters;
using FastGoat.UserGroup.DatabaseSmallGroups;
using FastGoat.UserGroup.Floats;
using FastGoat.UserGroup.GModuleN;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Matrix;
using FastGoat.UserGroup.GModuleN;
using FastGoat.UserGroup.Perms;
using FastGoat.UserGroup.Polynoms;
using FastGoat.UserGroup.Words;
using static FastGoat.Commons.IntExt;
using static FastGoat.Commons.EnumerableExt;
using FastGoat.UserGroup.Padic;
using FastGoat.UserGroup.Words.Tools;
using GroupRegX = System.Text.RegularExpressions;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");

KPoly<K> RecNewtonInverse<K>(KPoly<K> F, int N) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
{
    if (F.Degree >= N)
        throw new($"F={F} and {N}");

    if (N == 1)
        return F[0].Inv() * F.One;

    var mid = (N / 2) + (N % 2);
    var F0 = F.Div(F.X.Pow(mid)).rem;
    var G = RecNewtonInverse(F0, mid);
    var G0 = (G + (1 - F * G) * G).Div(F.X.Pow(N)).rem;
    return G0;
}

Polynomial<ZnInt, Xi>[] Firr(Polynomial<ZnInt, Xi> F)
{
    if (F.P == 0)
        throw new();

    var p = F.P;
    var (x, y) = F.ExtractAllIndeterminates.Deconstruct();
    var P0y = F.Substitute(F.Zero, x);
    var _P0y = P0y.ToKPoly(y);
    var k = SolveAll_k_pow_m_equal_one_mod_n_strict(p, p - 1).First();
    var a0 = k * _P0y.KOne;
    var firr = IntFactorisation.BerlekampProbabilisticAECF(_P0y, a0).Order().ToArray();
    return firr.Select(fi => fi.Substitute(F.X(y))).ToArray();
}

Polynomial<ZnInt, Xi>[] HenselLiftingStep(Polynomial<ZnInt, Xi> F, Polynomial<ZnInt, Xi>[] fi, Polynomial<ZnInt, Xi> I,
    int o)
{
    var xis = F.ExtractAllIndeterminates;
    if (xis.Length != 2)
        throw new();

    var (x, y) = xis.Deconstruct();
    var P0 = F;
    var P1 = F.D(y);
    var tmp = new List<Polynomial<ZnInt, Xi>>();
    foreach (var f in fi)
    {
        var df = f.D(y).Div(I).rem;
        var F0 = P0.Div(f).rem.Div(I).rem;
        var F1 = P1.Div(f).rem.Div(I).rem;
        var _F0 = F0.ToKPoly(x);
        var _F1 = F1.ToKPoly(x);
        var _x = _F0.X;
        var div = (RecNewtonInverse(_F1, int.Max(_F1.Degree, o) + 1) * _F0).Div(_x.Pow(o)).rem;
        var R1 = (df * div.ToPolynomial(F.Indeterminates, x)).Div(f).rem.Div(I).rem;
        var fr = f + R1;
        var c = fr[new(fr.Indeterminates, y)];
        tmp.Add(fr * c.Inv());
    }

    return tmp.Order().ToArray();
}

Polynomial<ZnInt, Xi>[] HenselLifting(Polynomial<ZnInt, Xi> F, Polynomial<ZnInt, Xi>[] firr)
{
    var (x, y) = F.ExtractAllIndeterminates.Deconstruct();
    var all = firr.ToArray();

    var o = F.DegreeOf(x) + 1;
    var I = F.X(x);
    while (I.Degree < o && all.Length > 1)
    {
        I = I.Pow(2);
        all = HenselLiftingStep(F, all, I, o);
    }

    return all;
}

Polynomial<ZnInt, Xi>[] Recombinaison(Polynomial<ZnInt, Xi> F, Polynomial<ZnInt, Xi>[] fi)
{
    var (x, y) = F.ExtractAllIndeterminates.Deconstruct();
    var o = F.DegreeOf(x) + 1;
    var xo = F.X(x).Pow(o);

    var facts = new List<Polynomial<ZnInt, Xi>>();
    var rem = new HashSet<Polynomial<ZnInt, Xi>>(fi);
    while (rem.Count != 0)
    {
        var sz = rem.Count;
        var combs = rem.AllCombinations();
        foreach (var l in combs.Where(l => l.Length != 0))
        {
            var fact = l.Aggregate(F.One, (acc, li) => acc * li).Div(xo).rem;
            if (F.Div(fact).rem.IsZero())
            {
                facts.Add(fact);
                rem.ExceptWith(l);
                break;
            }
        }

        if (rem.Count == sz)
            throw new();
    }

    return facts.Order().ToArray();
}

void FactorsFxy(Polynomial<Rational, Xi> F)
{
    if (F.Coefs.Any(e => !e.Value.IsInteger()))
        throw new();

    var (x, y) = F.ExtractAllIndeterminates.Deconstruct();
    var disc = Ring.Discriminant(F.Substitute(F.Zero, x), y);
    var decomp = PrimesDecompositionBigInt(disc.ConstTerm.Absolute.Num).Distinct();
    Console.WriteLine($"F({x},{y} = {F}");
    foreach (var p in Primes10000.Except(decomp).Where(p => p < 500))
    {
        var coefs = F.Coefs.ToDictionary(e => e.Key, e => new ZnInt(p, (int)e.Value.Num))
            .Where(e => !e.Value.IsZero())
            .ToDictionary(e => e.Key, e => e.Value);
        var Fp = new Polynomial<ZnInt, Xi>(F.Indeterminates, ZnInt.ZpZero(p), new(coefs));
        try
        {
            Console.WriteLine(new { Fp, p });
            var firr = Firr(Fp);
            var lifts = HenselLifting(Fp, firr);
            var facts = Recombinaison(Fp, lifts);

            var P0X2 = firr.Aggregate(Fp.One, (acc, c) => acc * c);
            var P0 = facts.Aggregate(Fp.One, (acc, c) => acc * c);
            if (!Fp.Equals(P0))
                throw new();

            var facts1 = facts.Select(f =>
            {
                var coefs1 = f.Coefs.ToDictionary(e => e.Key,
                        e => new Rational(e.Value.K * 2 <= p ? e.Value.K : e.Value.K - p))
                    .Where(e => !e.Value.IsZero())
                    .ToDictionary(e => e.Key, e => e.Value);
                return new Polynomial<Rational, Xi>(F.Indeterminates, Rational.KZero(), new(coefs1));
            }).ToArray();

            var P1 = facts1.Aggregate(F.One, (acc, c) => acc * c);
            if (!F.Equals(P1))
                throw new();

            Console.WriteLine($"P(X1,X2) = {Fp}");
            firr.Println($"P(0,X2) = {P0X2}");
            lifts.Println("Hensel Lifting");
            facts.Println($"Factors in F{p}[{x},{y}]");
            Console.WriteLine($"{facts.Glue(" * ", "({0})")} = {P0}");
            facts1.Println($"Factors in Q[{x},{y}]");
            Console.WriteLine($"{facts1.Glue(" * ", "({0})")} = {P1}");
            Console.WriteLine();
            break;
        }
        catch (Exception)
        {
            Console.WriteLine($"########### P = {p} wont work");
        }
    }
}

void Run1()
{
    Ring.DisplayPolynomial = MonomDisplay.StarCaret;
    var (X2, X1) = Ring.Polynomial(ZnInt.ZnZero(101), MonomOrder.Lex, "X2", "X1").Deconstruct();

    // (X2^2 + 100*X1^2 + 100) * (X2^2 +  99*X2 + X1^2)
    var P = (X2.Pow(2) + 100 * X1.Pow(2) + 100) * (X2.Pow(2) + 99 * X2 + X1.Pow(2));

    var firr = Firr(P);
    var lifts = HenselLifting(P, firr);
    var facts = Recombinaison(P, lifts);

    var P0X2 = firr.Aggregate(P.One, (acc, c) => acc * c);
    Console.WriteLine($"P(X1,X2) = {P}");
    firr.Println($"P(0,X2) = {P0X2}");
    lifts.Println("Hensel Lifting");
    facts.Println("Factors");

    var P0 = facts.Aggregate(P.One, (acc, c) => acc * c);

    Console.WriteLine($"{facts.Glue(" * ", "({0})")} = {P0}");
    Console.WriteLine($"Check:{P.Equals(P0)}");
    Console.WriteLine();
}

// AECF p400-402
// P(X1,X2) = X2^4 +  99*X2^3 + 100*X2^2 +   2*X1^2*X2 +   2*X2 + 100*X1^4 + 100*X1^2
// P(0,X2) = X2^4 +  99*X2^3 + 100*X2^2 +   2*X2
//     X2
//     X2 + 100
//     X2 +  99
//     X2 + 1
// Hensel Lifting
//     X2 +  38*X1^4 +  50*X1^2 + 100
//     X2 +  38*X1^4 +  51*X1^2 +  99
//     X2 +  63*X1^4 +  50*X1^2
//     X2 +  63*X1^4 +  51*X1^2 + 1
// Factors
//     X2^2 + 100*X1^2 + 100
//     X2^2 +  99*X2 + X1^2
// (X2^2 + 100*X1^2 + 100) * (X2^2 +  99*X2 + X1^2) = X2^4 +  99*X2^3 + 100*X2^2 +   2*X1^2*X2 +   2*X2 + 100*X1^4 + 100*X1^2
// Check:True
// 

{
    Ring.DisplayPolynomial = MonomDisplay.StarCaret;
    var (X2, X1) = Ring.Polynomial(Rational.KZero(), MonomOrder.Lex, "X2", "X1").Deconstruct();

    // All Working 
    var polys = new[]
    {
        (X2.Pow(2) - X1.Pow(2) - 1) * (X2.Pow(2) - 2 * X2 + X1.Pow(2)),
        (X2.Pow(2) - 3 * X2 - 4 * X1) * (X2.Pow(2) - 4 * X1 * X2 - 4),
        (X2 - 3 * X1.Pow(2) - 2 * X1) * (X2.Pow(2) - 3 * X1 * X2 - 4),
        (X2.Pow(2) - 3 * X1 - 3) * (X2.Pow(2) - 4 * X2 - 3 * X1.Pow(2)),
        (X2 - 2 * X1 - 3) * (X2.Pow(2) - X2 - 2 * X1),
        (X2 - 2 * X1 - 4) * (X2.Pow(2) - X1.Pow(2) - 3),
        (X2 - 4 * X1.Pow(2) - 2 * X1) * (X2.Pow(2) - 4 * X1.Pow(2) - 3)
    };

    foreach (var F in polys)
        FactorsFxy(F);
}
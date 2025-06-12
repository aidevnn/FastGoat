using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Floats;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Perms;
using FastGoat.UserGroup.Polynoms;

namespace FastGoat.Examples;

public static class ResolventTheory
{
    static ResolventTheory()
    {
        Ring.DisplayPolynomial = MonomDisplay.StarCaret;
    }

    static KPoly<Rational> ResQuad(KPoly<Rational> P) => P.X.Pow(2) - Ring.Discriminant(P);

    // JS Milne, FT v5.10, page 50-51
    static KPoly<Rational> ResCub(KPoly<Rational> P)
    {
        var (_, b, c, d, e) = P.Coefs.Reverse().Deconstruct();
        var x = P.X;
        return x.Pow(3) - c * x.Pow(2) + (b * d - 4 * e) * x - b.Pow(2) * e + 4 * c * e - d.Pow(2);
    }

    // JS Milne, FT v5.10, page 50-51
    static void GalGrQuarticPol(KPoly<Rational> P)
    {
        if (P.Degree != 4 || !P.LC.Equals(P.KOne))
            throw new();

        var Pfacts = IntFactorisation.FirrZ2(P);
        if (Pfacts.Length != 1)
            throw new();

        Console.WriteLine($"P = {P}");
        Console.WriteLine("Chebotarev Gal(P) = {0}", GaloisApplicationsPart2.GaloisGroupChebotarev(P).Name);

        var Res = ResCub(P);
        var Rfacts = IntFactorisation.FirrZ2(Res);
        if (Rfacts.Length != 1)
            Console.WriteLine($"Res(P) = {Res} = {Rfacts.Glue("*", "({0})")}");
        else
            Console.WriteLine($"Res(P) = {Res} is irreductible");

        var RfactsDeg2or3 = Rfacts.Where(a => a.Degree > 1).ToArray();
        if (RfactsDeg2or3.Length == 0)
        {
            Console.WriteLine("Gal(P) = V4");
            Console.WriteLine();
            return;
        }

        var R0 = RfactsDeg2or3[0];
        var R0facts = IntFactorisation.AlgebraicFactors(R0);
        Console.WriteLine($"{R0} = {R0facts.Glue("*", "({0})")} in Q(y)");

        if (R0facts.Count != R0.Degree)
            Console.WriteLine("Gal(P) = S4");
        else if (R0facts.Count == 3)
            Console.WriteLine("Gal(P) = A4");
        else if (R0facts.Count == 2)
        {
            var P0 = P.Substitute(R0facts[0].X);
            var P0facts = IntFactorisation.AlgebraicFactors(P0);
            if (P0facts.Count != 1)
            {
                Console.WriteLine($"{P0} = {P0facts.Glue("*", "({0})")} in Q(y)");
                Console.WriteLine("Gal(P) = C4");
            }
            else
            {
                Console.WriteLine($"{P0} is irreductible in Q(y)");
                Console.WriteLine("Gal(P) = D4");
            }
        }

        Console.WriteLine();
    }

    static T HResF<T>(Perm[] sg, T X, Func<T[], T> f0, params T[] ts) where T : struct, IFieldElt<T>, IRingElt<T>, IElt<T>
    {
        var n = sg[0].Sn.N;
        if (sg.Any(p => p.Table.Length != n))
            throw new();

        return sg.Aggregate(X.One, (acc, perm) => acc * (X - f0(perm.Apply(ts))));
    }

    static T dn<T>(T[] t) where T : struct, IFieldElt<T>, IRingElt<T>, IElt<T>
    {
        var n = t.Length;
        var one = t[0].One;
        return n.Range().Grid2D().Where(e => e.t1 < e.t2).Aggregate(one, (acc, a) => acc * (t[a.t1] - t[a.t2]));
    }

    static T F0<T>(T[] t) where T : struct, IFieldElt<T>, IRingElt<T>, IElt<T>
    {
        return t[0] * t[2] + t[1] * t[3];
    }

    static T F1<T>(T[] t) where T : struct, IFieldElt<T>, IRingElt<T>, IElt<T>
    {
        // return t[0] - t[1];
        // return t[0] + t[1];
        return t[0] * t[1].Pow(2) * t[2].Pow(3) * t[3].Pow(4);
        return t[0] * t[2] + t[1] * t[3];
        return t[0] * t[1].Pow(2) + t[1] * t[2].Pow(2) + t[2] * t[3].Pow(2) + t[3] * t[0].Pow(2);
    }

    static void HResolventSymbolic(KPoly<Rational> P, Perm[] H, Func<Polynomial<Rational, Xi>[], Polynomial<Rational, Xi>> f0)
    {
        var n = P.Degree;
        var (X, xis) = Ring.Polynomial(Rational.KZero(), MonomOrder.GrLex, (n, "x"), "X");
        var Xi = X.ExtractIndeterminate;

        var diff = xis.Aggregate(X.One, (acc, xi) => acc * (X - xi)) - P.ToPolynomial(X);
        var dec = Ring.Decompose(diff, "X");
        var bsx = dec.Item1.Values.ToArray();
        var Q0 = HResF(H, X, f0, xis);

        var bs = bsx.Prepend(Q0).ToArray();
        bs.Println("System of symmetric polynomials");

        var rgb = Ring.ReducedGrobnerBasis(bs);
        rgb.Println("Reduced system");
        var Res = rgb.First(p => p.NbIndeterminates == 1 && p.ExtractIndeterminate.Equals(Xi)).ToKPoly(X);
        Console.WriteLine($"Res(P) = {Res}");
    }

    static void HResolventNumeric(KPoly<Rational> P, Perm[] H, Func<KPoly<BigCplx>[], KPoly<BigCplx>> f0, int digits = 30)
    {
        var n = P.Degree;
        var P0 = P.ToBcPoly(digits);
        var roots = FG.NRoots(P0);
        roots.Println("Roots");

        var proots = roots.Select(r => P0.One * r).ToArray();
        var Q0 = HResF(H, P0.X, f0, proots).ToBcPoly(digits * 2 / 3);
        var Q1 = new KPoly<Rational>(Q0.x, Rational.KZero(), Q0.Coefs.Select(c => c.RealPart.RoundEven.ToRational).ToArray());
        var factors = IntFactorisation.FactorsQ(Q1).OrderBy(e => e.Item1.Degree).ThenBy(e => e.Item1.NormB(2))
            .Select(e =>
            {
                if (e.Item2 == 1)
                    return e.Item1.Degree == 0 || e.Item1.Equals(e.Item1.X) ? $"{e.Item1}" : $"({e.Item1})";
                return e.Item1.Equals(e.Item1.X) ? $"{e.Item1}^{e.Item2}" : $"({e.Item1})^{e.Item2}";
            }).Glue(" * ");
        Console.WriteLine($"Res(P) = {Q1}");
        Console.WriteLine($"       = {factors}");
    }

    static void AnResSymb(KPoly<Rational> P)
    {
        var n = P.Degree;
        if (n < 3 || n > 5)
            throw new();

        var sn = FG.Symmetric(n);
        var an = FG.Alternate(n);
        var sigAn = Group.Cosets(sn, an, CosetType.Left).GroupBy(e => e.Value).Select(e => e.Key.X).ToArray();
        Console.WriteLine($"P = {P}");
        Console.WriteLine($"A{n}-Resolvent");
        HResolventSymbolic(P, sigAn, dn);
        Console.WriteLine($"ResQuad(P) = {ResQuad(P)}");
        Console.WriteLine();
    }

    static void AnResNum(KPoly<Rational> P)
    {
        var n = P.Degree;
        if (n != 3 && n != 4)
            throw new();

        var sn = FG.Symmetric(n);
        var an = FG.Alternate(n);
        var sigAn = Group.Cosets(sn, an, CosetType.Left).GroupBy(e => e.Value).Select(e => e.Key.X).ToArray();
        Console.WriteLine($"P = {P}");
        Console.WriteLine($"A{n}-Resolvent");
        HResolventNumeric(P, sigAn, dn);
        Console.WriteLine($"ResQuad(P) = {ResQuad(P)}");
        Console.WriteLine();
    }

    static void D8ResSymb(KPoly<Rational> P)
    {
        if (P.Degree != 4)
            throw new();

        var s4 = FG.Symmetric(4);
        var d8 = Group.Generate("D8", s4, s4[(1, 2, 3, 4)], s4[(1, 3)]);
        var sigD8 = Group.Cosets(s4, d8, CosetType.Left).GroupBy(e => e.Value).Select(e => e.Key.X).ToArray();
        Console.WriteLine($"P = {P}");
        Console.WriteLine($"D8-Resolvent");
        HResolventSymbolic(P, sigD8, F0);
        Console.WriteLine($"ResCub(P) = {ResCub(P)}");
        Console.WriteLine();
    }

    static void D8ResNum(KPoly<Rational> P)
    {
        if (P.Degree != 4)
            throw new();

        var s4 = FG.Symmetric(4);
        var d8 = Group.Generate("D8", s4, s4[(1, 2, 3, 4)], s4[(1, 3)]);
        var sigD8 = Group.Cosets(s4, d8, CosetType.Left).GroupBy(e => e.Value).Select(e => e.Key.X).ToArray();
        Console.WriteLine($"P = {P}");
        Console.WriteLine($"D8-Resolvent");
        HResolventNumeric(P, sigD8, F0);
        Console.WriteLine($"ResCub(P) = {ResCub(P)}");
        Console.WriteLine();
    }

    static void V4ResNum(KPoly<Rational> P)
    {
        if (P.Degree != 4)
            throw new();

        var s4 = FG.Symmetric(4);
        var v4 = Group.Generate("V4", s4, s4[(1, 2), (3, 4)], s4[(1, 3), (2, 4)]);
        // var v4 = Group.Generate("V4", s4, s4[(1, 2, 3, 4)]);
        var sigV4 = Group.Cosets(s4, v4, CosetType.Left).GroupBy(e => e.Value).Select(e => e.Key.X).ToArray();
        Console.WriteLine($"P = {P}");
        Console.WriteLine($"V4-Resolvent");
        HResolventNumeric(P, sigV4, F1);
        // Console.WriteLine($"ResCub(P) = {ResCub(P)}");
        Console.WriteLine();
    }

    public static void JSMilneResolventCubic()
    {
        var x = FG.QPoly();
        GalGrQuarticPol(x.Pow(4) - 4 * x.Pow(2) + 1);
        GalGrQuarticPol(x.Pow(4) - 10 * x.Pow(2) + 4);
        GalGrQuarticPol(x.Pow(4) + 1);

        GalGrQuarticPol(x.Pow(4) + x.Pow(3) + x.Pow(2) + x + 1);
        GalGrQuarticPol(x.Pow(4) + 4 * x.Pow(2) + 2);
        GalGrQuarticPol(x.Pow(4) - 5 * x.Pow(2) + 5);

        GalGrQuarticPol(x.Pow(4) - 2);
        GalGrQuarticPol(x.Pow(4) - 2 * x.Pow(2) + 3);
        GalGrQuarticPol(x.Pow(4) + x.Pow(2) - 1);

        GalGrQuarticPol(x.Pow(4) + 8 * x + 12);

        GalGrQuarticPol(x.Pow(4) - 4 * x + 2);
        GalGrQuarticPol(x.Pow(4) + x + 1);
        GalGrQuarticPol(x.Pow(4) + x.Pow(3) - 5);
    }

    public static void ExamplesHResolventSymbolic()
    {
        var x = FG.QPoly();
        AnResSymb(x.Pow(3) + x + 1);
        AnResSymb(x.Pow(3) - 3 * x.Pow(2) + 5);

        D8ResSymb(x.Pow(4) - 4 * x.Pow(2) + 2);
        D8ResSymb(x.Pow(4) + x.Pow(2) + 1);
        D8ResSymb(x.Pow(4) + x + 1);

        // AnResSymb(x.Pow(4) + x.Pow(2) + 1); // very slow ~3min
        // AnResSymb(x.Pow(4) + x + 1); // very slow ~3min
    }

    public static void ExamplesHResolventNumeric()
    {
        var x = FG.QPoly();
        AnResNum(x.Pow(3) + x + 1);
        AnResNum(x.Pow(3) - 3 * x.Pow(2) + 5);

        AnResNum(x.Pow(4) + x.Pow(2) + 1);
        AnResNum(x.Pow(4) + 1);
        AnResNum(x.Pow(4) + x + 1);

        D8ResNum(x.Pow(4) - 4 * x.Pow(2) + 2);
        D8ResNum(x.Pow(4) + x.Pow(2) + 1);
        D8ResNum(x.Pow(4) + 1);
        D8ResNum(x.Pow(4) + x + 1);

        V4ResNum(x.Pow(4) + x.Pow(3) + x.Pow(2) + x + 1);
        V4ResNum(x.Pow(4) + 1);
        V4ResNum(x.Pow(4) - 2);
        V4ResNum(x.Pow(4) - 8 * x + 12);
        V4ResNum(x.Pow(4) + x + 1);
    }
}
using System.Diagnostics;
using System.Numerics;
using FastGoat.Commons;
using FastGoat.Examples;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Characters;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Matrix;
using FastGoat.UserGroup.Perms;
using FastGoat.UserGroup.Polynoms;
using FastGoat.UserGroup.Words;
using FastGoat.UserGroup.Words.ToddCoxeter;
using static FastGoat.Commons.IntExt;
using static FastGoat.Commons.EnumerableExt;
using FastGoat.UserGroup.Padic;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");
Ring.DisplayPolynomial = MonomDisplay.StarCaret;

KPoly<Rational> TschirnhausenTransformation(KPoly<Rational> P)
{
    var n = P.Degree;
    var y = FG.EPoly(P);
    for (int i = 0; i < 50; i++)
    {
        var Q = PolynomialFactorization.RandPoly(Rational.KZero(), 10, new[] { n - 1 });
        var U = GaloisTheory.MinPolynomial(Q.Substitute(y), P.x);
        if (Ring.Gcd(U, U.Derivative).Degree == 0)
            return Q;
    }

    throw new();
}

T HResF<T>(Perm[] sg, T X, Func<T[], T> f0, params T[] ts) where T : struct, IFieldElt<T>, IRingElt<T>, IElt<T>
{
    var n = sg[0].Sn.N;
    if (sg.Any(p => p.Table.Length != n))
        throw new();

    return sg.Aggregate(X.One, (acc, perm) => acc * (X - f0(perm.Apply(ts))));
}

(KPoly<Rational>, KPoly<BigCplx>[]) HResolventNumeric(KPoly<Rational> P, Perm[] H, Func<KPoly<BigCplx>[], KPoly<BigCplx>> f0,
    int digits = 40, int err = 10)
{
    var n = P.Degree;
    var P0 = P.ToBcPoly(digits);
    var roots = FG.NRoots(P0);
    roots.Println("Roots");

    var proots = roots.Select(r => P0.One * r).ToArray();
    var Q0 = HResF(H, P0.X, f0, proots).ToIntPoly(err);
    Console.WriteLine($"Res(P) = {Q0}");
    return (Q0, roots.Select(r => r * P0.One).ToArray());
}

void CohenDegree5(KPoly<Rational> P, int O = 40)
{
    if (P.Degree != 5 || !P.LT.Equals(P.KOne))
        throw new();

    KPoly<BigCplx> F1(KPoly<BigCplx>[] t)
    {
        var (x1, x2, x3, x4, x5) = t.Deconstruct();
        return x1.Pow(2) * (x2 * x5 + x3 * x4) + x2.Pow(2) * (x1 * x3 + x4 * x5) + x3.Pow(2) * (x1 * x5 + x2 * x4)
               + x4.Pow(2) * (x1 * x2 + x3 * x5) + x5.Pow(2) * (x1 * x4 + x2 * x3);
    }

    KPoly<BigCplx> F2(KPoly<BigCplx>[] t)
    {
        var (x1, x2, x3, x4, x5) = t.Deconstruct();
        return x1 * x2.Pow(2) + x2 * x3.Pow(2) + x3 * x4.Pow(2) + x4 * x5.Pow(2) + x5 * x1.Pow(2);
    }

    var s5 = FG.Symmetric(5);
    var c5 = Group.Generate("C5", s5, s5[(1, 2, 3, 4, 5)]);
    var m20 = Group.Generate("M20", s5, s5[(1, 2, 3, 4, 5)], s5[(2, 3, 5, 4)]);
    var sigM20 = Group.Cosets(s5, m20, Perm.OrbitsComparer, CosetType.Left).GroupBy(e => e.Value).Select(e => e.Key.X).ToArray();
    var sigD10 = Group.Cosets(FG.Dihedral(5), c5, Perm.OrbitsComparer, CosetType.Left).GroupBy(e => e.Value).Select(e => e.Key.X)
        .ToArray();

    Console.WriteLine($"P = {P}");
    var (R, roots) = HResolventNumeric(P, sigM20, F1, O, O / 2);
    var sqfrR = Ring.Gcd(R, R.Derivative).Degree == 0;
    Console.WriteLine($"P Is SquareFree {sqfrR}");
    if (!sqfrR)
        throw new($"Res isnt square free");

    var facts = IntFactorisation.FactorsQ(R);
    facts.Println("Factors");
    var discP = Ring.Discriminant(P);
    var discIsSquare = discP.IsSquare;
    Console.WriteLine($"Disc(P) = {discP} is square {discIsSquare}");

    if (facts.Length == 1)
    {
        if (!discIsSquare)
            Console.WriteLine("Gal(P) = S5");
        else
            Console.WriteLine("Gal(P) = A5");
    }
    else
    {
        if (!discIsSquare)
            Console.WriteLine("Gal(P) = M20");
        else
        {
            var d = roots[0].One;
            var r0 = -facts.First(e => e.Item1.Degree == 1).Item1.ToBcPoly(O / 2)[0];
            var perm = s5.First(perm => (F1(perm.Apply(roots))[0].ToBigCplx(O / 2) - r0).IsZero4d());
            Console.WriteLine(new { r0, perm });
            var roots2 = perm.Apply(roots);
            var (r1, r2, r3, r4, r5) = roots2.Deconstruct();
            for (int k = 0; k < 50; ++k)
            {
                d = (r1 * r2 * (r2 - r1) + r2 * r3 * (r3 - r2) + r3 * r4 * (r4 - r3) + r4 * r5 * (r5 - r4) + r5 * r1 * (r1 - r5))
                    .Pow(2);
                
                if (!d.IsZero())
                    break;

                var A = TschirnhausenTransformation(P);
                Console.WriteLine($"Tschirnhausen Transformation A = {A}");
                roots2 = roots2.Select(r => A.Substitute(r)).ToArray();
                (r1, r2, r3, r4, r5) = roots2.Deconstruct();
            }

            var d0 = d[0].RoundEven;
            var dr = d0.RealPart.ToRational;
            Console.WriteLine($"d = {d} ~ {d0} = {dr} is square {dr.IsSquare}");
            if(dr.IsSquare)
                Console.WriteLine("Gal(P) = C5");
            else
                Console.WriteLine("Gal(P) = D5");
        }
    }

    Console.WriteLine();
}

{
    var x = FG.QPoly();

    var polys = new[]
    {
        x.Pow(5) - x + 1, 
        x.Pow(5) + 20 * x + 16,
        x.Pow(5) + 2,
        x.Pow(5) - 5 * x + 12,
        x.Pow(5) + x.Pow(4) - 4 * x.Pow(3) - 3 * x.Pow(2) + 3 * x + 1
    };

    foreach (var P in polys)
    {
        CohenDegree5(P);
        var lt = 100.Range().Select(i => Rng.Next(1, 4) * (-1).Pow(Rng.Next(1, 3))).Distinct().Take(3).ToArray();
        foreach(var r in lt)
            CohenDegree5(P.Substitute(x - r));

        Console.WriteLine("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
    }
}

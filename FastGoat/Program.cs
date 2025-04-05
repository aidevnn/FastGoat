using System.Numerics;
using System.Reflection;
using FastGoat.Commons;
using FastGoat.Examples;
using FastGoat.Structures;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.EllCurve;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.LWE;
using FastGoat.UserGroup.Matrix;
using FastGoat.UserGroup.Padic;
using FastGoat.UserGroup.Polynoms;
using Rq = FastGoat.Structures.VecSpace.KPoly<FastGoat.UserGroup.Integers.Rational>;
using static FastGoat.Commons.IntExt;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");
Ring.DisplayPolynomial = MonomDisplay.StarCaret;
// RecomputeAllPrimesUpTo(5000000);

void AddNewIntegralPoint(EllGroup<Rational> g, HashSet<EllPt<Rational>> set, EllPt<Rational> pt)
{
    var sz = 0;
    while (set.Count != sz)
    {
        sz = set.Count;
        var tmp = set.Select(e => g.Op(e, pt)).Where(e => e.IsIntegral()).ToHashSet();
        set.UnionWith(tmp);
    }
}

(BigInteger M, BigInteger e) SolveEq(BigInteger b1, BigInteger b2, int nMax = 200)
{
    if (b1 < 0 && b2 < 0)
        return (0, 0);

    foreach (var e in nMax.SeqLazy(1).Where(e => GcdBigInt(e, b1) == 1))
    {
        foreach (var M in nMax.SeqLazy(1).Where(M => GcdBigInt(e, M) == 1 && GcdBigInt(M, b2) == 1))
        {
            var N2 = (b1 * BigInteger.Pow(M, 4) + b2 * BigInteger.Pow(e, 4));
            // Console.WriteLine(new { M, e, N2 });
            if (N2 < 0) continue;
            var N = SqrtBigInt(N2);
            if (N * N == N2 && GcdBigInt(N, e) == 1)
                return (M, e);
        }
    }

    return (0, 0);
}

Dictionary<EllPt<Rational>, HashSet<EllPt<Rational>>> EllGoverH(EllGroup<Rational> E, HashSet<EllPt<Rational>> G,
    HashSet<EllPt<Rational>> H)
{
    var cosets = new Dictionary<EllPt<Rational>, HashSet<EllPt<Rational>>>();
    cosets.Add(E.O, H);
    var maxH = G.Max(pt => pt.Height());
    var rem = G.Except(H).ToHashSet();
    while (rem.Count > 0)
    {
        var pt = rem.MinBy(pt => pt.Height());
        var cos = H.Select(pt1 => E.Op(pt, pt1)).Where(pt0 => pt0.Height() <= maxH).ToHashSet();
        cosets.Add(cos.MinBy(pt1 => pt1.Height()), cos);
        rem.ExceptWith(cos);
    }

    return cosets;
}

HashSet<EllPt<Rational>> Descent(EllGroup<Rational> E, HashSet<EllPt<Rational>> G, HashSet<EllPt<Rational>> H) =>
    EllGoverH(E, G, H).Keys.ToHashSet(EqualityComparer<EllPt<Rational>>.Create(
        (pt0, pt1) => pt0.X.Equals(pt1.X) && pt0.Y.Equals(-pt1.Y),
        pt => pt.X.GetHashCode())
    );

HashSet<EllPt<Rational>> Independants(ConcreteGroup<EllPt<Rational>> E, HashSet<EllPt<Rational>> pts, BigInteger maxH)
{
    var bsE = (EllGroup<Rational>)E.BaseGroup;
    var gens = E.GetGenerators().OrderBy(pt => pt.Height()).ThenBy(pt => !pt.IsO ? pt.X.Absolute : "0").ToHashSet();
    var setPts = E.ToHashSet();
    foreach (var pt in pts.OrderBy(pt => pt.Height()).ThenByDescending(pt => pt.IsIntegral()))
    {
        if (setPts.Add(pt))
        {
            gens.Add(pt);
            var pti = E.Invert(pt);
            if (pt.IsIntegral())
            {
                AddNewIntegralPoint(bsE, setPts, pt);
                AddNewIntegralPoint(bsE, setPts, pti);
            }

            var tmp = setPts.Concat([pt, pti]).ToArray();
            setPts.UnionWith(tmp.Grid2D(tmp).Select(e => E.Op(e.t1, e.t2)).Where(pt0 => pt0.Height() <= maxH)
                .ToHashSet());
        }
    }

    return gens;
}

(EllGroup<Rational> E, ConcreteGroup<EllPt<Rational>> gEll, HashSet<EllPt<Rational>> intPts, HashSet<EllPt<Rational>>
    gens)
    EllGenerators(BigInteger b, int nMax = 200)
{
    var ng = EllipticCurves.NagellLutzTorsionGroup(b, 0);
    var x = FG.QPoly();
    Console.WriteLine($"Elliptic curve y^2 = {x.Pow(3) + new Rational(b) * x}");
    var divs_b = DividorsBigInt(BigInteger.Abs(b)).Order().ToArray();
    var sols = divs_b.SelectMany(b1 => new[] { b1, -b1 }).Select(b1 => (b1, s: SolveEq(b1, b / b1, nMax)))
        .Where(e => e.s.M != 0)
        .Select(e => (e.b1, X: new Rational(e.b1 * BigInteger.Pow(e.s.M, 2), BigInteger.Pow(e.s.e, 2))))
        .Select(e => (e.b1, e.X, Y2: e.X.Pow(3) + b * e.X))
        .Select(e => (e.b1, e.X, e.Y2, Y: Rational.Sqrt(e.Y2)))
        .ToArray();
    var listPts = sols.Select(e => new EllPt<Rational>(e.X, e.Y)).SelectMany(e => new[] { e, ng.E.Invert(e) })
        .ToHashSet();

    foreach (var pt in listPts.Where(pt => pt.IsIntegral()).ToHashSet())
        AddNewIntegralPoint(ng.E, ng.intPts, pt);

    DisplayGroup.HeadElements(ng.gEll);

    var setG = listPts.Append(ng.E.O).Union(ng.intPts).ToHashSet();
    var set2G = setG.ToDictionary(pt => pt, pt => ng.E.Times(pt, 2)).Select(e => e.Value).ToHashSet();
    var desc1 = Descent(ng.E, setG, set2G);
    var desc2 = Descent(ng.E, desc1, ng.gEll.ToHashSet());
    var gens = Independants(ng.gEll, desc2, setG.Max(pt => pt.Height()));

    Console.WriteLine($"Start[{setG.Count}] Descent1[{desc1.Count}] Descent2[{desc2.Count}] Gens[{gens.Count}]");
    Console.WriteLine($"Elliptic curve y^2 = {x.Pow(3) + new Rational(b) * x}");
    var rank = gens.Count - ng.gEll.GetGenerators().Count();
    Console.WriteLine($"{ng.E} Rank = {rank} Gens = [{gens.Glue(", ")}]");
    Console.WriteLine();

    return (ng.E, ng.gEll, ng.intPts, gens);
}

void testRank1()
{
    Ring.DisplayPolynomial = MonomDisplay.StarCaret;
    var seq = new[] { 1, 5, 6, 7, 14, 15, 21, 22, 29, 30, 31, 34, 78, 210, 1254, 29274 }.Select(n => n * n).ToArray();

    GlobalStopWatch.Restart();
    var res = seq.Select(b => EllGenerators(-b)).ToArray();
    GlobalStopWatch.Show();

    Console.WriteLine();

    foreach (var ng in res)
        Console.WriteLine(
            $"{ng.E,-25} rank = {ng.gens.Count - ng.gEll.GetGenerators().Count()} gens = [{ng.gens.Glue(", ")}]");

    Console.WriteLine();

    GlobalStopWatch.Show(); // Time:36.398s
}

void testRank2()
{
    Ring.DisplayPolynomial = MonomDisplay.StarCaret;
    var seq = new[] { -1, 1, -3, 3, -5, 5, 7, 17, 73, -82 };

    GlobalStopWatch.Restart();
    var res = seq.Select(b => EllGenerators(b)).ToArray();
    GlobalStopWatch.Show();

    Console.WriteLine();

    foreach (var ng in res)
        Console.WriteLine(
            $"{ng.E,-25} rank = {ng.gens.Count - ng.gEll.GetGenerators().Count()} gens = [{ng.gens.Glue(", ")}]");

    Console.WriteLine();

    GlobalStopWatch.Show();
}

{
    // testRank1();
    testRank2();
}
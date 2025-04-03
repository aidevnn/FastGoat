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
RecomputeAllPrimesUpTo(5000000);

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

BigInteger SquareFree(BigInteger a)
{
    if (a == 0 || a * a.Sign == 1)
        return a;

    return a.Sign * PrimesDec(a.Sign * a).Select(e => BigInteger.Pow(e.Key, e.Value % 2))
        .Aggregate(BigInteger.One, (ai, aj) => ai * aj);
}

Rational SquareFreeRat(Rational a) => new(SquareFree(a.Num), SquareFree(a.Denom));

Rational Alpha(EllGroup<Rational> E, EllPt<Rational> pt)
{
    if (pt.IsO)
        return new(1);

    if (pt.X.IsZero())
        return SquareFreeRat(E.A);

    return SquareFreeRat(pt.X);
}

HashSet<Rational> AllAlpha(HashSet<Rational> gens)
{
    var sz = 0;
    var set = gens.ToHashSet();
    while (set.Count != sz)
    {
        sz = set.Count;
        var tmp = set.Grid2D(gens).Select(e => SquareFreeRat(e.t1 * e.t2)).ToHashSet();
        set.UnionWith(tmp);
    }

    return set;
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
            if (N * N == N2)
                return (M, e);
        }
    }

    return (0, 0);
}

List<EllPt<Rational>> Descent1(EllGroup<Rational> E, EllPt<Rational> P,
    Dictionary<EllPt<Rational>, HashSet<EllPt<Rational>>> cosets)
{
    var Q0 = P;
    var list = new List<EllPt<Rational>>();

    do
    {
        list.Add(Q0);
        var h = Q0.Height();
        var Q1 = cosets.Keys.ToDictionary(Qi => Qi, Qi => E.Op(Q0, E.Invert(Qi)))
            .Where(e => e.Value.Height() < h)
            .OrderBy(e => e.Key.Height())
            .Select(e => e.Key)
            .FirstOrDefault(new EllPt<Rational>());

        Q0 = cosets[Q1].MinBy(e => e.Height());
    } while (!Q0.IsO);

    return list;
}

HashSet<EllPt<Rational>> Descent2(ConcreteGroup<EllPt<Rational>> E, HashSet<EllPt<Rational>> pts)
{
    var setPts = new HashSet<EllPt<Rational>>();
    foreach (var pt in pts.Except(E).OrderBy(pt => pt.Height()).ThenBy(pt => !pt.IsIntegral()))
    {
        var pt1 = pt;
        while (true)
        {
            var pt2 = E.Select(pt0 => E.Op(pt0, pt1)).MinBy(pt0 => pt0.Height());
            if (pt2.Height() < pt1.Height())
            {
                pt1 = pt2;
            }
            else
            {
                setPts.Add(pt1);
                break;
            }
        }
    }

    return setPts;
}

(EllGroup<Rational> E, ConcreteGroup<EllPt<Rational>> gEll, HashSet<EllPt<Rational>> intPts, HashSet<EllPt<Rational>>)
    EllGenerators(BigInteger b)
{
    var ng = EllipticCurves.NagellLutzTorsionGroup(b, 0);
    var x = FG.QPoly();
    Console.WriteLine($"Elliptic curve y^2 = {x.Pow(3) + new Rational(b) * x}");
    var divs_b = DividorsBigInt(BigInteger.Abs(b)).Order().ToArray();
    var sols = divs_b.SelectMany(b1 => new[] { b1, -b1 }).Select(b1 => (b1, Me: SolveEq(b1, b / b1)))
        .Where(e => e.Me.M != 0)
        .Select(e => (e.b1, X: new Rational(e.b1 * BigInteger.Pow(e.Me.M, 2), BigInteger.Pow(e.Me.e, 2))))
        .Select(e => (e.b1, e.X, Y2: e.X.Pow(3) + b * e.X))
        .Select(e => (e.b1, e.X, e.Y2, Y: Rational.Sqrt(e.Y2)))
        .ToArray();
    var listPts = sols.Select(e => new EllPt<Rational>(e.X, e.Y)).SelectMany(e => new[] { e, ng.E.Invert(e) })
        .ToHashSet();

    foreach (var pt in listPts.Where(pt => pt.IsIntegral()).ToHashSet())
        AddNewIntegralPoint(ng.E, ng.intPts, pt);

    DisplayGroup.HeadElements(ng.gEll);
    // ng.intPts.OrderBy(e => e.Height()).ThenBy(e => e).Println("Integral points");

    var set0 = listPts.Append(ng.E.O).Union(ng.intPts).ToHashSet();
    var cosets = set0.Select(e => (e, e2: ng.E.Op(e, e))).GroupBy(e => e.e2)
        .ToDictionary(e => e.Key, e => e.Select(f => f.e).ToHashSet());
    var set1 = set0.OrderBy(pt => !pt.IsIntegral()).ThenBy(pt => pt.Height())
        .Select(pt => Descent1(ng.E, pt, cosets).Last()).ToHashSet();
    var set2 = Descent2(ng.gEll, set1);

    Console.WriteLine($"Start[{set0.Count}] Descent1[{set1.Count}] Descent2[{set2.Count}]");
    var setPts = ng.gEll.ToHashSet();
    var gens = ng.gEll.GetGenerators().ToHashSet();
    foreach (var pt in set2.OrderBy(pt => !pt.IsIntegral()).ThenBy(pt => pt.Height()))
    {
        if (setPts.Add(pt))
        {
            gens.Add(pt);
            var pti = ng.E.Invert(pt);
            if (pt.IsIntegral())
            {
                AddNewIntegralPoint(ng.E, setPts, pt);
                AddNewIntegralPoint(ng.E, setPts, pti);
            }

            setPts.UnionWith(setPts.Select(pt0 => ng.E.Op(pt0, pt)).ToHashSet());
            setPts.UnionWith(setPts.Select(pt0 => ng.E.Op(pt0, pti)).ToHashSet());
        }
    }

    var rank1 = gens.Count - ng.gEll.GetGenerators().Count();
    var rank2 = EllipticCurves.EllRank(b, 0); // instable
    Console.WriteLine($"{ng.E} Gens:[{gens.Glue(", ")}] Rank = {rank1} / {rank2} {(rank1 == rank2 ? "PASS" : "FAIL")}");
    Console.WriteLine();

    return (ng.E, ng.gEll, ng.intPts, gens);
}

HashSet<BigInteger> TorsionFree(ConcreteGroup<EllPt<Rational>> g, HashSet<EllPt<Rational>> set)
{
    if (!set.Any())
        return new();

    var setRemIntPts = set.OrderBy(e => g.Contains(e) ? 0 : 1).ThenBy(e => e).ToHashSet();
    var setRankGens = new List<EllPt<Rational>>();
    var setIntPts = g.ToHashSet();
    setRemIntPts.ExceptWith(g);
    var E = (EllGroup<Rational>)g.BaseGroup;
    while (setRemIntPts.Count != 0)
    {
        var e0 = setRemIntPts.MaxBy(t0 =>
        {
            var tmp = setIntPts.Union([t0, g.Invert(t0)]).ToHashSet();
            foreach (var pt in tmp.ToArray())
                AddNewIntegralPoint(E, tmp, pt);

            return tmp.Count;
        });

        var e1 = g.Invert(e0);
        if (!g.Contains(e0))
            setRankGens.Add(e0);

        setIntPts.UnionWith([e0, e1]);
        foreach (var pt in setIntPts.ToArray())
            AddNewIntegralPoint(E, setIntPts, pt);

        setRemIntPts.ExceptWith(setIntPts);
    }

    return setRankGens.Select(e => e.X.Num).ToHashSet();
}

HashSet<EllPt<Rational>> SearchTorsionFree(int b, int rk, HashSet<EllPt<Rational>> intPts, int nMax = 100000)
{
    var setGensInt = new HashSet<BigInteger>();
    var k = 0;
    BigInteger b1 = 0;
    foreach (var k0 in nMax.SeqLazy(1))
    {
        k = k0;
        b1 = b * BigInteger.Pow(k0, 4);
        var ng1 = EllipticCurves.NagellLutzTorsionGroup(b1, 0);
        var xs = intPts.Select(e => e.IsO ? e : new EllPt<Rational>(e.X * k0 * k0, e.Y * k0 * k0 * k0)).ToHashSet();
        foreach (var pt in xs.ToArray())
            AddNewIntegralPoint(ng1.E, xs, pt);

        var rks = TorsionFree(ng1.gEll, xs);

        var nMin = int.Min((int)SqrtBigInt(BigInteger.Abs(b1)), nMax / 2);
        Console.WriteLine($"{ng1.E} ---> {ng1.E} with k={k0}");
        foreach (var (x, y2) in nMin.SeqLazy(step: -1).Concat((nMax / 2).SeqLazy(nMin))
                     .Select(x => (x: new BigInteger(x), y2: BigInteger.Pow(x, 3) + b1 * x)))
        {
            if (rks.Count == rk)
                break;

            if (xs.Any(e => !e.IsO & e.X.Num == x) || y2 <= 0)
                continue;

            var y = SqrtBigInt(y2);
            if (y2 != y * y)
                continue;

            var pt0 = new EllPt<Rational>(new(x), new(y));
            var rks0 = TorsionFree(ng1.gEll, xs.Union([pt0, ng1.E.Invert(pt0)]).ToHashSet());
            if (rks0.Count == rks.Count)
                continue;

            AddNewIntegralPoint(ng1.E, xs, pt0);
            AddNewIntegralPoint(ng1.E, xs, ng1.E.Invert(pt0));
            rks.Clear();
            rks.UnionWith(rks0);
        }

        if (rks.Count == rk)
        {
            setGensInt = rks;
            break;
        }
    }

    var (k2, k3) = (new Rational(k).Pow(2), new Rational(k).Pow(3));
    return setGensInt.Select(x => (x, y: SqrtBigInt(BigInteger.Pow(x, 3) + b1 * x)))
        .Select(e => new EllPt<Rational>(new Rational(e.x) / k2, new Rational(e.y) / k3))
        .ToHashSet();
}

EllPt<Rational> PsiEll(EllPt<Rational> pt, Rational b, int k = 1)
{
    if (pt.IsO || pt.X.IsZero())
        return new();

    var (X2, Y2) = (pt.X.Pow(2), pt.Y.Pow(2));
    var (x, y) = (Y2 / (k * k * X2), pt.Y * (X2 - b) / (k.Pow(3) * X2));
    return new(x, y);
}

EllPt<Rational> PhiEll(EllPt<Rational> pt, Rational b) => PsiEll(pt, b, 2);

(EllGroup<Rational> E, ConcreteGroup<EllPt<Rational>> gEll, HashSet<EllPt<Rational>> setGens,
    HashSet<EllPt<Rational>> intPts) IntegralPoints(int b)
{
    var x = FG.QPoly();
    Console.WriteLine($"Elliptic curve y^2 = {x.Pow(3) + b * x}");
    var rank = EllipticCurves.EllRank(b, 0); // instable
    var ng = EllipticCurves.NagellLutzTorsionGroup(b, 0);
    DisplayGroup.HeadElements(ng.gEll);

    var intPts = ng.intPts.Union(ng.gEll).ToHashSet();
    foreach (var pt in ng.intPts)
        AddNewIntegralPoint(ng.E, intPts, pt);

    var setGens = SearchTorsionFree(b, rank, intPts);
    Console.WriteLine($"Elliptic curve y^2 = {x.Pow(3) + b * x}");
    Console.WriteLine($"{ng.E} pts ord infty [{setGens.Glue()}]");

    foreach (var pt in setGens)
        AddNewIntegralPoint(ng.E, intPts, pt);

    return (ng.E, ng.gEll, setGens, intPts);
}

void craft1()
{
    Ring.DisplayPolynomial = MonomDisplay.StarCaret;
    GlobalStopWatch.Restart();
    var seq = new[] { 1, -1, 5, 25, 34.Pow(2), 1254.Pow(2), 29274.Pow(2) };
    // var seq = new[] { 1, 5, 6, 7, 14, 15, 21, 22, 29, 30, 31, 34, 78, 210, 1254 }.Select(n => n * n)
    //     .Concat([-1, 5]).Order().ToArray();
    var x = FG.QPoly();
    foreach (var b in seq)
    {
        GlobalStopWatch.AddLap();
        var seqAlphas = new HashSet<Rational>();
        var seqAlphas_ = new HashSet<Rational>();
        var rank1 = 0;

        Console.WriteLine($"Start y^2 = {x.Pow(3) - b * x}");
        {
            var divs_b = DividorsBigInt(BigInteger.Abs(b)).Order().ToArray();
            var ng = IntegralPoints(-b);
            rank1 = ng.setGens.Count;

            var dico1 = ng.intPts.ToDictionary(a => a, a => PsiEll(a, ng.E.A));
            var dico2a = dico1.ToDictionary(a => a.Key, a => PhiEll(a.Value, -4 * ng.E.A));
            dico2a.Println("P -> 2P");
            if (dico2a.Any(e => !e.Value.Equals(ng.E.Times(e.Key, 2))))
            {
                dico2a.ToDictionary(e => e.Key, e => (e.Value, ng.E.Times(e.Key, 2))).Println();
                throw new();
            }

            var dicoAlphas = ng.intPts.Union(ng.setGens).ToDictionary(e => e, e => Alpha(ng.E, e));
            seqAlphas.UnionWith(AllAlpha(dicoAlphas.Values.ToHashSet()));
            var checkHom = dicoAlphas.Grid2D(dicoAlphas)
                .Select(e => ((e.t1.Key, e.t2.Key),
                    (Alpha(ng.E, ng.E.Op(e.t1.Key, e.t2.Key)), SquareFreeRat(e.t1.Value * e.t2.Value))))
                .ToArray();
            if (checkHom.Any(e => !e.Item2.Item1.Equals(e.Item2.Item2)))
            {
                checkHom.Where(e => !e.Item2.Item1.Equals(e.Item2.Item2)).Println("Homomorphism");
                throw new();
            }

            Console.WriteLine($"Divs:[{divs_b.Glue(", ")}]");
            var sols = divs_b.SelectMany(b1 => new[] { b1, -b1 }).Select(b1 => (b1, Me: SolveEq(b1, -b / b1)))
                .Where(e => e.Me.M != 0)
                .Select(e => (e.b1, X: new Rational(e.b1 * BigInteger.Pow(e.Me.M, 2), BigInteger.Pow(e.Me.e, 2))))
                .Select(e => (e.b1, e.X, Y2: e.X.Pow(3) - b * e.X))
                .Select(e => (e.b1, e.X, e.Y2, Y: Rational.Sqrt(e.Y2)))
                .ToArray();
            var sols2 = sols.Select(e => new EllPt<Rational>(e.X, e.Y)).ToHashSet();

            sols.Select(e => $"b1:{e.b1} X:{e.X} Y2:{e.Y2} {e.Y2.IsSquare}")
                .Println($"Check:{sols2.IsSupersetOf(ng.setGens)}");
        }

        Console.WriteLine($"Continue y^2 = {x.Pow(3) + 4 * new Rational(b) * x}");
        {
            var divs_b = DividorsBigInt(BigInteger.Abs(4 * b)).Order().ToArray();
            var ng = IntegralPoints(4 * b);

            var dico1 = ng.intPts.ToDictionary(a => a, a => PsiEll(a, ng.E.A));
            var dico2a = dico1.ToDictionary(a => a.Key, a => PhiEll(a.Value, -4 * ng.E.A));
            dico2a.Println("P -> 2P");
            if (dico2a.Any(e => !e.Value.Equals(ng.E.Times(e.Key, 2))))
            {
                dico2a.ToDictionary(e => e.Key, e => (e.Value, ng.E.Times(e.Key, 2))).Println();
                throw new();
            }

            var dicoAlphas = ng.intPts.Union(ng.setGens).ToDictionary(e => e, e => Alpha(ng.E, e));
            seqAlphas_.UnionWith(AllAlpha(dicoAlphas.Values.ToHashSet()));
            var checkHom = dicoAlphas.Grid2D(dicoAlphas)
                .Select(e => ((e.t1.Key, e.t2.Key),
                    (Alpha(ng.E, ng.E.Op(e.t1.Key, e.t2.Key)), SquareFreeRat(e.t1.Value * e.t2.Value))))
                .ToArray();
            if (checkHom.Any(e => !e.Item2.Item1.Equals(e.Item2.Item2)))
            {
                checkHom.Where(e => !e.Item2.Item1.Equals(e.Item2.Item2)).Println("Homomorphism");
                throw new();
            }

            Console.WriteLine($"Divs:[{divs_b.Glue(", ")}]");
            var sols = divs_b.SelectMany(b1 => new[] { b1, -b1 }).Select(b1 => (b1, Me: SolveEq(b1, 4 * b / b1)))
                .Where(e => e.Me.M != 0)
                .Select(e => (e.b1, X: new Rational(e.b1 * BigInteger.Pow(e.Me.M, 2), BigInteger.Pow(e.Me.e, 2))))
                .Select(e => (e.b1, e.X, Y2: e.X.Pow(3) + 4 * b * e.X))
                .Select(e => (e.b1, e.X, e.Y2, Y: Rational.Sqrt(e.Y2)))
                .ToArray();
            var sols2 = sols.Select(e => new EllPt<Rational>(e.X, e.Y)).ToHashSet();

            sols.Select(e => $"b1:{e.b1} X:{e.X} Y2:{e.Y2} {e.Y2.IsSquare}")
                .Println($"Check:{sols2.IsSupersetOf(ng.setGens)}");
        }

        Console.WriteLine($"#Alpha    = {seqAlphas.Count,3} [{seqAlphas.Order().Glue(", ")}]");
        Console.WriteLine($"#AlphaBar = {seqAlphas_.Count,3} [{seqAlphas_.Order().Glue(", ")}]");
        var rank2 = int.Log2(seqAlphas.Count * seqAlphas_.Count / 4);
        Console.WriteLine("rank = {0} / {1} {2}", rank1, rank2, rank1 == rank2 ? "PASS" : "FAIL");
        GlobalStopWatch.Show();
        Console.WriteLine();
    }

    GlobalStopWatch.Show();
    Console.WriteLine();
    Console.Beep();
}

// void craft2()
{
    Ring.DisplayPolynomial = MonomDisplay.StarCaret;
    GlobalStopWatch.Restart();
    var seq = new[] { 1, 5, 6, 7, 14, 15, 21, 22, 29, 30, 31, 34, 78, 210, 1254, 29274 }.Select(n => n * n)
        .Concat([-1, 5]).Order().ToArray();

    foreach (var b in seq)
        EllGenerators(-b);

    GlobalStopWatch.Show(); // Time:48.407s
}
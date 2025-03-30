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

HashSet<EllPt<Rational>> SearchTorsionFree(int n, int rk, HashSet<EllPt<Rational>> intPts, int nMax = 100000)
{
    var setGensInt = new HashSet<BigInteger>();
    var b = 0;
    BigInteger n1 = 0;
    foreach (var b0 in nMax.SeqLazy(1))
    {
        b = b0;
        n1 = BigInteger.Pow(n * b0 * b0, 2);
        var ng1 = EllipticCurves.NagellLutzTorsionGroup(-n1, 0);
        var xs = intPts.Select(e => e.IsO ? e : new EllPt<Rational>(e.X * b0 * b0, e.Y * b0 * b0 * b0)).ToHashSet();
        foreach (var pt in xs.ToArray())
            AddNewIntegralPoint(ng1.E, xs, pt);

        var rks = TorsionFree(ng1.gEll, xs);

        var nMin = int.Min((int)SqrtBigInt(n1), nMax / 2);
        Console.WriteLine($"{ng1.E} ---> {ng1.E} with b={b0}");
        foreach (var (x, y2) in nMin.SeqLazy(step: -1).Concat((nMax / 2).SeqLazy(nMin))
                     .Select(x => (x: new BigInteger(x), y2: BigInteger.Pow(x, 3) - n1 * x)))
        {
            if (rks.Count == rk)
                break;

            if (xs.Any(e => !e.IsO & e.X.Num == x) || y2 <= 0)
                continue;

            var y = SqrtBigInt(y2);
            if (y2 != y * y)
                continue;

            var pt0 = new EllPt<Rational>(new(x), new(y));

            AddNewIntegralPoint(ng1.E, xs, pt0);
            AddNewIntegralPoint(ng1.E, xs, ng1.E.Invert(pt0));
            rks.Add(x);
        }

        if (rks.Count == rk)
        {
            setGensInt = rks;
            break;
        }
    }

    var (b2, b3) = (new Rational(b).Pow(2), new Rational(b).Pow(3));
    return setGensInt.Select(x => (x, y: SqrtBigInt(BigInteger.Pow(x, 3) - n1 * x)))
        .Select(e => new EllPt<Rational>(new Rational(e.x) / b2, new Rational(e.y) / b3))
        .ToHashSet();
}

void test3Rank()
{
    Ring.DisplayPolynomial = MonomDisplay.StarCaret;
    GlobalStopWatch.Restart();
    var seq = new[] { 5, 6, 7, 14, 15, 21, 22, 23, 29, 30, 31, 34, 39, 41, 78, 210, 840, 1254, 29274 };
    foreach (var n in seq)
    {
        Console.WriteLine($"Elliptic curve y^2 = x^3 - {n}^2 * x");
        var rank = EllipticCurves.EllRank(-n * n, 0); // instable
        var ng = EllipticCurves.NagellLutzTorsionGroup(-n * n, 0);
        DisplayGroup.HeadElements(ng.gEll);

        var intPts = ng.intPts.Union(ng.gEll).ToHashSet();
        foreach (var pt in ng.intPts)
            AddNewIntegralPoint(ng.E, intPts, pt);

        // congruence problem for n
        // a/b^2 solution of Y^2=X^3 - n^2*X
        // (a^3 - n^2*b^4 * a)/b^6 = (c/b^3)^2
        // a^3 - n^2*b^4 * a = c^2
        // new congruence for n1 = (n*b^2)

        GlobalStopWatch.AddLap();
        var setGens = SearchTorsionFree(n, rank, intPts);
        Console.WriteLine($"Elliptic curve y^2 = x^3 - {n}^2 * x");
        Console.WriteLine($"{ng.E} pts ord infty [{setGens.Glue()}]");
        GlobalStopWatch.Show();
        Console.WriteLine();
    }

    GlobalStopWatch.Show(); // Time:1m41s
}

// void testNumField()
{
    var (X, Y, b, m, n, e) = Ring.Polynomial(Rational.KOne(), "X", "Y", "b", "m", "n", "e")
        .Select(hi => new EPolynomial<Rational>(hi, hi.One)).Deconstruct();
    var (x, y) = (m / e.Pow(2), n / e.Pow(3));
    var P = x.Pow(3) - b * x - y.Pow(2);
    Console.WriteLine(P);
    Console.WriteLine(P.Num / e.Pow(6).Num);
}

BigInteger SquareFree(BigInteger a)
{
    if (a == 0 || a * a.Sign == 1)
        return a;

    return a.Sign * PrimesDec(a.Sign * a).Select(e => BigInteger.Pow(e.Key, e.Value % 2))
        .Aggregate((ai, aj) => ai * aj);
}

BigInteger Alpha(EllGroup<Rational> E, EllPt<Rational> pt)
{
    if (pt.IsO)
        return 0;

    if (pt.X.IsZero() && pt.Y.IsZero())
        return SquareFree(E.A.Num);

    return SquareFree(pt.X.Num);
}

void craft1()
{
    var n = 5;
    Console.WriteLine($"Elliptic curve y^2 = x^3 - {n}^2 * x");
    // var rank = EllipticCurves.EllRank(-n * n, 0); // instable
    {
        var ng = EllipticCurves.NagellLutzTorsionGroup(-n * n, 0);
        DisplayGroup.HeadElements(ng.gEll);

        var intPts = ng.intPts.Union(ng.gEll).ToHashSet();
        foreach (var pt in ng.intPts)
            AddNewIntegralPoint(ng.E, intPts, pt);

        intPts.Println("intPts");

        var alphas = intPts.Select(pt => Alpha(ng.E, pt)).Where(e => e != 0)
            .SelectMany(e => DividorsBigInt(e))
            .SelectMany(e => new[] { e, -e }).ToHashSet();
        Console.WriteLine($"Alphas:[{alphas.Order().Glue(", ")}]");

        var b1 = alphas.Select(e => (e, -n / e)).ToHashSet();
        b1.Select(e => $"N^2 = {e.Item1}*M^4 + {e.Item2}*e^4").Println("System alpha");
        Console.WriteLine();
    }

    {
        var ng = EllipticCurves.NagellLutzTorsionGroup(4 * n * n, 0);
        DisplayGroup.HeadElements(ng.gEll);

        var intPts = ng.intPts.Union(ng.gEll).ToHashSet();
        foreach (var pt in ng.intPts)
            AddNewIntegralPoint(ng.E, intPts, pt);

        intPts.Println("intPts");

        var alphas = intPts.Select(pt => Alpha(ng.E, pt)).Where(e => e != 0)
            .SelectMany(e => DividorsBigInt(e))
            .SelectMany(e => new[] { e, -e }).ToHashSet();
        Console.WriteLine($"Alphas:[{alphas.Order().Glue(", ")}]");

        var b1 = alphas.Select(e => (e, 4 * n / e)).ToHashSet();
        b1.Select(e => $"N^2 = {e.Item1}*M^4 + {e.Item2}*e^4").Println("System alpha bar");
        Console.WriteLine();
    }
}

void craft2()
{
    var nMax = 50000;
    var seq = new[] { 5, 6, 7, 14, 15, 21, 22, 23, 29, 30, 31, 34, 39, 41, 78, 210, 840, 1254, 29274 };
    foreach (var n in seq)
    {
        nMax.SeqLazy(1)
            .Select(i => (i, d: (int)double.Round(double.Sqrt(i * i * n)), d2: i * i * n))
            .Where(e => e.d2 == e.d * e.d)
            .Take(5)
            .Println($"n = {n} (d, d2)");
    }
}

HashSet<(Rational x, Rational y, Rational z)> SolveEq(int a, int b, int c, int n)
{
    var m = double.Min(a, int.Min(b, c));
    var max = (int)double.Sqrt(n / m) + 1;
    var rg = max.Range().Select(i => new Rational(i)).ToArray().Grid2D();
    return rg.Select(e => (x: e.t1, y: e.t2)).Select(e => (e.x, e.y, az2: (n - a * e.x.Pow(2) - b * e.y.Pow(2)) / c))
        .Where(e => e.az2.Sign == 1 && e.az2.IsInteger() && e.az2.IsSquare)
        .Select(e => (e.x, e.y, z: Rational.Sqrt(e.az2)))
        .ToHashSet();
}

bool TunnelCriterion(int ni)
{
    var nf = PrimesDec(ni).Select(e => e.Key.Pow(e.Value % 2)).Aggregate(1, (acc, ai) => acc * ai);
    Console.WriteLine(new { ni, nf });

    var A = SolveEq(2, 1, 32, nf);
    A.Println("A:2x^2+y^2+32z^2=n");
    var B = SolveEq(2, 1, 8, nf);
    B.Println("B:2x^2+y^2+8z^2=n");
    var C = SolveEq(8, 2, 64, nf);
    C.Println("C:8x^2+2y^2+64z^2=n");
    var D = SolveEq(8, 2, 16, nf);
    D.Println("D:8x^2+2y^2+16z^2=n");

    if (nf % 2 == 1)
        return A.SetEquals(B);

    return C.SetEquals(D);
}

{
    GlobalStopWatch.Restart();
    var set = new HashSet<int>();
    HashSet<int> A003274 =
    [
        5, 6, 7, 13, 14, 15, 20, 21, 22, 23, 24, 28, 29, 30, 31, 34, 37, 38, 39, 41, 45, 46, 47, 52, 53, 54, 55, 56, 60,
        61, 62, 63, 65, 69, 70, 71, 77, 78, 79, 80, 84, 85, 86, 87, 88, 92, 93, 94, 95, 96, 101, 102, 103, 109, 110,
        111, 112, 116, 117, 118, 119, 120, 124, 125, 126
    ];

    var (min, max) = (2, 20);
    A003274.RemoveWhere(i => i < min || i > max);
    for (int n = min; n <= max; ++n)
    {
        Console.WriteLine($"Elliptic curve y^2 = x^3 - {n}^2 * x");
        // var rank = EllipticCurves.EllRank(-n * n, 0); // instable
        if (TunnelCriterion(n))
            set.Add(n);
    }

    Console.WriteLine($"Congruent Numbers upto {max} [{set.Glue(", ")}]");
    set.Except(A003274).Println("errors");
    A003274.Except(set).Println("missing");
    GlobalStopWatch.Show();
}
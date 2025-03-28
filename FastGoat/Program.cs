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

    var setRemIntsPts = set.OrderBy(e => g.Contains(e) ? 0 : 1).ThenBy(e => e).ToHashSet();
    var setRankGens = new List<EllPt<Rational>>();
    var setIntPts = g.ToHashSet();
    setRemIntsPts.ExceptWith(g);
    while (setRemIntsPts.Count != 0)
    {
        var e0 = setRemIntsPts.MaxBy(t0 =>
        {
            var tmp = setIntPts.Union([t0, g.Invert(t0)]).ToHashSet();
            foreach (var pt in tmp.ToArray())
                AddNewIntegralPoint((EllGroup<Rational>)g.BaseGroup, tmp, pt);

            return tmp.Count;
        });

        var e1 = g.Invert(e0);
        if (!g.Contains(e0))
            setRankGens.Add(e0);

        setIntPts.UnionWith([e0, e1]);
        foreach (var pt in setIntPts.ToArray())
            AddNewIntegralPoint((EllGroup<Rational>)g.BaseGroup, setIntPts, pt);

        setRemIntsPts.ExceptWith(setIntPts);
    }
    
    return setRankGens.Select(e => e.X.Num).ToHashSet();
}

HashSet<EllPt<Rational>> SearchTorsionFree(int n, int rk, HashSet<EllPt<Rational>> intPts, int nMax = 50000)
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

        var start = -(int)SqrtBigInt(n1);
        Console.WriteLine($"E={ng1.E} b={b0}");
        foreach (var (x, y2) in nMax.SeqLazy(start)
                     .Select(x => (x: new BigInteger(x), y2: BigInteger.Pow(x, 3) - n1 * x)))
        {
            if (rks.Count == rk)
                break;

            if (xs.Any(e => !e.IsO & e.X.Num == x) || y2 <= 0)
                continue;

            var dec = PrimesDec(y2);
            if (dec.Any(e => e.Value % 2 != 0))
                continue;

            var y = dec.Select(e => BigInteger.Pow(e.Key, e.Value / 2)).Aggregate((ai, aj) => ai * aj);
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

{
    Ring.DisplayPolynomial = MonomDisplay.StarCaret;
    GlobalStopWatch.Restart();
    var seq = new[] { 5, 6, 7, 14, 15, 21, 22, 23, 29, 30, 31, 34, 39, 41, 78, 210, 840, 1254, 29274 };
    foreach (var n in seq)
    {
        Console.WriteLine($"Elliptic curve y^2 = x^3 - {n}^2 * x");
        var rank = EllipticCurves.EllRank(-n * n, 0);
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
    
    GlobalStopWatch.Show();
}
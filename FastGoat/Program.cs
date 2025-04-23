using System.Numerics;
using FastGoat.Commons;
using FastGoat.Examples;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.EllCurve;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Polynoms;
using static FastGoat.Commons.IntExt;
using RegXGroup = System.Text.RegularExpressions.Group;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");
RecomputeAllPrimesUpTo(200000);

void testEllDB()
{
    Ring.DisplayPolynomial = MonomDisplay.StarCaret;
    GlobalStopWatch.Restart();

    var ellDB = EllipticExt.LMFDB_Ell_Q().Select(e => new EllDB(e.name, e.conductor, e.rank, e.torsType, e.model))
        .ToArray();
    foreach (var e in ellDB)
    {
        Console.WriteLine(e);
        var ellAnRank = EllipticCurvesPart2.EllAnalyticRank(e.model);
        var ngl = EC.NagellLutzTorsionGroup(ellAnRank.E.ToEllGroup());
        if (ellAnRank.rank != e.rank || ellAnRank.N.Num != e.conductor || !ngl.abType.SequenceEqual(e.torsType))
            throw new($"N={ellAnRank.N} rank={ellAnRank.rank} torsType=[{ngl.abType.Glue(", ")}]");
    }

    GlobalStopWatch.Show($"EllDB {ellDB.Length} curves");
    Console.WriteLine();
}

int FindN(int p, double pmin, double pmax, int L)
{
    if (p <= L)
        return 1;

    var n0 = (int)double.Round(p * 1.0 / L);
    if (n0 >= 1 && n0 * L >= pmin && n0 * L <= pmax && L % n0 == 0 && (p - 1) % n0 == 0)
        return n0;

    var (n1, n2) = (-1, -1);
    var _n1 = (int)double.Ceiling(p * 1.0 / L);
    if (_n1 >= 1 && _n1 * L >= pmin && _n1 * L <= pmax && L % _n1 == 0 && (p - 1) % _n1 == 0)
        n1 = _n1;
    
    var _n2 = (int)double.Floor(p * 1.0 / L);
    if (_n2 >= 1 && _n2 * L >= pmin && _n2 * L <= pmax && L % _n2 == 0 && (p - 1) % _n2 == 0)
        n2 = _n2;

    return int.Max(n1, n2);
}

int EllApBSGS(EllGroup<Rational> E, int p)
{
    var d = 2 * double.Sqrt(p);
    var (pmin, pmax) = (p + 1 - d, p + 1 + d);
    Console.WriteLine(new { p, d, pmin, pmax });
    if (p < 4 * d)
        return EC.EllAp(E, p);

    GlobalStopWatch.InfiniteLoopBreakerReset();

    var (L, N) = (1, -1);
    ZnBigInt xP, yP;
    var Ep = E.ToZnBigInt(p);
    var P = Ep.O;
    while (N == -1)
    {
        while (true)
        {
            xP = new ZnBigInt(p, Rng.Next(p));
            var y2P = xP.Pow(3) + Ep.ShortForm.A * xP + Ep.ShortForm.B;
            if (LegendreJacobiBigint(y2P.Unsigned, p) != 1)
                continue;

            yP = xP.One * NumberTheory.SqrtModANTV1(y2P.K, p);
            break;
        }

        var Q = Ep.ConvertFromShort(new EllPt<ZnBigInt>(xP, yP));
        P = Ep.Op(Q, P);
        var ordP = Group.BSGS(Ep, P, Ep.O, (int)pmax);

        var M = ordP;
        foreach (var (pi, _) in PrimesDec(ordP))
        {
            while (true)
            {
                if (Ep.Times(P, M / pi).IsO)
                    M /= pi;
                else
                    break;
            }
        }

        L = Lcm(M, L);
        N = FindN(p, pmin, pmax, L);
        Console.WriteLine(new { p, d, L, N, P, ordP, M });
        GlobalStopWatch.InfiniteLoopBreaker((int)double.Pow(p, 0.25) * 10);
    }

    var ap = -N * L + p + 1;
    // ap = new[] { ap, ap - L, ap + L }.MinBy(e => int.Abs(e));

    var ap0 = EC.EllAp(E, p);
    Console.WriteLine(new { L, ap, ap0 });
    Console.WriteLine();
    if (ap != ap0)
        throw new();

    return ap;
}

{
    foreach (var e in EllipticExt.LMFDB_Ell_Q().Where(e => e.conductor < 10000))
    {
        Console.WriteLine(new EllDB(e.name, e.conductor, e.rank, e.torsType, e.model));
        var E = EC.EllGroup(e.model);
        foreach (var p in Primes10000.Where(p => BigInteger.Pow(p, 2) < 4 * e.conductor && !E.Disc.Mod(p).IsZero()))
            EllApBSGS(E, p);
        
        Console.WriteLine();
    }
}
using System.Numerics;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.EllCurve;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Polynoms;
using GFelt = FastGoat.Structures.VecSpace.EPoly<FastGoat.UserGroup.Integers.ZnInt>;

namespace FastGoat.Examples;

public static class EllipticCurvesPart2
{
    static EllipticCurvesPart2()
    {
        GlobalStopWatch.Restart();
        Ring.DisplayPolynomial = MonomDisplay.StarCaret;
    }

    static (Rational N, Dictionary<int, int> ellAn, TateAlgo[] tate, EllCoefs<Rational>E, EllGroup<Rational> Ell)
        EllInfos(BigInteger[] curve)
    {
        var E = EC.EllCoefs(curve);
        var (Ell, N, tate) = EC.EllTateAlgorithm(E);
        var ellAn = EC.EllAn(Ell, N);

        E.Show();
        Console.WriteLine($"Model {E.ModelStr}");
        Console.WriteLine($"Kodaira=[{tate.Select(e => e.kp).Glue(", ")}] Cp=[{tate.Select(e => e.cp).Glue(", ")}]");
        Console.WriteLine($"Conductor={N}");
        Console.WriteLine($"EllAn = [{ellAn.AscendingByKey().GlueMap(", ", "{0}:{1}")}]");
        Console.WriteLine();

        return (N, ellAn, tate, E, Ell);
    }

    public static (int rank, double L, Rational N, Dictionary<int, int> ellAn, TateAlgo[] tate, EllCoefs<Rational> E)
        EllAnalyticRank(BigInteger[] curve)
    {
        var E = EC.EllCoefs(curve);
        var (rank, L, N, ellAn, tate, _) = EC.AnalyticRank(E);

        E.Show();
        Console.WriteLine($"Model {E.ModelStr}");
        Console.WriteLine($"Kodaira=[{tate.Select(e => e.kp).Glue(", ")}] Cp=[{tate.Select(e => e.cp).Glue(", ")}]");
        Console.WriteLine($"Conductor={N}");
        Console.WriteLine($"EllAn = [{ellAn.Where(e => e.Key <= 20).AscendingByKey().GlueMap(", ", "{0}:{1}")}]");
        Console.WriteLine($"Analytic Rank = {rank}");
        Console.WriteLine($"L^({rank})(E, 1) {(rank == 0 ? "" : $"/ {rank}!")} = {L:F6}");
        Console.WriteLine();

        return (rank, L, N, ellAn, tate, E);
    }

    static void EllFq(BigInteger[] curve, int p)
    {
        var n = (int)(double.Log(1000) / double.Log(p));
        var E = EC.EllGroup(curve);
        for (int i = 1; i <= n; i++)
        {
            var q = p.Pow(i);
            var card = q + 1 - EC.EllAq(E, q);
            Console.WriteLine($"q={q,-6} #{$"{E.ToGF(q)}",-25} = {card}");
            var gEll = EC.EllFq(curve, q);
            var abEll = Group.AbelianDecompositions(gEll);
            var abType = abEll.abType.Select(e => e.o).Glue(" x ");
            abEll.abType.Println(l => $"{l.g} of order {l.o}", $"Generators of {gEll.ShortName} ~ [{abType}]");
            if (card != gEll.Count())
                throw new();

            Console.WriteLine();
        }
    }

    static IEnumerable<GFelt> SetFq(GFelt g)
    {
        var (p, n) = (g.P, g.F.Degree);
        if (n == 1)
            g = g.One * NumberTheory.PrimitiveRootMod(p);
        
        yield return g.Zero;
        var x = g.One;
        for (int i = 0; i < p.Pow(n) - 1; i++)
        {
            yield return x;
            x *= g;
        }
    }

    static int EllCardBF(EllGroup<Rational> E0, int q)
    {
        var card = 1;
        var g = FG.FqX(q, 'a');
        EllGroup<GFelt> E = E0.ToGF(g);
        if (int.IsPow2(q))
        {
            // Mathematics of Public Key Cryptography. Version 2.0
            // Steven D Galbraith
            // CHAPTER 2. BASIC ALGORITHMIC NUMBER THEORY
            // 2.14.2 Solving Quadratic Equations in Finite Fields
            // Exercise 2.14.7. page 67
            foreach (var x in SetFq(g))
            {
                var a = E.a1 * x + E.a3;
                var b = x.Pow(3) + E.a2 * x.Pow(2) + E.a4 * x + E.a6;
                if (a.IsZero())
                {
                    if (b.IsZero() || NumberTheory.LegendreJacobiGf(b).IsOne())
                        ++card;
                }
                else
                {
                    if (IntFactorisation.Mk(b / a.Pow(2), 1, q).IsZero())
                        card += 2;
                }
            }
        }
        else
        {
            var Eql = E.ToLongWeierstrassForm();
            foreach (var x in SetFq(g))
            {
                var resQuad = NumberTheory.LegendreJacobiGf(x.Pow(3) + Eql.a2 * x.Pow(2) + Eql.a4 * x + Eql.a6);
                if (resQuad.IsOne())
                    card += 2;
                else if (resQuad.IsZero())
                    ++card;
            }
        }

        return card;
    }

    static void EllFqCard(BigInteger[] curve, int p, bool bf = false)
    {
        var n = (int)(double.Log(100000) / double.Log(p));
        var E = EC.EllGroup(curve);
        var ellAp = EC.EllAp(E, p);
        for (int i = 1; i <= n; i++)
        {
            var q = p.Pow(i);
            var aq = EC.EllAq(ellAp, p, i);
            var card = q + 1 - aq;
            if (q > 3000 || !bf)
            {
                Console.WriteLine($"q={q,-6} #{$"{E.ToGF(q)}",-25} = {card}");
            }
            else
            {
                var cardBf = EllCardBF(E, q);
                Console.WriteLine($"q={q,-6} #{$"{E.ToGF(q)}",-25} = {card,-7} cardBf = {cardBf}");
                if (card != cardBf)
                    throw new();
            }
        }

        Console.WriteLine();
    }

    public static void Example1TateAlgorithm()
    {
        EllInfos([-1, 0]);
        EllInfos([1, 0]);
        EllInfos([-5, 0]);
        EllInfos([5, 0]);
        EllInfos([-25, 0]);
        EllInfos([-961, 0]);
        EllInfos([-3, -18]);
        EllInfos([-123, -522]);
        EllInfos([0, -1, 1, 0, 0]);
        EllInfos([0, -1, 1, -10, -20]);
        EllInfos([0, 1, 0, 16, 180]);
        EllInfos([1, 1, 0, -22, -44]);
        EllInfos([1, -1, 1, -180, 1047]);
        EllInfos([0, -1, 1, 444, -826]);
        EllInfos([0, 1, 0, -5, 7]);
        EllInfos([1, -1, 0, -18, -81]);
        EllInfos([1, -1, 0, -17, 16]);
    }

    public static void Example2AnalyticRank()
    {
        GlobalStopWatch.AddLap();

        // rank 0
        EllAnalyticRank([-1, 0]);
        EllAnalyticRank([1, 0]);
        EllAnalyticRank([0, -1, 1, 0, 0]);
        EllAnalyticRank([0, -1, 1, -10, -20]);
        EllAnalyticRank([-3, -18]);
        EllAnalyticRank([-123, -522]);
        EllAnalyticRank([0, 1, 0, 16, 180]);
        EllAnalyticRank([0, 1, 1, 1, 0]);
        EllAnalyticRank([7, 0, 0, 16, 0]);
        EllAnalyticRank([1, -5, -5, 0, 0]);
        EllAnalyticRank([0, -1, -1, 0, 0]);
        EllAnalyticRank([1, -1, 1, -14, 29]);
        EllAnalyticRank([1, 0, 0, -45, 81]);
        EllAnalyticRank([43, -210, -210, 0, 0]);
        EllAnalyticRank([5, -3, -6, 0, 0]);
        EllAnalyticRank([17, -60, -120, 0, 0]);
        EllAnalyticRank([1, 0, 0, -1070, 7812]);
        EllAnalyticRank([-4, -4]);
        EllAnalyticRank([1, -1, 1, -29, -53]);
        EllAnalyticRank([-11, -14]);
        EllAnalyticRank([1, 1, 0, 1, 0]);
        EllAnalyticRank([1, 0, 1, -14, -64]);
        EllAnalyticRank([0, 4]);
        EllAnalyticRank([0, 1, 1, 1, -1]);
        EllAnalyticRank([1, 0, 0, 6, -28]);
        EllAnalyticRank([-7, -6]);
        EllAnalyticRank([-2, 1]);
        EllAnalyticRank([1, 1, 1, 0, 0]);
        EllAnalyticRank([1, -1, 1, -6, -4]);
        EllAnalyticRank([1, 0, 0, 15, 9]);
        EllAnalyticRank([1, 1, 1, 0, 1]);
        EllAnalyticRank([0, 1]);
        EllAnalyticRank([1, -1, 0, 6, 0]);
        EllAnalyticRank([1, 0, 1, -6, 4]);
        EllAnalyticRank([1, 0, 0, 159, 1737]);
        EllAnalyticRank([1, 0, 0, -1, 137]);
        EllAnalyticRank([1, -1, 1, -3, 3]);
        EllAnalyticRank([0, -1, 0, -4, 4]);
        EllAnalyticRank([1, 0, 0, -34, 68]);
        EllAnalyticRank([1, 0, 0, -4, -1]);
        EllAnalyticRank([1, 0, 0, 108, 11664]);
        EllAnalyticRank([1, 0, 0, 115, 561]);
        EllAnalyticRank([1, 0, 0, -828, 9072]);
        EllAnalyticRank([1, 0, 1, -19, 26]);
        EllAnalyticRank([1, -1, 1, -122, 1721]);
        EllAnalyticRank([1, 0, 0, -361, 2585]);
        EllAnalyticRank([1, 0, 1, 1922, 20756]);
        EllAnalyticRank([1, 0, 0, -1070, 7812]);
        EllAnalyticRank([1, 0, 0, -8696090, 9838496100]);
        // EllAnalyticRank([238, 952]); // slow, Conductor N=627162368

        // rank 1
        EllAnalyticRank([-5, 0]);
        EllAnalyticRank([-7, 0]);
        EllAnalyticRank([-25, 0]);
        EllAnalyticRank([-961, 0]);
        EllAnalyticRank([0, 0, 1, -1, 0]);
        EllAnalyticRank([0, 4, 0, -80, 400]);
        EllAnalyticRank([0, 0, 1, 1, 0]);
        EllAnalyticRank([1, 0, 0, -4767, 127449]);

        // rank 2
        EllAnalyticRank([-34 * 34, 0]);
        EllAnalyticRank([-41 * 41, 0]);
        EllAnalyticRank([-25, 25]);
        EllAnalyticRank([-81, 81]);

        // rank 3
        EllAnalyticRank([1, 1, 0, -36, 36]);
        EllAnalyticRank([1, 1, 0, -49, 49]);
        EllAnalyticRank([0, 0, 1, -7, 6]);

        // rank 4
        EllAnalyticRank([1, 1, 0, -87, 225]);
        EllAnalyticRank([1, 1, 0, -104, 276]);
        EllAnalyticRank([1, -1, 0, -79, 289]);

        // rank 5
        EllAnalyticRank([0, 0, 1, -79, 342]);

        GlobalStopWatch.Show();
        Console.WriteLine();
    }

    public static void Example3LMFDB()
    {
        GlobalStopWatch.AddLap();

        var ellDB = EllipticExt.LMFDB_Ell_Q().Select(e => new EllDB(e.name, e.conductor, e.rank, e.torsType, e.model))
            .ToArray();
        foreach (var e in ellDB)
        {
            Console.WriteLine(e);
            var ellAnRank = EllAnalyticRank(e.model);
            var ngl = EC.NagellLutzTorsionGroup(ellAnRank.E.ToEllGroup());
            if (ellAnRank.rank != e.rank || ellAnRank.N.Num != e.conductor || !ngl.abType.SequenceEqual(e.torsType))
                throw new($"N={ellAnRank.N} rank={ellAnRank.rank} torsType=[{ngl.abType.Glue(", ")}]");
        }

        GlobalStopWatch.Show($"EllDB {ellDB.Length} curves"); // # EllDB 166 curves Time:1m4s
        Console.WriteLine();
    }

    public static void Example4LargePrime()
    {
        var p1 = 10.Pow(6) + 3;
        var E = EC.EllGroup([1, 1]);
        EC.EllApBSGS(E, p1, LogLevel.Level1);

        var p2 = BigInteger.Pow(10, 10) + 19;
        EC.EllApBSGS(E, p2, LogLevel.Level1);

        GlobalStopWatch.Bench(5, $"BSGS p={p1}", () => EC.EllApBSGS(E, p1));
        GlobalStopWatch.Bench(5, $"BF   p={p1}", () => EC.EllAp(E, p1));

        GlobalStopWatch.Restart();
        EC.EllApBSGS(E, p2);
        GlobalStopWatch.Show($"BSGS p={p2}");
        /*
           Ep = Ell[1,1](Z/1000003Z) Ap = -723 Case = [1]#1
           Ep = Ell[1,1](Z/10000000019Z) Ap = 118240 Case = [1]#1
           # BSGS p=1000003 Avg Time:44 ms Dev:4.454
           # BF   p=1000003 Avg Time:595 ms Dev:10.423
           # BSGS p=10000000019 Time:6.573s
         */
    }

    public static void Example5CountingPointsFp()
    {
        var pmin = 64;
        foreach (var e in EllipticExt.LMFDB_Ell_Q().Where(e => e.conductor < 100000))
        {
            var E = EC.EllGroup(e.model);
            var pmax = (int)(2 * double.Sqrt((double)e.conductor));
            foreach (var p in IntExt.Primes10000.Where(p => p > pmin && p < pmax && e.conductor % p != 0))
            {
                Console.WriteLine(new EllDB(e.name, e.conductor, e.rank, e.torsType, e.model));

                var actualAp = EC.EllApBSGS(E, p, LogLevel.Level1);
                var expectedAp = EC.EllAp(E, p);
                if (actualAp != expectedAp)
                    throw new($"p={p} ap={actualAp} expected={expectedAp}");

                Console.WriteLine();
            }
        }
    }

    public static void Example6CountingPointsFp()
    {
        var E = EC.EllGroup([238, 952]);
        foreach (var p in IntExt.Primes10000.Where(p => p > 20000 && p < 26000))
        {
            var actualAp = EC.EllApBSGS(E, p, LogLevel.Level1);
            var expectedAp = EC.EllAp(E, p);
            if (actualAp != expectedAp)
                throw new($"p={p} ap={actualAp} expected={expectedAp}");
        }
    }

    public static void Example7EllFqCharacteristic_2_3()
    {
        BigInteger[][] curves =
        [
            [0, 1, 1, 0, 1],
            [0, 0, 1, 0, 1],
            [1, 0, 1, 0, 1],
            [1, 1, 0, 0, 1],
            [1, 0, 0, 1, 0],
            [1, 0, 2, 0, 1],
            [2, 0, 0, 1, 0]
        ];
        
        {
            var p = 2;
            foreach (var curve in curves.Where(e => !EC.EllCoefs(e).ToZnInt(p).Disc.IsZero()))
                EllFq(curve, p);
        }

        {
            var p = 3;
            foreach (var curve in curves.Where(e => !EC.EllCoefs(e).ToZnInt(p).Disc.IsZero()))
                EllFq(curve, p);
        }

        {
            var p = 7;
            foreach (var curve in curves.Where(e => !EC.EllCoefs(e).ToZnInt(p).Disc.IsZero()))
                EllFq(curve, p);
        }

        {
            var p = 17;
            foreach (var curve in curves.Where(e => !EC.EllCoefs(e).ToZnInt(p).Disc.IsZero()))
                EllFq(curve, p);
        }

    }

    public static void Example8CountingPointsFq()
    {
        BigInteger[][] curves =
        [
            [0, 1, 1, 0, 1],
            [0, 0, 1, 0, 1],
            [1, 0, 1, 0, 1],
            [1, 1, 0, 0, 1],
            [1, 0, 0, 1, 0],
            [1, 0, 2, 0, 1],
            [2, 0, 0, 1, 0]
        ];

        foreach (var p in IntExt.Primes10000.Take(10))
        {
            foreach (var curve in curves.Where(e => !EC.EllCoefs(e).ToZnInt(p).Disc.IsZero()))
                EllFqCard(curve, p, bf: true);
        }
    }

    public static void Example9EllApFrobTrace()
    {
        GlobalStopWatch.Restart();
        EC.EllApFrobTrace([1, 0]);
        EC.EllApFrobTrace([-1, 0]);
        EC.EllApFrobTrace([1, 1]);
        EC.EllApFrobTrace([-43, 166]);
        EC.EllApFrobTrace([0, 0, 1, 0, -7]);
        EC.EllApFrobTrace([1, 0, 0, 0, 1]);
        EC.EllApFrobTrace([1, -1, 0, -4, 4]);
        EC.EllApFrobTrace([1, -1, 1, -19353, 958713]);
        EC.EllApFrobTrace([1, 1, 1, -17714, 900047]);
        GlobalStopWatch.Show();
    }
}
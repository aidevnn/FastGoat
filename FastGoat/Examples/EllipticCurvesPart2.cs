using System.Numerics;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.EllCurve;
using FastGoat.UserGroup.Integers;

namespace FastGoat.Examples;

public static class EllipticCurvesPart2
{
    static (Rational N, Dictionary<int, int> ellAn, TateAlgo[] tate, EllCoefs<Rational>E, EllGroup<Rational> Ell)
        EllInfos(BigInteger[] curve)
    {
        var (a1, a2, a3, a4, a6) = EC.CurveArray(curve.Select(i => new Rational(i)).ToArray()).Deconstruct();
        var E = new EllCoefs<Rational>(a1, a2, a3, a4, a6);
        var Ell = new EllGroup<Rational>(a1, a2, a3, a4, a6);
        var dec = IntExt.PrimesDec(E.Disc.Absolute.Num);
        var tate = dec.Keys.Select(p => EC.TateAlgorithm(E, p)).ToArray();
        var N = tate.Select(e => new Rational(e.p).Pow(e.fp)).Aggregate((pi, pj) => pi * pj);
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
        var (a1, a2, a3, a4, a6) = EC.CurveArray(curve.Select(i => new Rational(i)).ToArray()).Deconstruct();
        var E = new EllCoefs<Rational>(a1, a2, a3, a4, a6);
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

    public static void Example1TateAlgorithm()
    {
        Ring.DisplayPolynomial = MonomDisplay.StarCaret;

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
        GlobalStopWatch.Restart();
        Ring.DisplayPolynomial = MonomDisplay.StarCaret;

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
}
using System.Numerics;
using System.Text.RegularExpressions;
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
    GlobalStopWatch.AddLap();

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

List<BigInteger> CandidatsN(double pmin, double pmax, BigInteger L)
{
    var n0 = (int)(pmin / (double)L);
    var n1 = (int)(pmax / (double)L) + 1;
    return (n1 - n0 + 1).SeqLazy(n0).Where(e => e > 0 && L % e == 0)
        .Select(e => e * L).Where(n => (double)n > pmin && (double)n < pmax).ToList();
}

EllPt<ZnBigInt> RandPt(EllGroup<ZnBigInt> Ep, HashSet<BigInteger> set)
{
    var p = Ep.ShortForm.A.Mod;
    while (true)
    {
        var xP = new ZnBigInt(p, DistributionExt.Dice(BigInteger.Zero, p));
        var y2P = xP.Pow(3) + Ep.ShortForm.A * xP + Ep.ShortForm.B;
        if (LegendreJacobiBigint(y2P.Unsigned, p) != 1)
            continue;

        var yP = xP.One * NumberTheory.SqrtModANTV1(y2P.K, p);
        if (set.Add(xP.Unsigned))
            return Ep.ConvertFromShort(new(xP, yP));
    }
}

long BSGSlong<T>(IGroup<T> g, T a, T b, double ord) where T : struct, IElt<T>
{
    var (m, tmp1) = ((long)Double.Sqrt(ord) + 1, a); // TODO: faster BSGS
    var L = new Dictionary<T, long>() { [a] = 1 };
    for (long i = 1; i < m; i++)
    {
        if (tmp1.Equals(b))
            return i;

        tmp1 = g.Op(tmp1, a);
        L[tmp1] = i + 1;
    }

    if (tmp1.Equals(b))
        return m;

    var (c, tmp2) = (g.Invert(tmp1), b);
    for (long j = 1; j < m; j++)
    {
        tmp2 = g.Op(tmp2, c);
        if (L.TryGetValue(tmp2, out long i))
            return j * m + i;
    }

    throw new($"{g.Name} ord={ord}; a={a} b={b}");
}

BigInteger EllApBSGS(EllGroup<Rational> E, BigInteger p, LogLevel log = LogLevel.Level2)
{
    var d = double.Sqrt((double)p);
    var (pmin, pmax) = ((double)p + 1 - 2 * d, (double)p + 1 + 2 * d);

    (BigInteger L, BigInteger N) = (1, -1);
    (BigInteger L_, BigInteger N_) = (1, -1);

    var g = 1000.SeqLazy().Select(_ => DistributionExt.Dice(BigInteger.One, p - 1))
        .Where(g => LegendreJacobiBigint(g, p) != 1).Select(g => new ZnBigInt(p, g)).First();
    var Ep = E.ToZnBigInt(p);
    // Mestre's Theorem, quadratic twist of Ep
    var Ep_ = new EllGroup<ZnBigInt>(Ep.ShortForm.A * g.Pow(2), Ep.ShortForm.B * g.Pow(3));
    int ct = 0, nbCands = 0, nbCands_ = 0;
    var set = new HashSet<BigInteger>();
    var set_ = new HashSet<BigInteger>();
    var cands = new List<BigInteger>();
    var cands_ = new List<BigInteger>();
    while (ct < 2 * d)
    {
        ++ct;
        var P = RandPt(Ep, set);
        var ordP = BSGSlong(Ep, P, Ep.O, pmax);

        L = LcmBigInt(ordP, L);
        cands = CandidatsN(pmin, pmax, L);
        if (cands.Count >= 1)
            N = cands[0] / L;

        var NL = N * L;
        var arrCands = $"[{cands.Glue(", ")}]";
        nbCands = cands.Count;

        var P_ = RandPt(Ep_, set_);
        var ordP_ = BSGSlong(Ep_, P_, Ep_.O, pmax);
        L_ = LcmBigInt(ordP_, L_);
        cands_ = CandidatsN(pmin, pmax, L_);
        if (cands_.Count >= 1)
            N_ = cands_[0] / L_;

        var NL_ = N_ * L_;
        var arrCands_ = $"[{cands_.Glue(", ")}]";
        nbCands_ = cands_.Count;

        if (log == LogLevel.Level2)
        {
            Console.WriteLine(new { p, d, pmin, pmax });
            Console.WriteLine(new { NL, N, L, arrCands, nbCands });
            Console.WriteLine(new { NL_, N_, L_, arrCands_, nbCands_ });
        }

        if (N == 1 || N_ == 1)
            break;

        if (nbCands == 1 && L % N == 0 && (p - 1) % N == 0 && cands_.Contains(2 * p + 2 - NL))
        {
            var N0 = (2 * p + 2 - NL) / L_;
            if (L_ % N0 == 0 && (p - 1) % N0 == 0)
                break;
        }

        if (nbCands_ == 1 && L_ % N_ == 0 && (p - 1) % N_ == 0 && cands.Contains(2 * p + 2 - NL_))
        {
            var N0_ = (2 * p + 2 - NL_) / L;
            if (L % N0_ == 0 && (p - 1) % N0_ == 0)
                break;
        }
    }

    int caseNum, outputNum;
    if (N == 1)
        (caseNum, outputNum) = (1, 1);
    else if (N_ == 1)
        (caseNum, outputNum) = (2, 2);
    else if (L % N == 0 && (p - 1) % N == 0 && nbCands == 1 && L_ < L)
        (caseNum, outputNum) = (3, 1);
    else if (L_ % N_ == 0 && (p - 1) % N_ == 0 && nbCands_ == 1 && L < L_)
        (caseNum, outputNum) = (4, 2);
    else if (nbCands == 1 && nbCands_ == 1 && 2 * (p + 1) == cands[0] + cands_[0])
        (caseNum, outputNum) = (5, 3);
    else if (nbCands == 2 && nbCands_ == 2)
    {
        if (L_ < L)
            (caseNum, outputNum) = (6, 1);
        else if (L < L_)
            (caseNum, outputNum) = (7, 2);
        else
            (caseNum, outputNum) = (8, 2);
    }
    else
        (caseNum, outputNum) = (9, 4);

    Stats.Update(caseNum, ct, Ep);
    if (outputNum == 1 || outputNum == 3)
    {
        var Ap = cands.Where(NL => (p - 1) % (NL / L) == 0).Select(NL => p + 1 - NL).First();
        if (log != LogLevel.Off)
            Console.WriteLine($"Ep = {Ep.ShortFormStr} Ap = {Ap} Case = [{caseNum}]#{outputNum}");

        return Ap;
    }
    else if (outputNum == 2)
    {
        var Ap = cands_.Where(NL => (p - 1) % (NL / L_) == 0).Select(NL => -(p + 1 - NL)).First();
        if (log != LogLevel.Off)
            Console.WriteLine($"Ep = {Ep.ShortFormStr} Ap = {Ap} Case = [{caseNum}]#2");

        return Ap;
    }
    else
        throw new();
}

void testEllApBSGS(int pmin = 64, LogLevel log = LogLevel.Level2)
{
    foreach (var e in EllipticExt.LMFDB_Ell_Q().Where(e => e.conductor < 100000))
    {
        var E = EC.EllGroup(e.model);
        var pmax = (int)(2 * double.Sqrt((double)e.conductor));
        foreach (var p in Primes10000.Where(p => p > pmin && p < pmax && e.conductor % p != 0))
        {
            if (log != LogLevel.Off)
                Console.WriteLine(new EllDB(e.name, e.conductor, e.rank, e.torsType, e.model));

            var actualAp = EllApBSGS(E, p, log);
            var expectedAp = EC.EllAp(E, p);
            if (actualAp != expectedAp)
                throw new($"p={p} ap={actualAp} expected={expectedAp}");

            if (log != LogLevel.Off)
                Console.WriteLine();
        }
    }
}

void testCase1()
{
    var p = 10.Pow(6) + 3;
    var E = EC.EllGroup([1, 1]);

    GlobalStopWatch.Bench(5, "BSGS", () => EllApBSGS(E, p));
    GlobalStopWatch.Bench(5, "BF", () => EC.EllAp(E, p));
    GlobalStopWatch.Bench(5, "BSGS", () => EllApBSGS(E, p, LogLevel.Off));
    GlobalStopWatch.Bench(5, "BF", () => EC.EllAp(E, p));
    GlobalStopWatch.Bench(50, "BSGS", () => EllApBSGS(E, p, LogLevel.Off)); // # BSGS Avg Time:44 ms Dev:6.730
    GlobalStopWatch.Bench(50, "BF", () => EC.EllAp(E, p)); // # BF Avg Time:589 ms Dev:10.685

    GlobalStopWatch.Restart();
    EllApBSGS(EC.EllGroup([1, 1]), BigInteger.Pow(10, 10) + 19);
    // Ep = Ell[1,1](Z/10000000019Z) Ap = 118240
    // #  Time:6.303s
    GlobalStopWatch.Show();
}

void testCase2()
{
    var E = EC.EllGroup([238, 952]);
    foreach (var p in Primes10000.Where(p => p > 20000 && p < 26000))
    {
        var actualAp = EllApBSGS(E, p, LogLevel.Level1);
        var expectedAp = EC.EllAp(E, p);
        if (actualAp != expectedAp)
            throw new($"p={p} ap={actualAp} expected={expectedAp}");
    }
}

void testBatch(int nb = 10)
{
    for (int k = 0; k < nb; ++k)
    {
        GlobalStopWatch.AddLap();
        testEllApBSGS(log: LogLevel.Level1);
        GlobalStopWatch.Show();
    }

    Console.Beep();
}

{
    // RngSeed(1251);
    GlobalStopWatch.Restart();
    // testEllApBSGS();

    // Stats.Run(testCase1);
    // Stats.Run(testCase2);
    // Stats.Run(() => testEllApBSGS(pmin: 64, log: LogLevel.Off));
    // Stats.Run(() => testEllApBSGS(pmin: 32, log: LogLevel.Off));
    Stats.Run(() => testBatch());
    // Stats.Run(() => testBatch(nb: 100));

    GlobalStopWatch.Show();
}

public static class Stats
{
    static Stats()
    {
        CasesCounter = 9.SeqLazy(1).ToDictionary(i => i, _ => 0);
        LoopCounter = new();
    }

    private static Dictionary<int, int> CasesCounter { get; }
    private static Dictionary<int, (int nb, EllGroup<ZnBigInt> E)> LoopCounter { get; }

    public static void Update(int caseNum, int ct, EllGroup<ZnBigInt> E)
    {
        CasesCounter[caseNum]++;

        if (LoopCounter.TryGetValue(ct, out var e))
            LoopCounter[ct] = (e.nb + 1, e.E);
        else
            LoopCounter[ct] = (1, E);
    }

    public static void Clear()
    {
        9.SeqLazy(1).ToList().ForEach(i => CasesCounter[i] = 0);
        LoopCounter.Clear();
    }

    public static void Show()
    {
        var total = CasesCounter.Values.Sum();
        CasesCounter.Println(e => $"Case[{e.Key:00}] = {e.Value,10}", $"Stats. Total = {total,10}");
        LoopCounter.AscendingByKey().Println(e => $"Ct = {e.Key,3} Nb = {e.Value.nb,7} First {e.Value.E}", $"Max loop");
        Console.WriteLine();
    }

    public static void Run(Action act)
    {
        Clear();
        act();
        Show();
    }
}
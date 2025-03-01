using System.Numerics;
using FastGoat.Commons;
using FastGoat.Examples;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.LWE;
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
// RngSeed(259128);
Ring.DisplayPolynomial = MonomDisplay.StarCaret;
RecomputeAllPrimesUpTo(1000000);

KMatrix<ZnBInt> Vandermonde(int n, ZnBInt w)
{
    var V = Ring.Matrix(w.One, n, n);
    var wPow = n.SeqLazy().Select(i => w.Pow(i)).ToArray();

    for (int i = 1; i < n; i++)
    for (int j = i; j < n; j++)
        V[i, j] = V[j, i] = wPow[j] * V[i - 1, j];

    return new(V);
}

Rq RqMulFFT(Rq a, Rq b, KMatrix<ZnBInt> Vw, KMatrix<ZnBInt> InvVw)
{
    var n = Vw.N;
    var A = a.CoefsExtended(n - 1).Select(e => e.Num * Vw.KOne).ToKMatrix(n);
    var _A = Vw * A;
    var B = b.CoefsExtended(n - 1).Select(e => e.Num * Vw.KOne).ToKMatrix(n);
    var _B = Vw * B;

    var _AB = _A.Zip(_B).Select(e => e.First * e.Second).ToKMatrix(n);
    return (InvVw * _AB).Select(e => new Rational(e.ToSignedBigInt)).ToKPoly();
}

KMatrix<ZnBInt> NTT(int n, ZnBInt w)
{
    var ntt = Ring.Matrix(w.One, n, n);
    var wPow = (2 * n).SeqLazy().Select(i => w.Pow(i)).ToArray();
    for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
        ntt[i, j] = wPow[(2 * i * j + j) % (2 * n)];

    return new(ntt);
}

Rq RqAddNTT(Rq a, Rq b, KMatrix<ZnBInt> ntt, KMatrix<ZnBInt> intt)
{
    var n = ntt.N;
    var A = a.CoefsExtended(n - 1).Select(e => e.Num * ntt.KOne).ToKMatrix(n);
    var _A = ntt * A;
    var B = b.CoefsExtended(n - 1).Select(e => e.Num * ntt.KOne).ToKMatrix(n);
    var _B = ntt * B;

    var _AB = _A.Zip(_B).Select(e => e.First + e.Second).ToKMatrix(n);
    return (intt * _AB).Select(e => new Rational(e.ToSignedBigInt)).ToKPoly();
}

Rq RqMulNTT(Rq a, Rq b, KMatrix<ZnBInt> ntt, KMatrix<ZnBInt> intt)
{
    var n = ntt.N;
    var A = a.CoefsExtended(n - 1).Select(e => e.Num * ntt.KOne).ToKMatrix(n);
    var _A = ntt * A;
    var B = b.CoefsExtended(n - 1).Select(e => e.Num * ntt.KOne).ToKMatrix(n);
    var _B = ntt * B;

    var _AB = _A.Zip(_B).Select(e => e.First * e.Second).ToKMatrix(n);
    return (intt * _AB).Select(e => new Rational(e.ToSignedBigInt)).ToKPoly();
}

void testTonelliShanks()
{
    {
        var (n, p) = (5, 41);
        var r = NumberTheory.SqrtModWP(n, p);
        var r2 = AmodPbigint(r * r, p);
        Console.WriteLine(new { n, p, r, r2 });
    }

    {
        var (n, p) = (1, 41);
        var r = NumberTheory.SqrtModWP(n, p);
        var r2 = AmodPbigint(r * r, p);
        Console.WriteLine(new { n, p, r, r2 });
    }

    {
        var (n, p) = (4, 241);
        Console.WriteLine($"{n}-thRoots(-1) mod {p}.");
        var seqTS = NumberTheory.Pow2NthRootsWP(-1, n, p).Order().ToArray();
        var seqAntv1 = NumberTheory.NthRootsANTV1(-1, n, p).Order().ToArray();
        var seqPM = (p - 2).SeqLazy(2).Where(i => PowMod(i, n, p) == p - 1).Order().ToArray();
        Console.WriteLine($"[{seqTS.Glue(", ")}]");
        Console.WriteLine($"[{seqAntv1.Glue(", ")}]");
        Console.WriteLine($"[{seqPM.Glue(", ")}]");
        Console.WriteLine();
    }

    {
        var (n, p) = (8, 241);
        Console.WriteLine($"{n}-thRoots(-1) mod {p}.");
        var seqTS = NumberTheory.Pow2NthRootsWP(-1, n, p).Order().ToArray();
        var seqAntv1 = NumberTheory.NthRootsANTV1(-1, n, p).Order().ToArray();
        var seqPM = (p - 2).SeqLazy(2).Where(i => PowMod(i, n, p) == p - 1).Order().ToArray();
        Console.WriteLine($"[{seqTS.Glue(", ")}]");
        Console.WriteLine($"[{seqAntv1.Glue(", ")}]");
        Console.WriteLine($"[{seqPM.Glue(", ")}]");
        Console.WriteLine();
    }

    {
        var N = 2048;
        var n = N / 2;
        var (primes0, sp) = RLWE.SequencePrimesBGV(N, level: 4);
        var primes = primes0.Append(sp).Select(e => e.Num).ToArray();
        foreach (var p in primes)
        {
            var log2p = BigInteger.Log2(p);
            Console.WriteLine(new { N, n, p, log2p });
            var seqTS = NumberTheory.Pow2NthRootsWP(-1, n, p).ToArray();
            if (seqTS.Length != n || seqTS.Any(e => !PowModEqualOnebigint(e, N, p)))
                throw new("Fail");

            Console.WriteLine($"{n}-thRoots(-1) mod {p}. Pass");
            GlobalStopWatch.Bench(nb: 50, "TS", () => NumberTheory.Pow2NthRootsWP(-1, n, p).ToArray());
            Console.WriteLine();
        }
        // { N = 2048, n = 1024, p = 23026821121, log2p = 34 }
        // 1024-thRoots(-1) mod 23026821121. Pass
        // # TS Avg Time:4 ms Dev:0.681
    }
}

void testFFT()
{
    var (N, p) = (16, 241);
    var n = N / 2;
    var w = new ZnBInt(p, NumberTheory.NthRootsANTV1(p - 1, n, p).First());
    var Ni = (N * w.One).Inv();
    Console.WriteLine($"{N}-thRoots of unity mod {p} is {w}");
    var Vw = Vandermonde(N, w);
    var InvVw = Ni * Vandermonde(N, w.Inv());
    Console.WriteLine("V(w)");
    Console.WriteLine(Vw);
    Console.WriteLine("V(iw)");
    Console.WriteLine(InvVw);
    Console.WriteLine("V * InvVw");
    Console.WriteLine(Vw * InvVw);
    Console.WriteLine();

    var m = N.SeqLazy().Select(_ => Rng.Next(p) * w.One).ToKMatrix(N);
    var _m = Vw * m;
    var __m = InvVw * _m;
    Console.WriteLine("m");
    Console.WriteLine(m.T);
    Console.WriteLine("_m");
    Console.WriteLine(_m.T);
    Console.WriteLine("__m");
    Console.WriteLine(__m.T);
    Console.WriteLine();

    var pm = FG.QPoly().Pow(n) + 1;
    var a = RLWE.GenUnif(n, p);
    var b = RLWE.GenUnif(n, p);
    var ab1 = (a * b).ResModSigned(pm, p);
    var ab2 = RqMulFFT(a, b, Vw, InvVw).ResModSigned(pm, p);
    Console.WriteLine($"a     = {a}");
    Console.WriteLine($"b     = {b}");
    Console.WriteLine($"a * b = {ab1}");
    Console.WriteLine($"      = {ab2}");
    Console.WriteLine();
}

void testNTT()
{
    var (N, p) = (16, 241);
    var n = N / 2;
    var w = new ZnBInt(p, NumberTheory.NthRootsANTV1(p - 1, n, p).First());
    var ni = (n * w.One).Inv();
    Console.WriteLine($"{n}-thRoots(-1) = {w} mod {p}");
    var ntt = NTT(n, w);
    var intt = ni * NTT(n, w.Inv()).T;
    Console.WriteLine("NTT(w)");
    Console.WriteLine(ntt);
    Console.WriteLine("INTT(w)");
    Console.WriteLine(intt);
    Console.WriteLine("NTT(w) * INTT(w)");
    Console.WriteLine(ntt * intt);
    Console.WriteLine();

    var m = n.SeqLazy(1).Select(_ => Rng.Next(p) * w.One).ToKMatrix(n);
    var _m = ntt * m;
    var __m = intt * _m;
    Console.WriteLine("m");
    Console.WriteLine(m.T);
    Console.WriteLine("_m");
    Console.WriteLine(_m.T);
    Console.WriteLine("__m");
    Console.WriteLine(__m.T);
    Console.WriteLine();

    var pm = FG.QPoly().Pow(n) + 1;

    for (int i = 0; i < 10; i++)
    {
        var a = RLWE.GenUnif(n, p);
        var b = RLWE.GenUnif(n, p);

        var add_ab1 = (a + b).CoefsModSigned(p);
        var add_ab2 = RqAddNTT(a, b, ntt, intt);
        var mul_ab1 = (a * b).ResModSigned(pm, p);
        var mul_ab2 = RqMulNTT(a, b, ntt, intt);

        Console.WriteLine($"a     = {a}");
        Console.WriteLine($"b     = {b}");
        Console.WriteLine($"a + b = {add_ab1}");
        Console.WriteLine($"      = {add_ab2}");
        if (!add_ab1.Equals(add_ab2))
            throw new();

        Console.WriteLine($"a * b = {mul_ab1}");
        Console.WriteLine($"      = {mul_ab2}");
        Console.WriteLine();
        if (!mul_ab1.Equals(mul_ab2))
            throw new();
    }
}

void testBenchNTT()
{
    var (N, p) = (1024, 12289);
    var n = N / 2;
    var w = new ZnBInt(p, NumberTheory.NthRootsANTV1(p - 1, n, p).First());
    var ni = (n * w.One).Inv();
    var ntt = NTT(n, w);
    var intt = ni * NTT(n, w.Inv()).T;
    var pm = FG.QPoly().Pow(n) + 1;
    var a = RLWE.GenUnif(n, p);
    var b = RLWE.GenUnif(n, p);

    GlobalStopWatch.Bench(5, "PolAdd", () => (a + b).ResModSigned(pm, p));
    GlobalStopWatch.Bench(5, "NTTadd", () => RqAddNTT(a, b, ntt, intt));
    GlobalStopWatch.Bench(5, "PolMul", () => (a * b).ResModSigned(pm, p));
    GlobalStopWatch.Bench(5, "NTT   ", () => RqMulNTT(a, b, ntt, intt));
    GlobalStopWatch.Bench(50, "PolAdd", () => (a + b).ResModSigned(pm, p));
    GlobalStopWatch.Bench(50, "NTTadd", () => RqAddNTT(a, b, ntt, intt));
    GlobalStopWatch.Bench(50, "PolMul", () => (a * b).ResModSigned(pm, p));
    GlobalStopWatch.Bench(50, "NTT   ", () => RqMulNTT(a, b, ntt, intt));
}

void testSqrtMod()
{
    var p = 12289;
    NumberTheory.PohligHellmanInfos? ph = NumberTheory.PreparePohligHellman(p);
    Console.WriteLine($"Sqrt mod {p}");
    var seq = p.SeqLazy().Where(m => LegendreJacobi(m, p) == 1).ToArray();

    var seqWP = seq.ToDictionary(m => m, m => NumberTheory.SqrtModWP(m, p))
        .ToDictionary(e => e.Key, e => int.Min(e.Value, p - e.Value));
    var seqANTV1 = seq.ToDictionary(m => m, m => NumberTheory.SqrtModANTV1(m, p))
        .ToDictionary(e => e.Key, e => int.Min(e.Value, p - e.Value));
    var seqMPKCV2 = seq.ToDictionary(m => m, m => NumberTheory.SqrtModMPKCV2(m, p))
        .ToDictionary(e => e.Key, e => int.Min(e.Value, p - e.Value));
    var seqMPKCV2ph = seq.ToDictionary(m => m, m => NumberTheory.SqrtModMPKCV2(m, p, ph))
        .ToDictionary(e => e.Key, e => int.Min(e.Value, p - e.Value));

    if (seqWP.Any(e => e.Value != seqANTV1[e.Key] || e.Value != seqMPKCV2[e.Key] || e.Value != seqMPKCV2ph[e.Key]))
        throw new();

    GlobalStopWatch.Bench(50, "WP             ", () => seq.Select(m => NumberTheory.SqrtModWP(m, p)).ToArray());
    GlobalStopWatch.Bench(50, "ANTV1          ", () => seq.Select(m => NumberTheory.SqrtModANTV1(m, p)).ToArray());
    GlobalStopWatch.Bench(5, "MPKCV2         ", () => seq.Select(m => NumberTheory.SqrtModMPKCV2(m, p)).ToArray());
    GlobalStopWatch.Bench(50, "MPKCV2 discrLog", () => seq.Select(m => NumberTheory.SqrtModMPKCV2(m, p, ph)).ToArray());
    GlobalStopWatch.Bench(50, "WP             ", () => seq.Select(m => NumberTheory.SqrtModWP(m, p)).ToArray());
    GlobalStopWatch.Bench(50, "ANTV1          ", () => seq.Select(m => NumberTheory.SqrtModANTV1(m, p)).ToArray());
    GlobalStopWatch.Bench(5, "MPKCV2         ", () => seq.Select(m => NumberTheory.SqrtModMPKCV2(m, p)).ToArray());
    GlobalStopWatch.Bench(50, "MPKCV2 discrLog", () => seq.Select(m => NumberTheory.SqrtModMPKCV2(m, p, ph)).ToArray());
    Console.WriteLine();
}

void testNthRootsMod(int n, int minP = 1)
{
    foreach (var p in Primes10000.Where(p => p < 47000 && p > minP && p % n == 1).Take(10))
    {
        var nbDigits = $"{p + 1}".Length;
        var fmt = $"{{0,{nbDigits}}}";
        var seqResidues = p.SeqLazy().Where(m => PowMod(m, (p - 1) / n, p) == 1).ToArray();
        var primRoots = NumberTheory.AllNthRootUnityMod(n, p).Order().ToArray();
        Console.WriteLine($"primitive {n}-throots of unity mod {p} and (w|{p}){n}");
        Console.WriteLine($"[{primRoots.ToDictionary(qi => qi, qi => PowMod(qi, (p - 1) / n, p)).GlueMap(", ")}]");
        if (primRoots.Length != n || primRoots.Any(w => w != 1 && !PowModEqualOne(w, n, p)))
            throw new($"errors roots of unity modulo");

        var (mul, rem) = FactorMultiplicity(n, p - 1);
        Console.WriteLine($"{p} - 1 = {n}^{mul} * {rem}");

        var seqAMM1 = seqResidues.ToDictionary(m => m, m => NumberTheory.NthRootANTV1(m, n, p))
            .ToDictionary(e => e.Key, e => primRoots.Select(f => e.Value * f % p).Min());
        var seqAMM2 = seqResidues.ToDictionary(m => m, m => NumberTheory.NthRootANTV1((BigInteger)m, n, p))
            .ToDictionary(e => e.Key, e => primRoots.Select(f => e.Value * f % p).Min());

        if (seqAMM1.Any(e => e.Value != seqAMM2[e.Key]))
        {
            var err = seqAMM1.Where(e => e.Value != seqAMM2[e.Key]).ToArray();
            foreach (var (e1, e2) in err)
            {
                Console.WriteLine($"{e1} => {e2} => {seqAMM2[e1]}");
                Console.WriteLine(PowModBigint(e2, n, p));
                Console.WriteLine(PowModBigint(seqAMM2[e1], n, p));
                Console.WriteLine();
            }

            Console.WriteLine($"nb err={err.Length}");
            throw new();
        }

        if (seqAMM1.Count < 50)
        {
            var fmtHead = $"{n}-thRoots({{0,{nbDigits}}}) mod {p}";
            seqAMM1.Select(e => (m: e.Key, root: primRoots.Select(o => o * e.Value % p).Min()))
                .Select(e =>
                    $"{string.Format(fmtHead, e.m)} = [{primRoots.Select(w => w * e.root % p).Glue(", ", fmt)}]")
                .Println($"Nb residues = {seqAMM1.Count}");
            Console.WriteLine();
        }

        GlobalStopWatch.Bench(50, "AMMint   ",
            () => seqResidues.Select(m => NumberTheory.NthRootANTV1(m, n, p)).ToArray());
        GlobalStopWatch.Bench(50, "AMMbigint",
            () => seqResidues.Select(m => NumberTheory.NthRootANTV1((BigInteger)m, n, p)).ToArray());
        GlobalStopWatch.Bench(50, "AMMint   ",
            () => seqResidues.Select(m => NumberTheory.NthRootANTV1(m, n, p)).ToArray());
        GlobalStopWatch.Bench(50, "AMMbigint",
            () => seqResidues.Select(m => NumberTheory.NthRootANTV1((BigInteger)m, n, p)).ToArray());

        Console.WriteLine("end");
        Console.WriteLine();
    }
}

void testCRT()
{
    var p = 151;
    var ph = NumberTheory.PreparePohligHellman(p);
    ph.crtTable.Println($"crt N={ph.N}");
    ph.primesDec.Println("Primes Decomp");
    ph.factors.Println("Factors Decomp");

    for (int a = 0; a < ph.N; ++a)
    {
        var ais = ph.crtTable.Select(e => a % e.mi).ToArray();
        var x = NumberTheory.CRT(ais, ph.crtTable, ph.N);
        Console.WriteLine($"a={a,3} => [{ais.Glue(", ", "{0,3}")}] => {x,3} {x == a}");
        if (x != a)
            throw new();
    }
    
    Console.WriteLine();
}

void testPH()
{
    var p = 97;
    var ph = NumberTheory.PreparePohligHellman(p);
    var g = NumberTheory.PrimitiveRootMod(p);
    for (int h = 1; h < ph.N; h++)
    {
        var a = NumberTheory.PohligHellman(h, g, ph);
        var ga = PowMod(g, a, p);
        Console.WriteLine($"g={g} h={h,3} a={a,3} g^a={ga,3}");
        if (h != ga)
            throw new();
    }
    
    Console.WriteLine();
}

{
    testTonelliShanks();

    testFFT();
    testNTT();
    testBenchNTT();
    
    testCRT();
    testPH();
    
    testSqrtMod();
    
    testNthRootsMod(3);
    testNthRootsMod(n: 3, minP: 10000);
    testNthRootsMod(5);
    testNthRootsMod(n: 5, minP: 10000);
    testNthRootsMod(7);
    testNthRootsMod(n: 7, minP: 10000);

    Console.WriteLine($"primRoot mean iter : {NumberTheory.nbLoopPrimRoot / NumberTheory.nbCallPrimRoot:F4}");
}
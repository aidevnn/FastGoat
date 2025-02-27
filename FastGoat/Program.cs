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
// RecomputeAllPrimesUpTo(1000000);

BigInteger AmodPSignedbigint(BigInteger m, BigInteger p)
{
    var m0 = AmodPbigint(m, p);
    return m0 > p / 2 ? m0 - p : m0;
}

// wikipedia Tonelli–Shanks algorithm
int SqrtMod(int m, int p)
{
    var S = PrimesDecomposition(p - 1).Count(pi => pi == 2);
    var Q = (p - 1) / (1 << S);
    var z = 0;
    for (var z0 = 1; z0 < p / 2; z0++)
    {
        var m0 = PowMod(z0, (p - 1) / 2, p);
        if (m0 == p - 1)
        {
            z = z0;
            break;
        }
    }

    var (M, c, t, R) = (S, PowMod(z, Q, p), PowMod(m, Q, p), PowMod(m, (Q + 1) / 2, p));
    var k = 0;
    while (true)
    {
        ++k;
        // Console.WriteLine(new { k, M, c, t, R });
        if (t == 0)
            return 0;
        if (t == 1)
            return R;

        var t0 = t;
        var i = (M + 1).SeqLazy().FirstOrDefault(i => PowMod(t0, 1 << i, p) == 1, -1);
        if (i == M && i > 1)
            return 0;

        if (i == -1)
            throw new();

        var b = PowMod(c, 1 << (M - i - 1), p);
        c = b * b % p;
        (M, t, R) = (i, t * c % p, R * b % p);
    }
}

BigInteger SqrtModBigint(BigInteger m, BigInteger p)
{
    m = AmodPSignedbigint(m, p);
    var S = PrimesDecompositionBigInt(p - 1).Count(pi => pi == 2);
    var Q = (p - 1) / (BigInteger.One << S);
    var z = BigInteger.Zero;
    for (BigInteger z0 = 1; z0 < p / 2; z0++)
    {
        var m0 = PowModBigint(z0, (p - 1) / 2, p);
        if (m0 == p - 1)
        {
            z = z0;
            break;
        }
    }

    var (M, c, t, R) = (S, PowModBigint(z, Q, p), PowModBigint(m, Q, p), PowModBigint(m, (Q + 1) / 2, p));
    var k = 0;
    while (true)
    {
        ++k;
        // Console.WriteLine(new { k, M, c, t, R });
        if (t == 0)
            return 0;
        if (t == 1)
            return AmodPSignedbigint(R, p);

        var t0 = t;
        var i = (M + 1).SeqLazy().FirstOrDefault(i => PowModBigint(t0, BigInteger.One << i, p) == 1, -1);
        if (i == M && i > 1)
            return 0;

        if (i == -1)
            throw new();

        var b = PowModBigint(c, BigInteger.One << (M - i - 1), p);
        c = b * b % p;
        (M, t, R) = (i, t * c % p, R * b % p);
    }
}

IEnumerable<BigInteger> Pow2NthRootMod(BigInteger m, BigInteger n, BigInteger p)
{
    if (!BigInteger.IsPow2(n))
        throw new();

    var r = SqrtModBigint(m, p);
    if (n == 2)
    {
        yield return -r;
        yield return r;
        yield break;
    }

    foreach (var r0 in Pow2NthRootMod(r, n / 2, p))
        yield return r0;

    foreach (var r0 in Pow2NthRootMod(-r, n / 2, p))
        yield return r0;
}

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
        var r = SqrtModBigint(n, p);
        var r2 = AmodPbigint(r * r, p);
        Console.WriteLine(new { n, p, r, r2 });
    }

    {
        var (n, p) = (1, 41);
        var r = SqrtModBigint(n, p);
        var r2 = AmodPbigint(r * r, p);
        Console.WriteLine(new { n, p, r, r2 });
    }

    {
        var (n, p) = (4, 241);
        Console.WriteLine($"{n}-thRoots(-1) mod {p}.");
        var seqTS = Pow2NthRootMod(-1, n, p).Select(i => AmodPbigint(i, p)).Order().ToArray();
        var seqPM = (p - 2).SeqLazy(2).Where(i => PowMod(i, n, p) == p - 1).Order().ToArray();
        Console.WriteLine($"[{seqTS.Glue(", ")}]");
        Console.WriteLine($"[{seqPM.Glue(", ")}]");
    }

    {
        var (n, p) = (8, 241);
        Console.WriteLine($"{n}-thRoots(-1) mod {p}.");
        var seqTS = Pow2NthRootMod(-1, n, p).Select(i => AmodPbigint(i, p)).Order().ToArray();
        var seqPM = (p - 2).SeqLazy(2).Where(i => PowMod(i, n, p) == p - 1).Order().ToArray();
        Console.WriteLine($"[{seqTS.Glue(", ")}]");
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
            var seqTS = Pow2NthRootMod(-1, n, p).ToArray();
            if (seqTS.Length != n || seqTS.Any(e => !PowModEqualOnebigint(e, N, p)))
                throw new("Fail");

            Console.WriteLine($"{n}-thRoots(-1) mod {p}. Pass");
            GlobalStopWatch.Bench(nb: 5, "TS first solution.", () => Pow2NthRootMod(-1, n, p).First());
            GlobalStopWatch.Bench(nb: 5, "TS all  solutions.", () => Pow2NthRootMod(-1, n, p).Last());
            Console.WriteLine();
        }
        // { N = 2048, n = 1024, p = 23026821121, log2p = 34 }
        // 1024-thRoots(-1) mod 23026821121. Pass
        // # TS first solution. Avg Time:6 ms Dev:1.200
        // # TS all  solutions. Avg Time:662 ms Dev:5.389
    }
}

void testFFT()
{
    var (N, p) = (16, 241);
    var n = N / 2;
    var w = new ZnBInt(p, Pow2NthRootMod(-1, n, p).First());
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
    var w = new ZnBInt(p, Pow2NthRootMod(-1, n, p).First());
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
    var a = RLWE.GenUnif(n, p);
    var b = RLWE.GenUnif(n, p);
    var ab1 = (a * b).ResModSigned(pm, p);
    var ab2 = RqMulNTT(a, b, ntt, intt);
    Console.WriteLine($"a     = {a}");
    Console.WriteLine($"b     = {b}");
    Console.WriteLine($"a * b = {ab1}");
    Console.WriteLine($"      = {ab2}");
    Console.WriteLine();
}

void testBenchNTT()
{
    var (N, p) = (1024, 12289);
    var n = N / 2;
    var w = new ZnBInt(p, Pow2NthRootMod(-1, n, p).First());
    var ni = (n * w.One).Inv();
    var ntt = NTT(n, w);
    var intt = ni * NTT(n, w.Inv()).T;
    var pm = FG.QPoly().Pow(n) + 1;
    var a = RLWE.GenUnif(n, p);
    var b = RLWE.GenUnif(n, p);

    GlobalStopWatch.Bench(5, "PolMul", () => (a * b).ResModSigned(pm, p));
    GlobalStopWatch.Bench(5, "NTT   ", () => RqMulNTT(a, b, ntt, intt));
    GlobalStopWatch.Bench(50, "PolMul", () => (a * b).ResModSigned(pm, p));
    GlobalStopWatch.Bench(50, "NTT   ", () => RqMulNTT(a, b, ntt, intt));
}

{
    // testFFT();
    // testNTT();
    // testBenchNTT();
}

int NthRoot(int a, int r, int p)
{
    var (S, q0) = FactorMultiplicity(r, p - 1);
    var q = (int)q0;
    if (q == p - 1)
        throw new($"r must divide q-1");

    if (PowMod(a, (p - 1) / r, p) != 1)
        throw new($"a must be {r}-pow residue mod p");

    var k = p.SeqLazy().First(u => (u * q + 1) % r == 0);
    foreach (var z in (p - r).SeqLazy(r).Where(z => PowMod(z, (p - 1) / r, p) != 1))
    {
        var (M, c, t, R) = (S, PowMod(z, q, p), PowMod(a, k * q, p), PowMod(a, (k * q + 1) / r, p));
        var z_ = PowMod(z, (p - 1) / r, p);
        var ct = 0;
        while (true)
        {
            ++ct;
            // var ti = PowMod(t, PowMod(r, M - 1, p), p);
            // var ci = PowMod(c, PowMod(r, M - 1, p), p);
            // var Ri = PowMod(R, r, p);
            // Console.WriteLine(new { k, u, z, z_, M, c, ci, t, ti, R, Ri, a });
            if (c == 1)
                break; // TODO: fix when r multiplicity in p-1 is >1 with residual primitive nthroots

            if (t == 0)
                return 0;
            if (t == 1)
                return R;

            var t0 = t;
            var i = (M + 1).SeqLazy().FirstOrDefault(i => PowMod(t0, PowMod(r, i, p), p) == 1, -1);
            if (i == M && i > 1)
                break; // TODO: fix when r multiplicity in p-1 is >1 with residual primitive nthroots

            if (i == -1)
                throw new();

            var b = PowMod(c, PowMod(r, M - i - 1, p), p);
            c = PowMod(b, r, p);
            (M, t, R) = (i, t * c % p, R * b % p);
        }
    }

    throw new();
}

void testSqrtMod()
{
    var (n, p) = (2, 12289);
    Console.WriteLine($"Sqrt mod {p}");
    var seq = p.SeqLazy().Where(m => LegendreJacobi(m, p) == 1).ToArray();

    var seqWP = seq.Select(m => SqrtMod(m, p)).SelectMany(e => new[] { e, p - e }).Distinct().Order().ToArray();
    var seqAMM = seq.Select(m => NthRoot(m, 2, p)).SelectMany(e => new[] { e, p - e }).Distinct().Order().ToArray();
    if (!seqWP.SequenceEqual(seqAMM))
    {
        Console.WriteLine(seqWP.Glue(", "));
        Console.WriteLine(seqAMM.Glue(", "));
        throw new();
    }

    GlobalStopWatch.Bench(50, "TS ", () => seq.Select(m => SqrtMod(m, p)).ToArray());
    GlobalStopWatch.Bench(50, "TSG", () => seq.Select(m => NthRoot(m, 2, p)).ToArray());
    GlobalStopWatch.Bench(50, "TS ", () => seq.Select(m => SqrtMod(m, p)).ToArray());
    GlobalStopWatch.Bench(50, "TSG", () => seq.Select(m => NthRoot(m, 2, p)).ToArray());
    Console.WriteLine();
}

void testNthRootsMod(int n, int minP = 1)
{
    foreach (var p in Primes10000.Where(p => p > minP && p % n == 1).Take(10))
    {
        var nbDigits = $"{p + 1}".Length;
        var fmt = $"{{0,{nbDigits}}}";
        var seq = p.SeqLazy().Where(m => PowMod(m, (p - 1) / n, p) == 1).ToArray();
        var q = p.SeqLazy().Where(m => PowModEqualOne(m, n, p)).Append(1).Order().ToArray();
        Console.WriteLine($"primitive {n}-throots of unity mod {p} and (w|{p}){n}");
        Console.WriteLine($"[{q.ToDictionary(qi => qi, qi => PowMod(qi, (p - 1) / n, p)).GlueMap(", ")}]");
        var (mul, rem) = FactorMultiplicity(n, p - 1);
        Console.WriteLine($"{p} - 1 = {n}^{mul} * {rem}");

        var seqRoots = seq.ToDictionary(m => m, m => NthRoot(m, n, p));

        if (seqRoots.Count < 50)
        {
            var fmtHead = $"{n}-thRoots({{0,{nbDigits}}}) mod {p}";
            seqRoots.Select(e => (m: e.Key, root: q.Select(o => o * e.Value % p).Min()))
                .Select(e => $"{string.Format(fmtHead, e.m)} = [{q.Select(w => w * e.root % p).Glue(", ", fmt)}]")
                .Println($"Nb residues = {seqRoots.Count}");
            if (seqRoots.Any(e => PowMod(e.Value, n, p) != e.Key))
            {
                Console.WriteLine(seqRoots.Where(e => PowMod(e.Value, n, p) != e.Key).GlueMap());
                throw new();
            }
        }

        Console.WriteLine("end");
        Console.WriteLine();
    }
}

{
    testSqrtMod();
    
    testNthRootsMod(3);
    testNthRootsMod(n: 3, minP: 10000);
    testNthRootsMod(5);
    testNthRootsMod(n: 5, minP: 10000);
    testNthRootsMod(7);
    testNthRootsMod(n: 7, minP: 10000);
}
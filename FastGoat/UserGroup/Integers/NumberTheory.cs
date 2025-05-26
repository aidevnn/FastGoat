using System.Numerics;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Padic;
using static FastGoat.Commons.IntExt;

namespace FastGoat.UserGroup.Integers;

public static class NumberTheory
{
    // Tonelli–Shanks generalisation algorithm
    // CHAPTER 2. BASIC ALGORITHMIC NUMBER THEORY. page 59
    // Mathematics of Public Key Cryptography. Version 2.0
    // Steven D.Galbraith
    public static int NthRootMPKCV2(int a, int r, int p, PohligHellmanInfos<int>? ph = null)
    {
        var q = (int)FactorMultiplicity(r, p - 1).rem;
        if (q == p - 1)
            throw new($"r must divide q-1");

        if (PowMod(a, (p - 1) / r, p) != 1)
            throw new($"a must be residue mod p");

        // Lemma 2.9.2. page58, [...] w is like a “ﬁrst approximation” to the 'rth' root of a modulo p
        var u = Bezout(-q, r).x % r;
        if (u < 0) u += r;
        var n = 100.SeqLazy().Select(_ => Rng.Next(2, p)).First(n => PowMod(n, (p - 1) / r, p) != 1);
        var (y, w, b) = (PowMod(n, q, p), PowMod(a, (u * q + 1) / r, p), PowMod(a, q, p));
        int j;
        if (ph is null)
            j = p.SeqLazy(1).First(k => PowMod(y, k * r, p) == b);
        else
            j = AmodP((1 - p) / r, p) * PohligHellman(b, y, ph) % p;

        return w * PowMod(y, -j, p) % p;
    }

    public static int SqrtModMPKCV2(int a, int p, PohligHellmanInfos<int>? ph = null) => NthRootMPKCV2(a, 2, p, ph);

    // Tonelli–Shanks generalisation algorithm
    // CHAPTER 2. BASIC ALGORITHMIC NUMBER THEORY. page 59
    // Mathematics of Public Key Cryptography. Version 2.0
    // Steven D.Galbraith
    public static BigInteger NthRootMPKCV2(BigInteger a, int r, BigInteger p, PohligHellmanInfos<BigInteger>? ph = null)
    {
        var q = FactorMultiplicity(r, p - 1).rem;
        if (q == p - 1)
            throw new($"r must divide q-1");

        if (PowModBigint(a, (p - 1) / r, p) != 1)
            throw new($"a must be residue mod p");

        // Lemma 2.9.2. page58, [...] w is like a “ﬁrst approximation” to the 'rth' root of a modulo p
        var u = BezoutBigInt(-q, r).Xa % r;
        if (u < 0) u += r;
        var n = 100.SeqLazy().Select(_ => DistributionExt.Dice(2, p - 1))
            .First(n => PowModBigint(n, (p - 1) / r, p) != 1);
        var (y, w, b) = (PowModBigint(n, q, p), PowModBigint(a, (u * q + 1) / r, p), PowModBigint(a, q, p));
        BigInteger j;
        if (ph is null)
            j = ((int)p).SeqLazy(1).First(k => PowModBigint(y, k * r, p) == b);
        else
            j = AmodPbigint((1 - p) / r, p) * PohligHellman(b, y, ph) % p;

        return w * PowModBigint(y, -j, p) % p;
    }

    public static BigInteger SqrtModMPKCV2(BigInteger a, BigInteger p, PohligHellmanInfos<BigInteger>? ph = null) =>
        NthRootMPKCV2(a, 2, p, ph);

    // L.Adleman, K.Manders, and G.Miller generalization of Tonelli's algorithm
    // Chapter 7. Solving Equations over Finite Fields. page 160
    // Algorithmic Number Theory Volume 1: Efficient Algorithms
    // Eric Bach and Jeffrey Shallit
    public static int NthRootANTV1(int a, int r, int p)
    {
        var (s, t0) = FactorMultiplicity(r, p - 1);
        var t = (int)t0;
        if (t == p - 1)
            throw new($"r must divide p-1");

        if (PowMod(a, (p - 1) / r, p) != 1)
            throw new($"a must be residue mod p");

        var rs = (p - 1) / t;
        var h = 1000.SeqLazy().Select(_ => Rng.Next(2, p - 1)).First(i => PowMod(i, (p - 1) / r, p) != 1);
        var (ar, at, g) = (PowMod(a, t, p), PowMod(a, rs, p), PowMod(h, t, p));
        var e = 0;
        var ri = 1;
        var rsi = (p - 1) / (t * r);
        for (int i = 0; i < s; i++)
        {
            for (var ei = 0; ei < r; ei++)
            {
                var ga = (PowMod(g, e + ei * ri, p) * ar) % p;
                if (PowMod(ga, rsi, p) == 1)
                {
                    e += ei * ri;
                    break;
                }
            }

            ri *= r;
            rsi /= r;
        }

        var r_ = InvModPbez(r, t);
        var (br, bt) = (PowMod(g, -e / r, p), PowMod(at, r_, p));
        var (alpha, beta) = Bezout(t, rs);
        var b = (PowMod(br, alpha, p) * PowMod(bt, beta, p)) % p;
        return b;
    }

    // L.Adleman, K.Manders, and G.Miller generalization of Tonelli's algorithm
    // Chapter 7. Solving Equations over Finite Fields. page 160
    // Algorithmic Number Theory Volume 1: Efficient Algorithms
    // Eric Bach and Jeffrey Shallit
    public static BigInteger NthRootANTV1(BigInteger a, int r, BigInteger p)
    {
        var (s, t) = FactorMultiplicity(r, p - 1);
        if (t == p - 1)
            throw new($"r must divide p-1");

        if (PowModBigint(a, (p - 1) / r, p) != 1)
            throw new($"a must be residue mod p");

        var rs = (p - 1) / t;
        var h = 1000.SeqLazy().Select(_ => DistributionExt.Dice(2, p - 1))
            .First(i => PowModBigint(i, (p - 1) / r, p) != 1);
        var (ar, at, g) = (PowModBigint(a, t, p), PowModBigint(a, rs, p), PowModBigint(h, t, p));
        var e = 0;
        var ri = 1;
        var rsi = (p - 1) / (t * r);
        for (int i = 0; i < s; i++)
        {
            for (var ei = 0; ei < r; ei++)
            {
                var ga = (PowModBigint(g, e + ei * ri, p) * ar) % p;
                if (PowModBigint(ga, rsi, p) == 1)
                {
                    e += ei * ri;
                    break;
                }
            }

            ri *= r;
            rsi /= r;
        }

        var r_ = InvModPbezbigint(r, t);
        var (br, bt) = (PowModBigint(g, -e / r, p), PowModBigint(at, r_, p));
        var (alpha, beta, _) = BezoutBigInt(t, rs);
        var b = (PowModBigint(br, alpha, p) * PowModBigint(bt, beta, p)) % p;
        return b;
    }

    public static int SqrtModANTV1(int a, int p)
    {
        var g = p.SeqLazy().Select(_ => Rng.Next(2, p)).First(i => LegendreJacobi(i, p) != 1);
        var (s, t0) = FactorMultiplicity(2, p - 1);
        var (t, e) = ((int)t0, 0);
        for (int i = 2; i <= s; i++)
        {
            var ag = (a * PowMod(g, -e, p)) % p;
            var agp = PowMod(ag, (p - 1) / (1 << i), p);
            if (agp != 1)
                e += 1 << (i - 1);
        }

        var h = (a * PowMod(g, -e, p)) % p;
        var b = (PowMod(g, e / 2, p) * PowMod(h, (t + 1) / 2, p)) % p;
        return b;
    }

    public static BigInteger SqrtModANTV1(BigInteger a, BigInteger p)
    {
        var g = 1000.SeqLazy().Select(_ => DistributionExt.Dice(2, p - 1)).First(i => LegendreJacobiBigint(i, p) != 1);
        var (s, t) = FactorMultiplicity(2, p - 1);
        var e = 0;
        for (int i = 2; i <= s; i++)
        {
            var ag = (a * PowModBigint(g, -e, p)) % p;
            var agp = PowModBigint(ag, (p - 1) / (1 << i), p);
            if (agp != 1)
                e += 1 << (i - 1);
        }

        var h = (a * PowModBigint(g, -e, p)) % p;
        var b = (PowModBigint(g, e / 2, p) * PowModBigint(h, (t + 1) / 2, p)) % p;
        return b;
    }

    public static EPoly<ZnInt> LegendreJacobiGf(EPoly<ZnInt> a)
    {
        var (p, n) = (a.P, a.F.Degree);
        var q = BigInteger.Pow(p, n);
        return a.FastPow((q - 1) / 2);
    }

    public static EPoly<ZnInt> SqrtFqANTV1(EPoly<ZnInt> a, EPoly<ZnInt> d)
    {
        if (a.IsZero() || a.IsOne())
            return a;
        
        var (p, n) = (a.P, a.F.Degree);
        var q = BigInteger.Pow(p, n);
        if (n == 1)
            return SqrtModANTV1(a[0].K, p) * a.One;
        
        var g = 1000.SeqLazy().Select(_ => DistributionExt.Dice(2, q - 1)).Select(i => d.FastPow(i))
            .First(g => !LegendreJacobiGf(g).IsOne());
        var (s, t) = FactorMultiplicity(2, q - 1);
        var e = 0;
        for (int i = 2; i <= s; i++)
        {
            var ag = a * g.Pow(-e);
            var agp = ag.FastPow((q - 1) / (1 << i));
            if (!agp.IsOne())
                e += 1 << (i - 1);
        }

        var h = a * g.Pow(-e);
        var b = g.Pow(e / 2) * h.FastPow((t + 1) / 2);
        return b;
    }

    public static EPoly<ZnInt> NthRootANTV1(EPoly<ZnInt> a, int r, EPoly<ZnInt> d)
    {
        if (r < 1 || !IsPrime(r))
            throw new($"r must be prime");

        if (r == 1)
            return a;

        var q = BigInteger.Pow(a.P, a.F.Degree);
        var (s, t0) = FactorMultiplicity(r, q - 1);
        var t = (int)t0;
        if (t == q - 1)
            throw new($"r must divide q - 1");

        if (!a.FastPow((q - 1) / r).IsOne())
            throw new($"a must be residue mod q");

        var rs = (q - 1) / t;
        var h = 1000.SeqLazy().Select(_ => DistributionExt.Dice(2, q - 1)).Select(i => d.FastPow(i))
            .First(di => !di.FastPow((q - 1) / r).IsOne());
        var (ar, at, g) = (a.FastPow(t), a.FastPow(rs), h.FastPow(t));
        var e = 0;
        var ri = 1;
        var rsi = (q - 1) / (t * r);
        for (int i = 0; i < s; i++)
        {
            for (var ei = 0; ei < r; ei++)
            {
                var ga = ar * g.Pow(e + ei * ri);
                if (ga.FastPow(rsi).IsOne())
                {
                    e += ei * ri;
                    break;
                }
            }

            ri *= r;
            rsi /= r;
        }

        var r_ = InvModPbez(r, t);
        var (br, bt) = (g.Pow(-e / r), at.Pow(r_));
        var (alpha, beta, _) = BezoutBigInt(t, rs);
        var b = br.FastPow(alpha) * bt.FastPow(beta);
        return b;
    }

    public static EPoly<ZnInt> NthRoot(EPoly<ZnInt> a, int r, EPoly<ZnInt> d)
    {
        return PrimesDecomposition(r).Aggregate(a, (acc, ri) => NthRootANTV1(acc, ri, d));
    }

    public static EPoly<ZnInt> PrimitiveRoot(EPoly<ZnInt> a)
    {
        var p = a.P;
        var n = a.F.Degree;
        var q = BigInteger.Pow(p, n);
        var dec = PrimesDecUnsafe(q - 1);
        var pis = dec.ToDictionary(e => e.Key, e => BigInteger.Pow(e.Key, e.Value));
        ++nbCallPrimRoot;
        var t = a.One;
        GlobalStopWatch.InfiniteLoopBreakerReset();
        while (pis.Count != 0)
        {
            GlobalStopWatch.InfiniteLoopBreaker(100, "Infinity loop PrimitiveRoot");
            ++nbLoopPrimRoot;
            var g = n.SeqLazy().Select(i => Rng.Next(p) * a.X.Pow(i)).ToVec().Sum();
            if (g.IsZero())
                continue;

            var S = pis.Keys.ToHashSet();
            foreach (var pi in S)
            {
                var gi = g.FastPow((q - 1) / pis[pi]);
                if (!gi.FastPow(pis[pi] / pi).IsOne())
                {
                    pis.Remove(pi);
                    t = t * gi;
                }
            }
        }

        return t;
    }

    // Tonelli-Shanks algorithm from wikipedia
    public static int SqrtModWP(int a, int p)
    {
        var (S, q0) = FactorMultiplicity(2, p - 1);
        var q = (int)q0;
        if (q == p - 1)
            throw new($"r must divide p-1");

        if (PowMod(a, (p - 1) / 2, p) != 1)
            throw new($"a must be quadratic residue mod p");

        var z = p.SeqLazy().Select(_ => Rng.Next(2, p)).First(i => LegendreJacobi(i, p) != 1);
        var (M, c, t, R) = (S, PowMod(z, q, p), PowMod(a, q, p), PowMod(a, (q + 1) / 2, p));
        while (true)
        {
            if (t == 0)
                return 0;

            if (t == 1)
                return R;

            var t0 = t;
            var i = (M + 1).SeqLazy().FirstOrDefault(i => PowMod(t0, 1 << i, p) == 1, -1);
            if (i == M && i > 1 || i == -1)
                break;

            var b = PowMod(c, 1 << (M - i - 1), p);
            c = b * b % p;
            (M, t, R) = (i, t * c % p, R * b % p);
        }

        throw new();
    }

    // Tonelli-Shanks algorithm from wikipedia
    public static BigInteger SqrtModWP(BigInteger a, BigInteger p)
    {
        var (S, q) = FactorMultiplicity(2, p - 1);
        if (q == p - 1)
            throw new($"r must divide p-1");

        if (PowModBigint(a, (p - 1) / 2, p) != 1)
            throw new($"a must be quadratic residue mod p");

        var one = BigInteger.One;
        var z = 1000.SeqLazy().Select(_ => DistributionExt.Dice(2, p - 1)).First(i => LegendreJacobiBigint(i, p) != 1);
        var (M, c, t, R) = (S, PowModBigint(z, q, p), PowModBigint(a, q, p), PowModBigint(a, (q + 1) / 2, p));
        while (true)
        {
            if (t == 0)
                return 0;

            if (t == 1)
                return R;

            var t0 = t;
            var i = (M + 1).SeqLazy().FirstOrDefault(i => PowModBigint(t0, one << i, p) == 1, -1);
            if (i == M && i > 1 || i == -1)
                break;

            if (i == -1)
                throw new();

            var b = PowModBigint(c, one << (M - i - 1), p);
            c = b * b % p;
            (M, t, R) = (i, t * c % p, R * b % p);
        }

        throw new();
    }

    public static double nbCallPrimRoot = 0;
    public static double nbLoopPrimRoot = 0;

    public static int PrimitiveRootMod(int p)
    {
        if (!IsPrime(p))
            throw new("p must be prime");

        var dec = PrimesDec(p - 1);
        var pis = dec.ToDictionary(e => e.Key, e => PowMod(e.Key, e.Value, p));
        ++nbCallPrimRoot;
        var t = 1;
        while (pis.Count != 0)
        {
            ++nbLoopPrimRoot;
            var g = Rng.Next(2, p);
            var S = pis.Keys.ToHashSet();
            foreach (var pi in S)
            {
                var gi = PowMod(g, (p - 1) / pis[pi], p);
                if (PowMod(gi, pis[pi] / pi, p) != 1)
                {
                    pis.Remove(pi);
                    t = t * gi % p;
                }
            }
        }

        return t;
    }

    public static BigInteger PrimitiveRootMod(BigInteger p)
    {
        if (!IsPrime(p))
            throw new("p must be prime");

        var dec = PrimesDec(p - 1);
        var pis = dec.ToDictionary(e => e.Key, e => PowModBigint(e.Key, e.Value, p));
        ++nbCallPrimRoot;
        BigInteger t = 1;
        while (pis.Count != 0)
        {
            ++nbLoopPrimRoot;
            var g = DistributionExt.Dice(2, p - 1);
            var S = pis.Keys.ToHashSet();
            foreach (var pi in S)
            {
                var gi = PowModBigint(g, (p - 1) / pis[pi], p);
                if (PowModBigint(gi, pis[pi] / pi, p) != 1)
                {
                    pis.Remove(pi);
                    t = t * gi % p;
                }
            }
        }

        return t;
    }

    public static int NthRootUnityMod(int n, int p) => PowMod(PrimitiveRootMod(p), (p - 1) / n, p);

    public static BigInteger NthRootUnityMod(BigInteger n, BigInteger p) =>
        PowModBigint(PrimitiveRootMod(p), (p - 1) / n, p);

    public static IEnumerable<int> AllNthRootUnityMod(int n, int p)
    {
        var w = NthRootUnityMod(n, p);
        var wi = w;
        for (int i = 0; i < n; i++)
        {
            yield return wi;
            wi = wi * w % p;
        }
    }

    public static IEnumerable<BigInteger> AllNthRootUnityMod(int n, BigInteger p)
    {
        var w = NthRootUnityMod(n, p);
        var wi = w;
        for (int i = 0; i < n; i++)
        {
            yield return wi;
            wi = wi * w % p;
        }
    }

    public static (int mi, int quo, int inv)[] CrtTable(int[] m)
    {
        var mod = m.Aggregate((mi, mj) => mi * mj);
        return m.Select(mi => (mi, quo: mod / mi))
            .Select(e => (e.mi, e.quo, InvModPbez(e.quo, e.mi)))
            .ToArray();
    }

    public static (BigInteger mi, BigInteger quo, BigInteger inv)[] CrtTable(BigInteger[] m)
    {
        var mod = m.Aggregate((mi, mj) => mi * mj);
        return m.Select(mi => (mi, quo: mod / mi))
            .Select(e => (e.mi, e.quo, InvModPbezbigint(e.quo, e.mi)))
            .ToArray();
    }

    public static int CRT(int[] a, (int mi, int quo, int inv)[] crtTable, int mod)
    {
        var x = 0;
        foreach (var (i, ai) in a.Index())
        {
            var (mi, quo, inv) = crtTable[i];
            x = (x + (ai % mi) * quo * inv) % mod;
        }

        return x;
    }

    public static BigInteger CRT(BigInteger[] a, (BigInteger mi, BigInteger quo, BigInteger inv)[] crtTable,
        BigInteger mod)
    {
        BigInteger x = 0;
        foreach (var (i, ai) in a.Index())
        {
            var (mi, quo, inv) = crtTable[i];
            x = (x + (ai % mi) * quo * inv) % mod;
        }

        return x;
    }

    public record PohligHellmanInfos<T>(
        T p,
        T N,
        T g,
        Dictionary<int, int> primesDec,
        Dictionary<(int pi, int j), T> factors,
        (T mi, T quo, T inv)[] crtTable);

    public static PohligHellmanInfos<int> PreparePohligHellman(int p)
    {
        if (!IsPrime(p))
            throw new("p must be prime");

        var N = p - 1;
        var g = PrimitiveRootMod(p);
        var dec = PrimesDec(N);
        var factors = dec.OrderBy(e => e.Key)
            .SelectMany(e => (e.Value + 1).SeqLazy().Select(j => ((e.Key, j), e.Key.Pow(j))))
            .ToDictionary(e => e.Item1, e => e.Item2);
        var crtTable = CrtTable(dec.OrderBy(e => e.Key).Select(e => factors[(e.Key, e.Value)]).ToArray());
        return new(p, N, g, dec, factors, crtTable);
    }

    public static PohligHellmanInfos<BigInteger> PreparePohligHellman(BigInteger p)
    {
        if (!IsPrime(p))
            throw new("p must be prime");

        var N = p - 1;
        var g = PrimitiveRootMod(p);
        var dec = PrimesDec(N);
        var factors = dec.OrderBy(e => e.Key)
            .SelectMany(e => (e.Value + 1).SeqLazy().Select(j => ((e.Key, j), BigInteger.Pow(e.Key, j))))
            .ToDictionary(e => e.Item1, e => e.Item2);
        var crtTable = CrtTable(dec.OrderBy(e => e.Key).Select(e => factors[(e.Key, e.Value)]).ToArray());
        return new(p, N, g, dec, factors, crtTable);
    }

    public static int PohligHellman(int h, int g, PohligHellmanInfos<int> ph)
    {
        var gPow = ph.factors.ToDictionary(e => e.Key, e => PowMod(g, ph.N / e.Value, ph.p));
        var hPow = ph.factors.ToDictionary(e => e.Key, e => PowMod(h, ph.N / e.Value, ph.p));
        var listAi = new int[ph.primesDec.Count];
        foreach (var (i, (pi, ei)) in ph.primesDec.Index())
        {
            var ai = 0;
            foreach (var j in ei.SeqLazy(1))
            {
                var (g0, h0) = (gPow[(pi, j)], hPow[(pi, j)]);
                var u = PowMod(g0, -ai, ph.p);
                h0 = h0 * u % ph.p;
                if (h0 != 1)
                {
                    g0 = gPow[(pi, 1)];
                    var b = 1;
                    var T = g0;
                    while (h0 != T)
                    {
                        ++b;
                        T = T * g0 % ph.p;
                    }

                    ai = (ai + b * ph.factors[(pi, j - 1)]) % ph.p;
                }
            }

            listAi[i] = ai;
        }

        return CRT(listAi, ph.crtTable, ph.N);
    }

    public static BigInteger PohligHellman(BigInteger h, BigInteger g, PohligHellmanInfos<BigInteger> ph)
    {
        var gPow = ph.factors.ToDictionary(e => e.Key, e => PowModBigint(g, ph.N / e.Value, ph.p));
        var hPow = ph.factors.ToDictionary(e => e.Key, e => PowModBigint(h, ph.N / e.Value, ph.p));
        var listAi = new BigInteger[ph.primesDec.Count];
        foreach (var (i, (pi, ei)) in ph.primesDec.Index())
        {
            BigInteger ai = 0;
            foreach (var j in ei.SeqLazy(1))
            {
                var (g0, h0) = (gPow[(pi, j)], hPow[(pi, j)]);
                var u = PowModBigint(g0, -ai, ph.p);
                h0 = h0 * u % ph.p;
                if (h0 != 1)
                {
                    g0 = gPow[(pi, 1)];
                    var b = 1;
                    var T = g0;
                    while (h0 != T)
                    {
                        ++b;
                        T = T * g0 % ph.p;
                    }

                    ai = (ai + b * ph.factors[(pi, j - 1)]) % ph.p;
                }
            }

            listAi[i] = ai;
        }

        return CRT(listAi, ph.crtTable, ph.N);
    }

    public static IEnumerable<int> Pow2NthRootsWP(int a, int r, int p)
    {
        var primRoots = AllNthRootUnityMod(r, p);
        var g = int.Log2(r).SeqLazy().Aggregate(a, (acc, _) => SqrtModWP(acc, p));
        return primRoots.Select(w => w * g % p);
    }

    public static IEnumerable<BigInteger> Pow2NthRootsWP(BigInteger a, int r, BigInteger p)
    {
        var primRoots = AllNthRootUnityMod(r, p);
        var g = int.Log2(r).SeqLazy().Aggregate(a, (acc, _) => SqrtModWP(acc, p));
        return primRoots.Select(w => w * g % p);
    }
}
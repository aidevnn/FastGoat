using System.Numerics;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Integers;
using Rq = FastGoat.Structures.VecSpace.KPoly<FastGoat.UserGroup.Integers.Rational>;

namespace FastGoat.UserGroup.LWE;

public partial class RLWE
{
    public static bool NoiseMode { get; set; } = true;

    public static void NoiseOn() => NoiseMode = true;

    public static void NoiseOff() => NoiseMode = false;

    public static Rq GenDiscrGauss(int n, double sigma)
    {
        return DistributionExt.DiscreteGaussianSample(n, sigma).ToKPoly(Rational.KZero());
    }

    public static Rq GenTernary(int n)
    {
        return DistributionExt.DiceSample(n, [-1, 0, 1]).ToKPoly(Rational.KZero());
    }

    public static Rq GenUnif(int n, int q)
    {
        return DistributionExt.DiceSample(n, -q / 2 + 1, q / 2).Select(e => new Rational(e)).ToKPoly();
    }

    public static Rq GenUnif(int n, BigInteger q)
    {
        return DistributionExt.DiceSampleBigInt(n, -q / 2 + 1, q / 2).Select(e => new Rational(e)).ToKPoly();
    }

    public static Rq GenUnif(int n, Rational q) => GenUnif(n, q.Num);

    public static Rq GenXpow(int n) => IntExt.RngSign * FG.QPoly().Pow(IntExt.Rng.Next(n));
    public static double Omega(int n) => double.Sqrt(n) / 2;

    public static double Alpha(int n) => 1.0 / (2 * double.Log2(n) * double.Log2(n) * Omega(n));

    public static double Sigma(int n, int p) => Alpha(n) * p;

    public static double Sigma(int n, Rational p) => Alpha(n) * p;

    public static Rational FirstPrimeEqualOneMod(Rational p, Rational bound)
    {
        var t = (int)double.Sqrt(p);
        var k = (bound / p).Trunc;
        var start = k * p + 1;
        return t.SeqLazy().Select(j => start + j * p)
            .FirstOrDefault(p0 => p0 > bound && IntExt.IsPrime(p0.Num), k.One);
    }

    public static Rational FirstPrimeEqualOneMod(int p, int bound) =>
        FirstPrimeEqualOneMod(new Rational(p), new Rational(bound));

    public static Rational FirstPrimeEqualOneMod(Rational p) => FirstPrimeEqualOneMod(p, "1");
    
    public static Rational FirstPrimeEqualOneMod(int p) => FirstPrimeEqualOneMod(new Rational(p), "1");

    public static Rational RlwePrime(int n)
    {
        var a = Alpha(n);
        var c = 2 * Omega(n);
        return FirstPrimeEqualOneMod(2 * n, (int)double.Ceiling(c / a));
    }

    public static Rq IntVecToRq(Vec<ZnInt64> v) => v.Select(e => new Rational(e.Signed)).ToKPoly();

    public static Rq SKBGV(int n)
    {
        var o = Rational.KOne();
        var arr = DistributionExt.DiceSample(n, [-o, o]).ToArray();
        var nb = (int)double.Sqrt(n);
        var zeros = DistributionExt.DiceSample(nb, (n - 1).Range()).ToArray();
        foreach (var i in zeros) arr[i] *= 0;
        return arr.ToKPoly();
    }

    public static RLWECipher PKBGV(Rq pm, Rq sk, Rational t, Rational q)
    {
        var noiseMode = NoiseMode ? 1 : 0;
        var n = pm.Degree;
        var sigma = Sigma(n, t);
        var e = GenDiscrGauss(n, sigma) * noiseMode;
        while (true)
        {
            var pkb = GenUnif(n, q);
            var pka = (t * e + pkb * sk).ResMod(pm);
            if (pka.NormInf() * 2 < q)
                continue;

            return new(pka.CoefsModSigned(q), pkb, pm, t, q);
        }
    }

    public static RLWECipher SWKBGV(Rq pm, Rq sk1, Rq sk2, Rational t, Rational q, Rational sp)
    {
        var noiseMode = NoiseMode ? 1 : 0;
        var n = pm.Degree;
        var sigma = Sigma(n, t);
        var eswk = GenDiscrGauss(n, sigma) * noiseMode;
        while (true)
        {
            var swkb = GenUnif(n, sp * q);
            var swka = (t * eswk + swkb * sk1 + sp * sk2).ResMod(pm);
            if (swka.NormInf() * 2 < sp * q)
                continue;

            return new(swka.CoefsModSigned(sp * q), swkb, pm, t, sp * q);
        }
    }

    public static Dictionary<Rational, (Rational nextMod, RLWECipher rlk)> LeveledSwitchKeyGenBGV(int level, Rq pm,
        Rq sk, Rational t, Rational[] seqMods, Rational sp)
    {
        var seqRlks = new Dictionary<Rational, (Rational nextMod, RLWECipher rlk)>();
        for (int i = 0; i <= level; i++)
        {
            var qi = seqMods[i];
            if (i == 0)
            {
                var cz = new RLWECipher(pm.Zero, pm.Zero, pm, t, qi);
                seqRlks[qi] = (t.One, cz);
                continue;
            }

            var spi = sp.Pow(i);
            var rlk = SWKBGV(pm, sk, sk.Pow(2), t, qi, spi);
            seqRlks[qi] = (seqMods[i - 1], rlk);
        }

        return seqRlks;
    }

    public static (Rq pm, Rq sk, Rational t, Rational[] primes, Rational sp, RLWECipher pk,
        Dictionary<Rational, (Rational nextMod, RLWECipher rlk)> rlks)
        SetupBGV(int N, int t0, int level, Rq sk, bool differentPrimes = true)
    {
        if (!int.IsPow2(N))
            throw new($"N = {N} must be 2^k");

        if (t0 <= 2 || !IntExt.Primes10000.Contains(t0))
            throw new($"T = {t0} must be an odd prime");

        if (sk.Degree + 1 != N / 2 || sk.Coefs.Any(e => e.Absolute > 1))
            throw new($"sk = {sk} must be a ternary vector {{-1,0,1}} of length N/2");

        var pm = FG.CyclotomicPolynomial(N);
        var t = new Rational(t0);

        Rational[] primes;

        var q = FirstPrimeEqualOneMod(N * t);
        if (differentPrimes)
        {
            var seqPrimes = new List<Rational>() { q };
            for (int i = 0; i < level; i++)
                seqPrimes.Add(FirstPrimeEqualOneMod(N * t, seqPrimes.Last()));

            primes = seqPrimes.Order().ToArray();
        }
        else
        {
            primes = Enumerable.Repeat(q, level + 1).ToArray();
        }

        if (primes.Length != level + 1 || primes[0].IsOne())
            throw new($"sequence moduli N={N} t={t} q={q} [{primes.Glue(", ")}]");

        var qL = primes.Aggregate((pi, pj) => pi * pj);
        var pk = PKBGV(pm, sk, t, qL);

        var seqMods = (level + 1).SeqLazy(1).Select(i => primes.Take(i).Aggregate((pi, pj) => pi * pj)).ToArray();
        var sp = FirstPrimeEqualOneMod(N * t, primes.Last());
        var seqRlks = LeveledSwitchKeyGenBGV(level, pm, sk, t, seqMods, sp);

        return (pm, sk, t, primes, sp, pk, seqRlks);
    }

    public static (Rq pm, Rq sk, Rational t, Rational[] primes, Rational sp, RLWECipher pk,
        Dictionary<Rational, (Rational nextMod, RLWECipher rlk)> rlks)
        SetupBGV(int N, int t, int level, bool differentPrimes)
    {
        return SetupBGV(N, t, level, SKBGV(N / 2), differentPrimes);
    }

    public static (Rq pm, Rq sk, Rational t, Rational[] primes, Rational sp, RLWECipher pk,
        Dictionary<Rational, (Rational nextMod, RLWECipher rlk)> rlks)
        SetupBGV(int N, int level, bool differentPrimes = true)
    {
        return SetupBGV(N, (int)RlwePrime(N / 2).Num, level, SKBGV(N / 2), differentPrimes);
    }

    public static RLWECipher EncryptBGV(Rq m, RLWECipher pk, bool noise = true)
    {
        var noiseMode = NoiseMode && noise ? 1 : 0;
        var (pm, t, q) = pk.PM_T_Q;
        var n = pm.Degree;
        var sigma = Sigma(n, t);
        var ea = GenDiscrGauss(n, sigma) * noiseMode;
        var eb = GenDiscrGauss(n, sigma) * noiseMode;
        var u = GenTernary(n) * noiseMode + pm.One * (1 - noiseMode);

        var m0 = m.ResModSigned(pm, t);
        var a = (u * pk.A + m0 + t * ea).ResModSigned(pm, q);
        var b = (u * pk.B + t * eb).ResModSigned(pm, q);
        return new RLWECipher(a, b, pm, t, q);
    }

    public static Rq DecryptBGV(RLWECipher cipher, Rq sk)
    {
        var (pm, t, q) = cipher.PM_T_Q;
        return (cipher.A - sk * cipher.B).ResModSigned(pm, q).CoefsModSigned(t);
    }

    public static Rq ErrorsBGV(RLWECipher cipher, Rq sk)
    {
        var (pm, _, q) = cipher.PM_T_Q;
        return (cipher.A - sk * cipher.B).ResModSigned(pm, q) - DecryptBGV(cipher, sk);
    }

    public static RLWECipher SwitchKeyBGV(RLWECipher cipher, RLWECipher swk)
    {
        var spi = swk.Q / cipher.Q;
        return spi * cipher.A - cipher.B * swk;
    }

    public static RLWECipher MulRelinBGV(RLWECipher cipher1, RLWECipher cipher2, RLWECipher rlk)
    {
        var (pm, t, q) = cipher1.PM_T_Q;
        var spi = rlk.Q / cipher1.Q;
        var c0 = (cipher1.A * cipher2.A).ResModSigned(pm, q);
        var c1 = (cipher1.B * cipher2.A + cipher1.A * cipher2.B).ResModSigned(pm, q);
        var c2 = (cipher1.B * cipher2.B).ResModSigned(pm, q);
        var c = new RLWECipher(c0, c1, pm, t, rlk.Q);
        return spi * c + c2 * rlk;
    }

    public static RLWECipher[] AKBGV(Rq sk, RLWECipher pk, Rational spi)
    {
        var (pm, t, q) = pk.PM_T_Q;
        var N = 2 * pm.Degree;
        return N.SeqLazy()
            .Select(j => SWKBGV(pm, sk, sk.Substitute(pm.X.Pow(j)).ResModSigned(pm, t), t, q, spi))
            .ToArray();
    }

    public static RLWECipher AutoMorphBGV(RLWECipher cipher, int k, RLWECipher ak)
    {
        var (pm, _, q) = cipher.PM_T_Q;
        var xk = pm.X.Pow(k);
        var spi = ak.Q / q;
        var a = cipher.A.Substitute(xk).ResModSigned(pm, q);
        var b = cipher.B.Substitute(xk).ResModSigned(pm, q);
        return spi * a - b * ak;
    }
}
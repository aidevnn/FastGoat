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

    public static Rq GenDiscrGauss(int n, double s = 3.2)
    {
        return DistributionExt.DiscreteGaussianSample(n, s).ToKPoly(Rational.KZero());
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
        var e = GenDiscrGauss(n) * noiseMode;
        while (true)
        {
            var c1 = GenUnif(n, q);
            var c0 = (t * e + c1 * sk).ResMod(pm);
            if (c0.NormInf() * 2 < q)
                continue;

            return new(c0.CoefsModSigned(q), c1, pm, t, q);
        }
    }

    public static RLWECipher SWKBGV(Rq pm, Rq sk1, Rq sk2, Rational t, Rational q, Rational sp)
    {
        var noiseMode = NoiseMode ? 1 : 0;
        var n = pm.Degree;
        var eswk = GenDiscrGauss(n) * noiseMode;
        while (true)
        {
            var swkb = GenUnif(n, sp * q);
            var swka = (t * eswk + swkb * sk1 + sp * sk2).ResMod(pm);
            if (swka.NormInf() * 2 < sp * q)
                continue;

            return new(swka.CoefsModSigned(sp * q), swkb, pm, t, sp * q);
        }
    }

    public static Dictionary<Rational, (Rational nextMod, RLWECipher[] skPow)>
        LeveledSwitchKeyGenBGV(int level, Rq pm, Rq sk, Rational t, Rational[] seqMods, Rational sp)
    {
        var n = pm.Degree;
        var seqSwks = new Dictionary<Rational, (Rational nextMod, RLWECipher[] skPow)>();
        for (int i = 0; i <= level; i++)
        {
            var qi = seqMods[i];
            if (i == 0)
            {
                var cz = new RLWECipher(pm.Zero, pm.Zero, pm, t, qi);
                seqSwks[qi] = (t.One, n.SeqLazy().Select(_ => cz).ToArray());
                continue;
            }

            var spi = sp.Pow(i);
            var swks = n.SeqLazy().Select(j => SWKBGV(pm, sk, sk.Pow(j), t, qi, spi)).ToArray();
            seqSwks[qi] = (seqMods[i - 1], swks);
        }

        return seqSwks;
    }

    public static Dictionary<Rational, (Rational nextMod, RLWECipher[] autSk)>
        LeveledAutMorphSwitchKeyGenBGV(int level, Rq pm, Rq sk, Rational t, Rational[] seqMods, Rational sp)
    {
        var n = pm.Degree;
        var N = 2 * n;
        var seqSwks = new Dictionary<Rational, (Rational nextMod, RLWECipher[] autSk)>();
        for (int i = 0; i <= level; i++)
        {
            var qi = seqMods[i];
            if (i == 0)
            {
                var cz = new RLWECipher(pm.Zero, pm.Zero, pm, t, qi);
                seqSwks[qi] = (t.One, N.SeqLazy().Select(_ => cz).ToArray());
                continue;
            }

            var spi = sp.Pow(i);
            var swks = N.SeqLazy()
                .Select(j => SWKBGV(pm, sk, sk.Substitute(pm.X.Pow(j)).ResModSigned(pm, t), t, qi, spi)).ToArray();
            seqSwks[qi] = (seqMods[i - 1], swks);
        }

        return seqSwks;
    }

    public static (Rq pm, Rq sk, Rational t, Rational[] primes, RLWECipher pk,
        Dictionary<Rational, (Rational nextMod, RLWECipher[] skPow)> swks)
        SetupBGV(int N, int t0, int level, Rq sk, bool differentPrimes = true)
    {
        if (IntExt.PrimesDecomposition(t0).Distinct().Count() != 1)
            throw new($"T = {t0} must be prime p^e");

        var n = N / 2;
        var pm = FG.QPoly().Pow(n) + 1;
        var t = new Rational(t0);

        Rational[] primes;
        if (differentPrimes)
            primes = IntExt.Primes10000.Where(t1 => t1 % N == 1 && t1 % t0 == 1).Take(level + 1)
                .Select(pi => new Rational(pi)).ToArray();
        else
            primes = Enumerable.Repeat(IntExt.Primes10000.First(t1 => t1 % N == 1 && t1 % t0 == 1) * t.One, level + 1)
                .ToArray();

        if (primes.Length != level + 1)
            throw new($"sequence moduli");

        var qL = primes.Aggregate((pi, pj) => pi * pj);
        var pk = PKBGV(pm, sk, t, qL);

        var seqMods = (level + 1).SeqLazy(1).Select(i => primes.Take(i).Aggregate((pi, pj) => pi * pj)).ToArray();
        var sp = new Rational(IntExt.Primes10000.FirstOrDefault(t1 => t1 > primes.Last() && t1 % t0 == 1, 1));
        var seqSwks = LeveledSwitchKeyGenBGV(level, pm, sk, t, seqMods, sp);

        return (pm, sk, t, primes, pk, seqSwks);
    }

    public static (Rq pm, Rq sk, Rational t, Rational[] primes, RLWECipher pk,
        Dictionary<Rational, (Rational nextMod, RLWECipher[] skPow)> swks)
        SetupBGV(int N, int t0, int level, bool differentPrimes)
    {
        return SetupBGV(N, t0, level, SKBGV(N / 2), differentPrimes);
    }

    public static (Rq pm, Rq sk, Rational t, Rational[] primes, RLWECipher pk,
        Dictionary<Rational, (Rational nextMod, RLWECipher[] skPow)> swks) SetupBGV(int N, int level,
            bool differentPrimes = true)
    {
        var t = IntExt.Primes10000.First(t1 => t1 % N == 1);
        return SetupBGV(N, t, level, SKBGV(N / 2), differentPrimes);
    }

    public static (Rq pm, Rq sk, Rational t, Rational[] primes, RLWECipher pk,
        Dictionary<Rational, (Rational nextMod, RLWECipher[] skPow)> swks) SetupBGV(int N, int t, int level)
    {
        return SetupBGV(N, t, level, SKBGV(N / 2));
    }

    public static (Rq pm, Rq sk, Rational t, Rational q, RLWECipher pk, RLWECipher rlk)
        KeyGenBGV(int n, int t0, Rq sk)
    {
        if (!IntExt.Primes10000.Contains(t0))
            throw new($"T = {t0} must be prime");

        var N = 2 * n;
        var (pm, _, t, primes, pk, swks) = SetupBGV(N, t0, level: 1, sk, differentPrimes: true);
        var q = primes[0];
        return (pm, sk, t, q, pk, swks[pk.Q].skPow[2]);
    }

    public static RLWECipher EncryptBGV(Rq m, RLWECipher pk, bool noise = true)
    {
        var noiseMode = NoiseMode && noise ? 1 : 0;
        var (pm, t, q) = pk.PM_T_Q;
        var n = pm.Degree;
        var ea = GenDiscrGauss(n) * noiseMode;
        var eb = GenDiscrGauss(n) * noiseMode;
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
        var (pm, t, q) = cipher.PM_T_Q;
        return (cipher.A - sk * cipher.B).ResModSigned(pm, q) - DecryptBGV(cipher, sk);
    }

    public static RLWECipher SwitchKeyBGV(RLWECipher cipher, RLWECipher swk)
    {
        var spi = swk.Q / cipher.Q;
        return spi * cipher.A - cipher.B * swk;
    }

    public static RLWECipher MulRelinBGV(RLWECipher cipher0, RLWECipher cipher1, RLWECipher rlk)
    {
        var (pm, t, q) = cipher0.PM_T_Q;
        var spi = rlk.Q / q;
        var d0 = (cipher0.A * cipher1.A).ResModSigned(pm, q);
        var d1 = (cipher0.A * cipher1.B + cipher0.B * cipher1.A).ResModSigned(pm, q);
        var d2 = (cipher0.B * cipher1.B).ResModSigned(pm, q);
        var c = new RLWECipher(d0, d1, pm, t, rlk.Q);
        return spi * c + d2 * rlk;
    }

    public static RLWECipher AutoMorphBGV(RLWECipher cipher, int k, RLWECipher ak)
    {
        var (pm, t, q) = cipher.PM_T_Q;
        var spi = ak.Q / q;
        var xk = pm.X.Pow(k);
        var a = cipher.A.Substitute(xk).ResMod(pm);
        var b = cipher.B.Substitute(xk).ResMod(pm);
        return spi * a - b * ak;
    }

    public static Rq ExtractVec(Vec<ZnInt64> v, int i = 0)
    {
        var x = FG.QPoly();
        var n = v.Length;
        return v.Select(e => new Rational(e.Signed))
            .Select((e, j) => i - j >= 0 ? e * x.Pow(i - j) : -e * x.Pow(n + i - j))
            .Aggregate((xi, xj) => xi + xj);
    }

    public static Rational[] ExtractArr(Rq poly, int n, int i = 0)
    {
        return n.SeqLazy(start: i, step: -1).Select(j => j >= 0 ? poly[j] : -poly[n + j]).ToArray();
    }

    public static (Rational ai, Rational[] bi)[] Extract(RLWECipher e)
    {
        var n = e.PM.Degree;
        return n.SeqLazy().Select(i => (ai: e.A[i], bi: ExtractArr(e.B, n, i).ToArray())).ToArray();
    }

    public static RLWECipher[] EXSK(Rq sk, RLWECipher pk)
    {
        var (pm, t, q) = pk.PM_T_Q;
        var n = pm.Degree;
        return n.SeqLazy().Select(i => EncryptBGV(sk[i] * pm.One, pk)).ToArray();
    }

    public static RLWECipher[] ExtractCoefs(RLWECipher cipher, RLWECipher[] exsk)
    {
        var extr = Extract(cipher);
        var x = cipher.PM.X;
        return extr.Select((e, i) => (e.ai - e.bi.Zip(exsk).Select(f => f.First * f.Second).ToVec().Sum()) * x.Pow(i))
            .ToArray();
    }
}
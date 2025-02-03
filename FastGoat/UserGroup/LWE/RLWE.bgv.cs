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
    
    public static (Rq pm, Rq sk, Rational t, Rational[] primes, RLWECipher pk, Dictionary<Rational, (Rational nextMod, RLWECipher rlk )>
        rlks)
        SetupBGV(int N, int t0, int level, Rq sk, bool differentPrimes = true)
    {
        if (!IntExt.Primes10000.Contains(t0))
            throw new($"T = {t0} must be prime");

        var n = N / 2;
        var pm = FG.QPoly().Pow(n) + 1;
        var t = new Rational(t0);

        Rational[] primes;
        if (differentPrimes)
            primes = IntExt.Primes10000.Where(t1 => t1 % N == 1 && t1 % t0 == 1).Take(level + 1)
                .Select(pi => new Rational(pi)).ToArray();
        else
            primes = Enumerable.Repeat(IntExt.Primes10000.First(t1 => t1 % N == 1 && t1 % t0 == 1) * t.One, level + 1).ToArray();

        if (primes.Length != level + 1)
            throw new($"sequence moduli");

        var qL = primes.Aggregate((pi, pj) => pi * pj);
        var epk = GenDiscrGauss(n);
        var c1pk = GenUnif(n, qL);
        var c0pk = (t * epk + c1pk * sk).ResModSigned(pm, qL);
        var pk = new RLWECipher(c0pk, c1pk, pm, t, qL);

        var seqMods = (level + 1).SeqLazy(1).Select(i => primes.Take(i).Aggregate((pi, pj) => pi * pj)).ToArray();

        var seqRlks = new Dictionary<Rational, (Rational nextMod, RLWECipher rlk )>();
        var sp = new Rational(IntExt.Primes10000.First(t1 => t1 > primes.Last() && t1 % N == 1 && t1 % t0 == 1));
        for (int i = 0; i <= level; i++)
        {
            var qi = seqMods[i];
            if (i == 0)
            {
                seqRlks[qi] = (t.One, new(pm.Zero, pm.Zero, pm, t, qi));
                continue;
            }

            var spi = sp.Pow(i);
            var erlk = GenDiscrGauss(n);
            var c1rlk = GenUnif(n, spi * qi);
            var c0rlk = (t * erlk + c1rlk * sk - spi * sk.Pow(2)).ResModSigned(pm, spi * qi);
            var rlk = new RLWECipher(c0rlk, c1rlk, pm, t, spi * qi);
            seqRlks[qi] = (seqMods[i - 1], rlk);
        }

        return (pm, sk, t, primes, pk, seqRlks);
    }

    public static (Rq pm, Rq sk, Rational t, Rational[] primes, RLWECipher pk, Dictionary<Rational, (Rational nextMod, RLWECipher rlk )>
        rlks) SetupBGV(int N, int level)
    {
        var t = IntExt.Primes10000.First(t1 => t1 % N == 1);
        return SetupBGV(N, t, level, SKBGV(N / 2));
    }

    public static (Rq pm, Rq sk, Rational t, Rational q, RLWECipher pk, RLWECipher rlk)
        KeyGenBGV(int n, int t0, Rq sk)
    {
        if (!IntExt.Primes10000.Contains(t0))
            throw new($"T = {t0} must be prime");

        var N = 2 * n;
        var (pm, _, t, primes, pk, rlks) = SetupBGV(N, t0, level: 1, sk, differentPrimes: true);
        var q = primes[0];
        return (pm, sk, t, q, pk, rlks[pk.Q].rlk);
    }

    public static RLWECipher EncryptBGV(Rq m, RLWECipher pk, bool noise = true)
    {
        var (pm, t, q) = pk.PM_T_Q;
        var n = pm.Degree;
        var ea = GenDiscrGauss(n);
        var eb = GenDiscrGauss(n).Zero;
        var u = GenTernary(n);
        if (!NoiseMode || !noise)
        {
            ea = eb = eb.Zero;
            u = u.One;
        }

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

    public static RLWECipher MulRelinBGV(RLWECipher cipher0, RLWECipher cipher1, RLWECipher rlk)
    {
        var (pm, t, q) = cipher0.PM_T_Q;
        var spi = rlk.Q / q;
        var d0 = (cipher0.A * cipher1.A).ResModSigned(pm, q);
        var d1 = (cipher0.A * cipher1.B + cipher0.B * cipher1.A).ResModSigned(pm, q);
        var d2 = (cipher0.B * cipher1.B).ResModSigned(pm, q);
        var c = new RLWECipher(d0, d1, pm, t, rlk.Q);
        return spi * c - d2 * rlk;
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

    public static (Rational ai, Rational[] bi)[] Extract(RLWECipher e, int n)
    {
        return n.SeqLazy().Select(i => (ai: e.A[i], bi: ExtractArr(e.B, n, i).ToArray())).ToArray();
    }

    public static (RLWECipher[] encXpow, RLWECipher[] exsk) EXSK(Rq pm, Rq sk, Rational t, Rational q, RLWECipher pk)
    {
        var n = pm.Degree;
        var encXpow = n.SeqLazy().Select(i => EncryptBGV(pm.X.Pow(i), pk)).ToArray();
        var exsk = n.SeqLazy().Select(i => EncryptBGV(sk[i].Signed(t) * pm.One, pk)).ToArray();
        return (encXpow, exsk);
    }

    public static RLWECipher[] ExtractCoefs(RLWECipher cipher, Rq pm, RLWECipher[] exsk)
    {
        var extr = Extract(cipher, pm.Degree);
        return extr.Select(e => e.ai - e.bi.ToVec().MulA(exsk.ToVec()).Sum()).ToArray();
    }
}
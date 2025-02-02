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

    public static (int t, Rational q0, Rational q1, Rational q2) SetupBGV(int N, int t, int level)
    {
        var tmp1 = IntExt.Primes10000.Where(t1 => t1 % N == 1 && t1 % t == 1).Take(1).ToArray();
        if (tmp1.Length == 0)
            throw new($"T = {t} and T^2 = {t * t}");

        var p2 = tmp1[0];
        var tmp2 = IntExt.Primes10000.Where(t1 => t1 % p2 == 1).Take(level + 1).ToArray();
        if (tmp2.Length <= level)
            throw new($"T = {t} and T^2 = {t * t} and P = {p2}");

        var (q0, p1) = tmp2.Select(e => new Rational(e)).Deconstruct();
        return (t, q0, q0 * p1, q0 * p1 * p2);
    }

    public static (int t, Rational q0, Rational q1, Rational q2) SetupBGV(int N, int level = 1)
    {
        var t = IntExt.Primes10000.First(t1 => t1 % N == 1);
        return SetupBGV(N, t, level);
    }

    public static int CheckParametersBGV(int N, int t, Rational q0, Rational q1, Rational q2)
    {
        // limit of precomputed primes
        var plast = IntExt.Primes10000.Last();
        if (q0 > new Rational(plast).Pow(2))
            return 1;

        // ciphertext modulus level are integers
        var (p0, p1, p2) = (q0, q1 / q0, q2 / q1);
        if (!p1.IsInteger() || !p2.IsInteger())
            return 2;

        // ciphertext modulus level0 is prime
        var dp0 = IntExt.PrimesDec(p0.Num);
        if (dp0.Count > 1 || dp0.Any(e => e.Value > 1))
            return 3;

        // ciphertext modulus level1 factor is prime
        var dp1 = IntExt.PrimesDec(p1.Num);
        if (dp1.Count > 1 || dp1.Any(e => e.Value > 1))
            return 4;

        // ciphertext modulus level2 factor is prime
        var dp2 = IntExt.PrimesDec(p2.Num);
        if (dp2.Count > 1 || dp2.Any(e => e.Value > 1))
            return 5;

        // ciphertext modulus factors equal one modulus plaintext modulus
        var one = Rational.KOne();
        // if (!p0.Mod(t).Equals(one) || !p0.Mod(N).Equals(one) || !p1.Mod(N).Equals(one) || !p2.Mod(N).Equals(one))
        if (!p0.Mod(p2).Equals(one) || !p1.Mod(p2).Equals(one) || !p2.Mod(t).Equals(one))
            return 6;

        // homomorphic multiplication
        if (p2.Num > p0.Num || p2.Num > p1.Num)
            return 7;

        return 0;
    }

    public static Rq SKBGV(int n)
    {
        var o = Rational.KOne();
        var arr = DistributionExt.DiceSample(n, [-o, o]).ToArray();
        var nb = (int)double.Sqrt(n);
        var zeros = DistributionExt.DiceSample(nb, (n - 1).Range()).ToArray();
        foreach (var i in zeros) arr[i] *= 0;
        return arr.ToKPoly();
    }

    public static (Rq pm, Rq sk, Rational t, Rational q, RLWECipher pk, RLWECipher rlk)
        KeyGenBGV(int n, int t0, Rq sk)
    {
        if (!IntExt.Primes10000.Contains(t0))
            throw new($"T = {t0} must be prime");

        // var check = CheckParametersBGV(2 * n, t0, q0, q1, q2);
        // if (check != 0)
        //     throw new($"Invalid parameters [{check}] {new { N = 2 * n, t0, q0, q1, q2 }}");

        var pm = FG.QPoly().Pow(n) + 1;
        var t = new Rational(t0);
        var q = new Rational(IntExt.Primes10000.First(t1 => t1 % (2 * n) == 1 && t1 % t0 == 1));
        var sp = new Rational(IntExt.Primes10000.First(t1 => t1 > q && t1 % (2 * n) == 1 && t1 % t0 == 1));

        var epk = GenDiscrGauss(n);
        var c1pk = GenUnif(n, sp * q);
        var c0pk = (t * epk + c1pk * sk).ResModSigned(pm, sp * q);
        var pk = new RLWECipher(c0pk, c1pk, pm, t, sp * q);

        var erlk = GenDiscrGauss(n);
        var c1rlk = GenUnif(n, sp * sp * q);
        var c0rlk = (t * erlk + c1rlk * sk - sp * sk.Pow(2)).ResModSigned(pm, sp * sp * q);
        var rlk = new RLWECipher(c0rlk, c1rlk, pm, t, sp * sp * q);

        return (pm, sk, t, q, pk, rlk);
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
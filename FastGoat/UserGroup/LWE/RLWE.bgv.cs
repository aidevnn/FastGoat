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
        return DistributionExt.DiceSample(n, 0, q - 1).Select(e => new Rational(e)).ToKPoly();
    }

    public static Rq GenUnif(int n, BigInteger q)
    {
        return DistributionExt.DiceSampleBigInt(n, 0, q - 1).Select(e => new Rational(e)).ToKPoly();
    }

    public static Rq GenUnif(int n, Rational q) => GenUnif(n, q.Num);
    public static Rq IntVecToRq(Vec<ZnInt64> v) => v.Select(e => new Rational(e.Signed)).ToKPoly();

    public static (Rq pm, Rq sk, Rational t, Rational q0, Rational q1, Rational q2, RLWECipher pk, RLWECipher rlk) 
        KeyGenBGV(int n, int t0, Rational q0, Rational q1, Rational q2, Rq sk)
    {
        var pm = FG.QPoly().Pow(n) + 1;
        var t = new Rational(t0);
        
        var epk = GenDiscrGauss(n);
        var c1pk = GenUnif(n, q1);
        var c0pk = (t * epk + c1pk * sk).ResModSigned(pm, q1);
        var pk = new RLWECipher(c0pk, c1pk, pm, t, q0, q1, q2);

        var p2 = q2 / q1;
        var erlk = GenDiscrGauss(n);
        var c1rlk = GenUnif(n, q1);
        var c0rlk = (t * erlk + c1rlk * sk + p2 * sk.Pow(2)).ResModSigned(pm, q2);
        var rlk = new RLWECipher(c0rlk, c1rlk, pm, t, q0, q1, q2);
        
        return (pm, sk, t, q0, q1, q2, pk, rlk);
    }

    public static (Rq pm, Rq sk, Rational t, Rational q0, Rational q1, Rational q2, RLWECipher pk, RLWECipher rlk) 
        KeyGenBGV(int n, int t, Rational q0, Rational q1, Rational q2)
    {
        var sk = 10000.SeqLazy().Select(_ => GenTernary(n))
            .First(s => !s[n - 1].IsZero() && s.Coefs.Count(e => e.IsZero()) <= n / 4);

        return KeyGenBGV(n, t, q0, q1, q2, sk);
    }

    public static RLWECipher EncryptBGV(Rq m, RLWECipher pk, bool noise = true)
    {
        var (pm, t, q0, q1, q2) = pk.PM_T_Q;
        var n = pm.Degree;
        var ea = GenDiscrGauss(n);
        var eb = GenDiscrGauss(n);
        var u = GenTernary(n);
        if (!NoiseMode || !noise)
        {
            ea = eb = eb.Zero;
            u = u.One;
        }

        var m0 = m.ResModSigned(pm, t);
        var a = (u * pk.A + m0 + t * ea).ResModSigned(pm, q1);
        var b = (u * pk.B + t * eb).ResModSigned(pm, q1);
        return new RLWECipher(a, b, pm, t, q0, q1, q2);
    }

    public static Rq DecryptBGV(RLWECipher cipher, Rq sk)
    {
        var (pm, t, q0, q1, _) = cipher.PM_T_Q;
        return (cipher.A - sk * cipher.B).ResModSigned(pm, q1).CoefsMod(t);
    }

    public static Rq ErrorsBGV(RLWECipher cipher, Rq sk)
    {
        var (pm, t, q0, q1, _) = cipher.PM_T_Q;
        return (cipher.A - sk * cipher.B).ResModSigned(pm, q1) - DecryptBGV(cipher, sk);
    }

    public static RLWECipher MulRelinBGV(RLWECipher cipher0, RLWECipher cipher1, RLWECipher rlk)
    {
        var (pm, t, q0, q1, q2) = cipher0.PM_T_Q;
        var p2 = q2 / q1;
        var d0 = (cipher0.A * cipher1.A).ResModSigned(pm, q1);
        var d1 = (cipher0.A * cipher1.B + cipher0.B * cipher1.A).ResModSigned(pm, q1);
        var d2 = (cipher0.B * cipher1.B).ResModSigned(pm, q1);

        var a = (p2 * d0 + d2 * rlk.A).ResModSigned(pm, q2);
        var b = (p2 * d1 + d2 * rlk.B).ResModSigned(pm, q2);
        return new RLWECipher(a, b, pm, t, q0, q1, q2).ModSwitch(q2, q1);
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
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

    public static (Rq pm, Rq sk, Rational t, Rational q, RLWECipher pk, RLWECipher rlk) KeyGenBGV(int n, int t, int q)
    {
        var pm = FG.QPoly().Pow(n) + 1;
        var sk = 10000.SeqLazy().Select(_ => GenTernary(n))
            .First(s => !s[n - 1].IsZero() && s.Coefs.Count(e => e.IsZero()) <= n / 4);
        
        var (T, Q) = (new Rational(t), new Rational(q));
        
        var epk = GenDiscrGauss(n);
        var c1pk = GenUnif(n, q);
        var c0pk = (t * epk + c1pk * sk).ResMod(pm, Q);
        var pk = new RLWECipher(c0pk, c1pk, pm, Q);
        
        var erlk = GenDiscrGauss(n);
        var c1rlk = GenUnif(n, q);
        var c0rlk = (t * erlk + c1rlk * sk + sk.Pow(2)).ResMod(pm, Q);
        var rlk = new RLWECipher(c0rlk, c1rlk, pm, Q);
        
        return (pm, sk, T, Q, pk, rlk);
    }

    public static (Rq pm, Rq sk, Rational t, Rational q, RLWECipher pk, RLWECipher rlk) KeyGenBGV(int n, int t, int q, Vec<ZnInt64> _sk)
    {
        var pm = FG.QPoly().Pow(n) + 1;
        var sk = _sk.Select(e => new Rational(e.Signed)).ToKPoly();
        
        var (T, Q) = (new Rational(t), new Rational(q));
        
        var epk = GenDiscrGauss(n);
        var c1pk = GenUnif(n, q);
        var c0pk = (t * epk + c1pk * sk).ResMod(pm, Q);
        var pk = new RLWECipher(c0pk, c1pk, pm, Q);
        
        var erlk = GenDiscrGauss(n);
        var c1rlk = GenUnif(n, q);
        var c0rlk = (t * erlk + c1rlk * sk + sk.Pow(2)).ResMod(pm, Q);
        var rlk = new RLWECipher(c0rlk, c1rlk, pm, Q);
        
        return (pm, sk, T, Q, pk, rlk);
    }

    public static RLWECipher EncryptBGV(Rq m, Rq pm, Rational t, Rational q, RLWECipher pk, bool noise = true)
    {
        var n = pm.Degree;
        var ea = GenDiscrGauss(n);
        var eb = GenDiscrGauss(n);
        var u = GenTernary(n);
        if (!NoiseMode || !noise)
        {
            ea = eb = eb.Zero;
            u = u.One;
        }

        var a = (u * pk.A + m + t * ea).ResMod(pm, q);
        var b = (u * pk.B + t * eb).ResMod(pm, q);
        return (a, b, pm, q);
    }

    public static Rq DecryptBGV(RLWECipher cipher, Rq pm, Rq sk, Rational t)
    {
        return (cipher.A - sk * cipher.B).ResMod(pm, t);
    }

    public static Rq ErrorsBGV(RLWECipher cipher, Rq pm, Rq sk, Rational t)
    {
        var diff = cipher.A - sk * cipher.B;
        var d = diff.ResMod(pm, t);
        return (diff - d).ResMod(pm, t);
    }

    public static RLWECipher MulRelinBGV(RLWECipher cipher0, RLWECipher cipher1, Rq pm, Rational q, RLWECipher rlk)
    {
        var d0 = (cipher0.A * cipher1.A).ResMod(pm, q);
        var d1 = (cipher0.A * cipher1.B + cipher0.B * cipher1.A).ResMod(pm, q);
        var d2 = (cipher0.B * cipher1.B).ResMod(pm, q);

        var a = (d0 + d2 * rlk.A).ResMod(pm, q);
        var b = (d1 + d2 * rlk.B).ResMod(pm, q);
        return (a, b, pm, q);
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
        var encXpow = n.SeqLazy().Select(i => EncryptBGV(pm.X.Pow(i), pm, t, q, pk)).ToArray();
        var exsk = n.SeqLazy().Select(i => EncryptBGV(sk[i].Signed(t) * pm.One, pm, t, q, pk)).ToArray();
        return (encXpow, exsk);
    }

    public static RLWECipher[] ExtractCoefs(RLWECipher cipher, Rq pm, RLWECipher[] exsk)
    {
        var extr = Extract(cipher, pm.Degree);
        return extr.Select(e => (e.ai - e.bi.ToVec().MulA(exsk.ToVec()).Sum()).CoefsMod(cipher.Q)).ToArray();
    }
}
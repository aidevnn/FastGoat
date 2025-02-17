using System.Numerics;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Integers;
using Rq = FastGoat.Structures.VecSpace.KPoly<FastGoat.UserGroup.Integers.Rational>;

namespace FastGoat.UserGroup.LWE;

public partial class RLWE
{
    public static Rq XpowA(int a, Rq pm, Rational q)
    {
        var x = pm.X;
        if (a == 0)
            return x.One;

        var n = pm.Degree;
        var sgn = a > 0 || a % n == 0 ? 0 : 1;
        return (x.Pow(IntExt.AmodP(a, n)) * (-1).Pow((a / n + sgn) % 2)).CoefsModSigned(q);
    }

    public static Vec<Rq> DecompRq(Rq a, Rational bs, Rational mod)
    {
        var size = (int)(BigInteger.Log10(mod.Num) / BigInteger.Log10(bs.Num)) + 1;
        var a0 = a.CoefsModSigned(mod);
        var queue = new Queue<Rq>();
        for (int i = 0; i < size; ++i)
        {
            var a1 = (a0 / bs).TruncPoly();
            var r = a0 - bs * a1;
            a0 = a1;
            queue.Enqueue(r);
        }

        return queue.ToVec();
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

    public static Rational GadgetBase(Rational t)
    {
        var u = (int)BigInteger.Log2(t.Num) + 1;
        return new Rational(BigInteger.Pow(2, u));
    }

    public static (Vec<RLWECipher> csm, Vec<RLWECipher> cm) EncryptRgswBGV(Rq m, RLWECipher pk, bool noiseMode = true)
    {
        var (pm, t, q) = pk.PM_T_Q;
        var B = GadgetBase(t);
        var size = (int)(BigInteger.Log10(q.Num) / BigInteger.Log10(B.Num)) + 1;
        var z = pm.Zero;
        RLWECipher ctZero() => EncryptBGV(z, pk, noiseMode);
        var cm = size.SeqLazy().Select(i => ctZero() + (m * B.Pow(i), z, pm, t, q)).ToVec();
        var csm = size.SeqLazy().Select(i => ctZero() + (z, -m * B.Pow(i), pm, t, q)).ToVec();
        return (csm, cm);
    }

    public static RLWECipher MulRgsw(RLWECipher c1, Vec<RLWECipher> cm, Vec<RLWECipher> csm)
    {
        var (_, t, q) = c1.PM_T_Q;
        var B = GadgetBase(t);
        var c1ac2m = DecompRq(c1.A, B, q).Zip(cm).Select(e => e.First * e.Second).ToVec().Sum();
        var c1bc2sm = DecompRq(c1.B, B, q).Zip(csm).Select(e => e.First * e.Second).ToVec().Sum();
        return c1ac2m - c1bc2sm;
    }

    public static (RLWECipher ct1, RLWECipher ctprep) CtPrep(RLWECipher ct)
    {
        var (pm, t, q) = ct.PM_T_Q;
        var N = new Rational(2 * pm.Degree);
        var a1 = (N * ct.A).CoefsModSigned(q);
        var b1 = (N * ct.B).CoefsModSigned(q);
        var a2 = ((N * ct.A - a1) / q).CoefsModSigned(N);
        var b2 = ((N * ct.B - b1) / q).CoefsModSigned(N);
        return ((a1, b1, pm, t, q), (a2, b2, pm, t, N));
    }

    public static ((Vec<RLWECipher> csm, Vec<RLWECipher> cm) plus, (Vec<RLWECipher> csm, Vec<RLWECipher> cm) minus)[]
        BRKgswBGV(Rq sk, RLWECipher pk)
    {
        var pm = pk.PM;
        var enc = (int k) => EncryptRgswBGV(pm.One * k, pk);
        var n = pm.Degree;
        return n.SeqLazy().Select(i => sk[i])
            .Select(c => (plus: c.IsOne() ? enc(1) : enc(0), minus: (-c).IsOne() ? enc(1) : enc(0)))
            .ToArray();
    }

    public static RLWECipher BlindRotategswBGV((Rational ai, Rational[] bi) ab, Rq f, 
        (Vec<RLWECipher> csm, Vec<RLWECipher> cm) rlwe0,
        ((Vec<RLWECipher> csm, Vec<RLWECipher> cm) plus, (Vec<RLWECipher> csm, Vec<RLWECipher> cm) minus)[] brk)
    {
        var (pm, t, qL) = brk[0].minus.cm[0].PM_T_Q;
        var B = GadgetBase(t);
        var n = pm.Degree;
        var x = pm.X;

        var beta = ab.bi;
        var alpha = (int)ab.ai.Num;
        var xalpha = XpowA(alpha, pm, qL);

        var (encSOne0, encOne0) = rlwe0;
        RLWECipher acc = ((f * xalpha).ResModSigned(pm, qL), x.Zero, pm, t, qL);
        
        for (int i = 0; i < n; i++)
        {
            var (encSi_plus, encSi_minus) = brk[i];
            var ai = (int)beta[i].Opp().Num;

            var exai = (XpowA(ai, pm, qL) - 1);
            var cxai = encSi_plus.cm.Select(e => exai * e).ToVec();
            var csxai = encSi_plus.csm.Select(e => exai * e).ToVec();

            var ex_ai = (XpowA(-ai, pm, qL) - 1);
            var cx_ai = encSi_minus.cm.Select(e => ex_ai * e).ToVec();
            var csx_ai = encSi_minus.csm.Select(e => ex_ai * e).ToVec();

            var acci = encOne0 + cxai + cx_ai;
            var sacci = encSOne0 + csxai + csx_ai;

            acc = MulRgsw(acc, acci, sacci);
        }

        return acc;
    }

    public static RLWECipher BlindRotategswBGV((Rational ai, Rational[] bi) ab,
        (Vec<RLWECipher> csm, Vec<RLWECipher> cm) rlwe0,
        ((Vec<RLWECipher> csm, Vec<RLWECipher> cm) plus, (Vec<RLWECipher> csm, Vec<RLWECipher> cm) minus)[] brk)
    {
        var (pm, t, qL) = brk[0].minus.cm[0].PM_T_Q;
        var x = pm.X;
        var c = pm.Degree / 2 - 1;
        var f = (2 * c + 1).SeqLazy(-c).Select(j => j * XpowA(j, pm, qL))
            .Aggregate(x.Zero, (acc, v) => acc + v).ResModSigned(pm, qL);

        return BlindRotategswBGV(ab, f, rlwe0, brk);
    }

    public static RLWECipher RepackingBGV(int n, RLWECipher[] accs, RLWECipher[] autSk)
    {
        var acc0 = accs[0];
        var (pm, t, q) = acc0.PM_T_Q;
        var d = pm.Degree;
        var x = pm.X;
        var CT = Ring.Matrix(acc0.One, d, d + 1);
        for (int i = 0; i < n; i++)
            CT[i, n] = accs[i];

        for (int k = n; k > 1; k /= 2)
        {
            for (int i = 0; i < k / 2; i++)
            {
                CT[i, k / 2] = CT[i, k] + x.Pow(k / 2) * CT[i + k / 2, k];
                var K = 1 + 2 * d / k;
                var crot = AutoMorphBGV(CT[i, k] - x.Pow(k / 2) * CT[i + k / 2, k], K, autSk[K]).ModSwitch(q);
                CT[i, k / 2] = CT[i, k / 2] + crot;
            }
        }

        return CT[0, 1];
    }

    public static (RLWECipher ctboot, RLWECipher ctsm) Bootstrapping(RLWECipher ct, RLWECipher pk, RLWECipher[] skAut,
        ((Vec<RLWECipher> csm, Vec<RLWECipher> cm) plus, (Vec<RLWECipher> csm, Vec<RLWECipher> cm) minus)[] brk)
    {
        var (pm, t, qL) = pk.PM_T_Q;
        var B = GadgetBase(t);
        var n = pm.Degree;

        // 1. Extract
        var (ct1, ctprep) = CtPrep(ct);
        var extract = Extract(ctprep);

        // 2. BlindRotate
        var rlwe0 = EncryptRgswBGV(pm.One, pk, noiseMode: false);
        var ni = (1 - qL) / n;
        var seqBR = new List<RLWECipher>();
        foreach (var ab in extract)
            seqBR.Add(ni * BlindRotategswBGV(ab, rlwe0, brk));

        // Step 3. Repacking
        var ctsm = RepackingBGV(n, seqBR.ToArray(), skAut);
        var ctboot = -ctsm + ct1.CoefsModSigned(qL);
        var Ni = (1 - t) / (2 * n);
        return (ctboot * Ni, ctsm);
    }
}
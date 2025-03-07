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

    public static Rq RotateStep(Rq f, Rational q, int n, int u, bool diff = true)
    {
        var k = diff ? 1 : 0;
        var p0 = n.SeqLazy().Select(_ => f.KZero).ToArray();
        for (int i = 0; i < n; i++)
        {
            var e = f[i];
            var idx = (int)IntExt.AmodPbigint(i + u, n);
            var fsgn = i + u > 0 || (i + u) % n == 0 ? 0 : 1;
            var sgn = (-1).Pow(((i + u) / n + fsgn) % 2);
            p0[idx] = (sgn * e - k * f[idx]).Signed(q);
        }
    
        return p0.ToKPoly();
    }

    public static Vec<Rq> DecompRq(Rq a, Rational[] primes)
    {
        return primes.Select(pi => a.CoefsModSigned(pi)).ToVec();
    }

    public static RLWECipher RotateStep(RLWECipher cipher, int u, bool diff = true)
    {
        var (pm, t, q) = cipher.PM_T_Q;
        var a = RotateStep(cipher.A, q, pm.Degree, u, diff);
        var b = RotateStep(cipher.B, q, pm.Degree, u, diff);
        return new(a, b, pm, t, q);
    }

    public static Vec<RLWECipher> DecompRNS(RLWECipher cipher, Rational[] primes)
    {
        var (pm, t, qL) = cipher.PM_T_Q;
        return primes.Select(pi => new RLWECipher(cipher.A.CoefsModSigned(pi), cipher.B.CoefsModSigned(pi), pm, t, qL))
            .ToVec();
    }

    public static Rq ExtractVec(Vec<ZnBigInt> v, int i = 0)
    {
        var x = FG.QPoly();
        var n = v.Length;
        return v.Select(e => new Rational(e.K))
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
        var pm = pk.PM;
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

    public static Rational[] RNSGadgetBase(Rational[] primes)
    {
        var qL = primes.Aggregate((pi, pj) => pi * pj);
        return primes
            .Select(e => (new Rational(IntExt.InvModPbezbigint((qL / e).Num, e.Num)).Signed(e) * (qL / e)).Signed(qL))
            .ToArray();
    }

    public static (Vec<RLWECipher> csm, Vec<RLWECipher> cm) EncryptRgswBGV(Rq m, RLWECipher pk, Rational[] B, 
        bool noiseMode = true)
    {
        var (pm, t, q) = pk.PM_T_Q;
        var z = pm.Zero;
        RLWECipher ctZero() => EncryptBGV(z, pk, noiseMode);
        var cm = B.Select(e => ctZero() + (m * e, z, pm, t, q)).ToVec();
        var csm = B.Select(e => ctZero() + (z, -m * e, pm, t, q)).ToVec();
        return (csm, cm);
    }

    public static RLWECipher MulRgsw(Vec<RLWECipher> c1, Vec<RLWECipher> cm, Vec<RLWECipher> csm)
    {
        var c1ac2m = c1.Zip(cm).Select(e => e.First.A * e.Second).ToVec().Sum();
        var c1bc2sm = c1.Zip(csm).Select(e => e.First.B * e.Second).ToVec().Sum();
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
        BRKgswBGV(Rq sk, RLWECipher pk, Rational[] B)
    {
        var pm = pk.PM;
        var enc = (int k) => EncryptRgswBGV(pm.One * k, pk, B);
        var n = pm.Degree;
        return n.SeqLazy().Select(i => sk[i])
            .Select(c => (plus: c.IsOne() ? enc(1) : enc(0), minus: (-c).IsOne() ? enc(1) : enc(0)))
            .ToArray();
    }

    public static RLWECipher BlindRotategswBGVslow((Rational ai, Rational[] bi) ab, Rq f,
        ((Vec<RLWECipher> csm, Vec<RLWECipher> cm) plus, (Vec<RLWECipher> csm, Vec<RLWECipher> cm) minus)[] brk, 
        Rational[] B, Rational[] primes)
    {
        var (pm, t, qL) = brk[0].minus.cm[0].PM_T_Q;
        var n = pm.Degree;
        var x = pm.X;

        var beta = ab.bi;
        var alpha = (int)ab.ai.Num;
        var xalpha = XpowA(alpha, pm, qL);

        var z = pm.Zero;
        var encOne0 = B.Select(e => new RLWECipher(pm.One * e, z, pm, t, qL)).ToVec();
        var encSOne0 = B.Select(e => new RLWECipher(z, -pm.One * e, pm, t, qL)).ToVec();
        
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

            acc = MulRgsw(DecompRNS(acc, primes), acci, sacci);
        }

        return acc;
    }

    public static RLWECipher BlindRotategswBGV((Rational ai, Rational[] bi) ab, Rq f,
        ((Vec<RLWECipher> csm, Vec<RLWECipher> cm) plus, (Vec<RLWECipher> csm, Vec<RLWECipher> cm) minus)[] brk,
        Rational[] B, Rational[] primes)
    {
        var (pm, t, qL) = brk[0].minus.cm[0].PM_T_Q;
        var n = pm.Degree;
        var x = pm.X;

        var beta = ab.bi;
        var alpha = (int)ab.ai.Num;

        var z = pm.Zero;
        var encOne0 = B.Select(e => new RLWECipher(pm.One * e, z, pm, t, qL)).ToVec();
        var encSOne0 = B.Select(e => new RLWECipher(z, -pm.One * e, pm, t, qL)).ToVec();

        var acc = RotateStep((f, x.Zero, pm, t, qL), alpha, diff: false);

        for (int i = 0; i < n; i++)
        {
            var (encSi_plus, encSi_minus) = brk[i];
            var ai = (int)beta[i].Opp().Num;

            var cxai = encSi_plus.cm.Select(e => RotateStep(e, ai)).ToVec();
            var csxai = encSi_plus.csm.Select(e => RotateStep(e, ai)).ToVec();

            var cx_ai = encSi_minus.cm.Select(e => RotateStep(e, -ai)).ToVec();
            var csx_ai = encSi_minus.csm.Select(e => RotateStep(e, -ai)).ToVec();

            var acci = encOne0 + cxai + cx_ai;
            var sacci = encSOne0 + csxai + csx_ai;

            acc = MulRgsw(DecompRNS(acc, primes), acci, sacci);
        }

        return acc;
    }

    public static RLWECipher BlindRotategswBGV((Rational ai, Rational[] bi) ab,
        ((Vec<RLWECipher> csm, Vec<RLWECipher> cm) plus, (Vec<RLWECipher> csm, Vec<RLWECipher> cm) minus)[] brk, 
        Rational[] B, Rational[] primes)
    {
        var (pm, _, qL) = brk[0].minus.cm[0].PM_T_Q;
        var x = pm.X;
        var c = pm.Degree / 2 - 1;
        var f = (2 * c + 1).SeqLazy(-c).Select(j => j * XpowA(j, pm, qL))
            .Aggregate(x.Zero, (acc, v) => acc + v).ResModSigned(pm, qL);

        return BlindRotategswBGV(ab, f, brk, B, primes);
    }

    public static RLWECipher RepackingBGV(RLWECipher[] accs, RLWECipher[] autSk)
    {
        var acc0 = accs[0];
        var (pm, _, q) = acc0.PM_T_Q;
        var n = pm.Degree;
        var x = pm.X;
        var CT = Ring.Matrix(acc0.One, n, n + 1);
        for (int i = 0; i < n; i++)
            CT[i, n] = accs[i];

        for (int k = n; k > 1; k /= 2)
        {
            for (int i = 0; i < k / 2; i++)
            {
                CT[i, k / 2] = CT[i, k] + x.Pow(k / 2) * CT[i + k / 2, k];
                var K = 1 + 2 * n / k;
                var crot = AutoMorphBGV(CT[i, k] - x.Pow(k / 2) * CT[i + k / 2, k], K, autSk[K]).ModSwitch(q);
                CT[i, k / 2] = CT[i, k / 2] + crot;
            }
        }

        return CT[0, 1];
    }

    public static RLWECipher Bootstrapping(RLWECipher cipher, RLWECipher pk,
        RLWECipher[] skAut,
        ((Vec<RLWECipher> csm, Vec<RLWECipher> cm) plus, (Vec<RLWECipher> csm, Vec<RLWECipher> cm) minus)[] brk, 
        Rational[] B, Rational[] primes)
    {
        var (pm, t, qL) = pk.PM_T_Q;
        var n = pm.Degree;

        // 1. Extract
        var (ct1, ctprep) = CtPrep(cipher);
        var extract = Extract(-ctprep);

        // 2. BlindRotate
        var ni = (1 - qL) / n;
        var seqBR = extract.Select(ab => ni * BlindRotategswBGV(ab, brk, B, primes)).ToArray();

        // Step 3. Repacking
        var Ni = (1 - t) / (2 * n);
        var ctsm = RepackingBGV(seqBR, skAut);
        return Ni * (ctsm + ct1.CoefsModSigned(qL));
    }
}
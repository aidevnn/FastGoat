using System.Collections.Concurrent;
using System.Numerics;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Integers;
using Rq = FastGoat.Structures.VecSpace.KPoly<FastGoat.UserGroup.Integers.Rational>;

namespace FastGoat.UserGroup.LWE;

public partial class RLWE
{
    public static KMatrix<ZnBigInt> NTT(int n, ZnBigInt w)
    {
        var ntt = Ring.Matrix(w.One, n, n);
        var wPow = (2 * n).SeqLazy().Select(i => w.Pow(i)).ToArray();
        for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            ntt[i, j] = wPow[(2 * i * j + j) % (2 * n)];

        return new(ntt);
    }

    public static KMatrix<ZnBigInt> CooleyTukey(Rq m, NTTInfos nttInfos)
    {
        var n = nttInfos.n;
        var z = nttInfos.w.Zero;
        var E = n.SeqLazy().Select(_ => z).ToArray();
        for (int j = 0; j < n / 2; j++)
        {
            var Aj = z;
            var Bj = z;
            for (int i = 0; i < n / 2; i++)
            {
                var w = nttInfos.wPows[(4 * i * j + 2 * i) % (2 * n)];
                Aj += w * m[2 * i].Num;
                Bj += w * m[2 * i + 1].Num;
            }

            var wj = nttInfos.wPows[2 * j + 1];
            E[j] = Aj + wj * Bj;
            E[j + n / 2] = Aj - wj * Bj;
        }

        return E.ToKMatrix(n);
    }

    public static Rq GentlemanSande(KMatrix<ZnBigInt> a, NTTInfos nttInfos)
    {
        var n = nttInfos.n;
        var z = nttInfos.w.Zero;
        var ni = (n * z.One).Inv();
        var E = n.SeqLazy().Select(_ => ni.Zero).ToArray();
        for (int i = 0; i < n / 2; i++)
        {
            var Ai = ni.Zero;
            var Bi = ni.Zero;
            for (int j = 0; j < n / 2; j++)
            {
                var (aj, aj_) = (a[j, 0], a[j + n / 2, 0]);
                Ai += (aj + aj_) * nttInfos.iwPows[2 * i * (2 * j + 1) % (2 * n)];
                Bi += (aj - aj_) * nttInfos.iwPows[(2 * i + 1) * (2 * j + 1) % (2 * n)];
            }

            E[2 * i] = ni * Ai;
            E[2 * i + 1] = ni * Bi;
        }

        return E.Select(e => new Rational(e.K)).ToKPoly();
    }

    public static NTTInfos PrepareNTT(int n, Rational t, ZnBigInt w)
    {
        var ntt = NTT(n, w);
        var ni = (n * w.One).Inv();
        var intt = ni * NTT(n, w.Inv()).T;
        var wPows = (2 * n).SeqLazy().Select(i => w.Pow(i)).ToArray();
        var iwPows = (2 * n).SeqLazy().Select(i => w.Inv().Pow(i)).ToArray();
        return new(n, w, ntt, wPows, t, intt, iwPows);
    }

    public static NTTInfos PrepareNTT(int n, Rational t, Rational[] primes)
    {
        var intPrimes = primes.Select(e => e.Num).ToArray();
        var q = intPrimes.Aggregate((pi, pj) => pi * pj);
        var crt = NumberTheory.CrtTable(intPrimes);
        var seq = intPrimes.Select(pi => NumberTheory.NthRootUnityMod(2 * n, pi)).ToArray();
        var w = new ZnBigInt(q, NumberTheory.CRT(seq, crt, q));
        return PrepareNTT(n, t, w);
    }

    public static KMatrix<ZnBigInt> XpowA(ZnBigInt a, NTTInfos nttInfos)
    {
        if (a.IsZero())
            return nttInfos.ntt.GetCol(0);

        var n = nttInfos.n;
        var a0 = a.K;
        var sgn = a0 > 0 || a0 % n == 0 ? 0 : 1;
        var e = (int)IntExt.AmodPbigint(a0, n);
        var s = (-1).Pow((int)((a0 / n + sgn) % 2));
        return s * nttInfos.ntt.GetCol(e);
    }

    public static KMatrix<ZnBigInt> Rq2NTT(Rq a, NTTInfos nttInfos)
    {
        var n = nttInfos.n;
        var A = a.CoefsExtended(n - 1).Select(e => e.Num * nttInfos.ntt.KOne).ToKMatrix(n);
        return nttInfos.ntt * A;
    }

    public static KMatrix<ZnBigInt> MulNTT(KMatrix<ZnBigInt> a, KMatrix<ZnBigInt> b)
    {
        return a.Zip(b).Select(e => e.First * e.Second).ToKMatrix(a.M);
    }

    public static Rq NTT2Rq(KMatrix<ZnBigInt> a, NTTInfos nttInfos)
    {
        return (nttInfos.intt * a).Select(e => new Rational(e.K)).ToKPoly();
    }

    public static NTTCipher RLWECipher2NTT(RLWECipher cipher, NTTInfos nttInfos)
    {
        return new(Rq2NTT(cipher.A, nttInfos), Rq2NTT(cipher.B, nttInfos), nttInfos);
    }

    public static RLWECipher NTTCipher2RLWE(NTTCipher cipher, NTTInfos nttInfos)
    {
        var n = nttInfos.n;
        var pm = FG.QPoly().Pow(n) + 1;
        var qL = new Rational(cipher.A.KOne.Mod);
        return new(NTT2Rq(cipher.A, nttInfos), NTT2Rq(cipher.B, nttInfos), pm, nttInfos.t, qL);
    }

    public static Vec<NTTCipher> RLWECipher2NTT(Vec<RLWECipher> ciphers, NTTInfos nttInfos)
    {
        return ciphers.Select(cipher => RLWECipher2NTT(cipher, nttInfos)).ToVec();
    }

    public static KMatrix<ZnBigInt> KMatMod(KMatrix<ZnBigInt> a, BigInteger mod1, BigInteger mod2) =>
        a.Select(e => new ZnBigInt(mod1, new ZnBigInt(mod2, e.K).K)).ToKMatrix(a.M);

    public static KMatrix<ZnBigInt> KMatMod(KMatrix<ZnBigInt> a, Rational pi) => KMatMod(a, a.KOne.Mod, pi.Num);

    public static NTTCipher NTTCipherMod(NTTCipher cipher, Rational pi)
    {
        return new NTTCipher(KMatMod(cipher.A, pi), KMatMod(cipher.B, pi), cipher.NttInfos);
    }

    public static Vec<NTTCipher> DecompRNS(NTTCipher cipher, Rational[] primes)
    {
        var a = GentlemanSande(cipher.A, cipher.NttInfos);
        var b = GentlemanSande(cipher.B, cipher.NttInfos);
        return primes.Select(pi => (A: a.CoefsModSigned(pi), B: b.CoefsModSigned(pi)))
            .Select(c => new NTTCipher(CooleyTukey(c.A, cipher.NttInfos), CooleyTukey(c.B, cipher.NttInfos),
                cipher.NttInfos)).ToVec();

        // var a = cipher.NttInfos.intt * cipher.A;
        // var b = cipher.NttInfos.intt * cipher.B;
        // return primes.Select(pi => (A: KMatMod(a, pi), B: KMatMod(b, pi)))
        //     .Select(c => new NTTCipher(cipher.NttInfos.ntt * c.A, cipher.NttInfos.ntt * c.B, cipher.NttInfos))
        //     .ToVec();
    }

    public static Vec<KMatrix<ZnBigInt>> DecompRNS(KMatrix<ZnBigInt> a, Rational[] primes, NTTInfos nttInfos)
    {
        var a0 = NTT2Rq(a, nttInfos);
        return primes.Select(pi => Rq2NTT(a0.CoefsModSigned(pi), nttInfos)).ToVec();
    }

    public static Vec<KMatrix<ZnBigInt>> DecompRq(KMatrix<ZnBigInt> a, Rational[] primes)
    {
        return primes.Select(pi => KMatMod(a, pi)).ToVec();
    }

    public static KMatrix<ZnBigInt> RotateStepNTT(KMatrix<ZnBigInt> f, NTTInfos nttInfos, int u, bool diff = true)
    {
        var k = diff ? 1 : 0;
        var n = nttInfos.n;
        var u0 = (u * f.KOne).K;
        var sgn = u0 > 0 || u0 % n == 0 ? 0 : 1;
        var idx = (int)IntExt.AmodPbigint(u0, n);
        var s = (-1).Pow((int)((u0 / n + sgn) % 2));
        return f.Zip(nttInfos.ntt.GetCol(idx)).Select(e => s * e.First * e.Second - k * e.First).ToKMatrix(n);
    }

    public static NTTCipher RotateStepNTT(NTTCipher cipher, int u, bool diff = true)
    {
        var a = RotateStepNTT(cipher.A, cipher.NttInfos, u, diff);
        var b = RotateStepNTT(cipher.B, cipher.NttInfos, u, diff);
        return (a, b, cipher.NttInfos);
    }

    public static NTTCipher MulRgsw(Vec<NTTCipher> c1, Vec<NTTCipher> cm, Vec<NTTCipher> csm)
    {
        var c1ac2m = c1.Zip(cm).Select(e => e.First.A * e.Second).ToVec().Sum();
        var c1bc2sm = c1.Zip(csm).Select(e => e.First.B * e.Second).ToVec().Sum();
        return c1ac2m - c1bc2sm;
    }

    public static (Vec<NTTCipher> csm, Vec<NTTCipher> cm) EncryptRgswNTTBGV(Rq m, RLWECipher pk, Rational[] B,
        NTTInfos nttInfos, bool noiseMode = true)
    {
        var (pm, t, q) = pk.PM_T_Q;
        var z = pm.Zero;
        RLWECipher ctZero() => EncryptBGV(z, pk, noiseMode);
        var cm = B.Select(e => RLWECipher2NTT(ctZero() + (m * e, z, pm, t, q), nttInfos)).ToVec();
        var csm = B.Select(e => RLWECipher2NTT(ctZero() + (z, -m * e, pm, t, q), nttInfos)).ToVec();
        return (csm, cm);
    }

    public static ((Vec<NTTCipher> csm, Vec<NTTCipher> cm) plus, (Vec<NTTCipher> csm, Vec<NTTCipher> cm) minus)[]
        BRKgswNTTBGV(Rq sk, RLWECipher pk, Rational[] B, NTTInfos nttInfos)
    {
        var pm = pk.PM;
        var enc = (int k) => EncryptRgswNTTBGV(pm.One * k, pk, B, nttInfos);
        var n = pm.Degree;
        return n.SeqLazy().Select(i => sk[i])
            .Select(c => (plus: c.IsOne() ? enc(1) : enc(0), minus: (-c).IsOne() ? enc(1) : enc(0)))
            .ToArray();
    }

    public static RLWECipher BlindRotateNTTgswBGV((Rational ai, Rational[] bi) ab, Rq f,
        ((Vec<NTTCipher> csm, Vec<NTTCipher> cm) plus, (Vec<NTTCipher> csm, Vec<NTTCipher> cm) minus)[] brk,
        Rational[] B, Rational[] primes)
    {
        var cipher = brk[0].minus.cm[0].Zero;
        var nttInfos = cipher.NttInfos;
        var n = nttInfos.n;
        var one = nttInfos.w.One;

        var beta = ab.bi;
        var alpha = ab.ai.Num * one;
        var xalpha = XpowA(alpha, nttInfos);

        var encOne0 = B.Select(e => e.Num * one * cipher.One).ToVec();
        var encSOne0 = B.Select(e => e.Num * one * cipher.ESK).ToVec();

        var acc = cipher.One * xalpha * Rq2NTT(f, nttInfos);
        for (int i = 0; i < n; i++)
        {
            var (encSi_plus, encSi_minus) = brk[i];
            var ai = beta[i].Opp().Num * one;

            var exai = XpowA(ai, nttInfos) - cipher.KOne;
            var cxai = encSi_plus.cm.Select(e => exai * e).ToVec();
            var csxai = encSi_plus.csm.Select(e => exai * e).ToVec();

            var ex_ai = XpowA(-ai, nttInfos) - cipher.KOne;
            var cx_ai = encSi_minus.cm.Select(e => ex_ai * e).ToVec();
            var csx_ai = encSi_minus.csm.Select(e => ex_ai * e).ToVec();

            var acci = encOne0 + cxai + cx_ai;
            var sacci = encSOne0 + csxai + csx_ai;

            acc = MulRgsw(DecompRNS(acc, primes), acci, sacci);
        }

        return NTTCipher2RLWE(acc, nttInfos);
    }

    public static RLWECipher BlindRotateNTTgswBGV((Rational ai, Rational[] bi) ab,
        ((Vec<NTTCipher> csm, Vec<NTTCipher> cm) plus, (Vec<NTTCipher> csm, Vec<NTTCipher> cm) minus)[] brk,
        Rational[] B, Rational[] primes)
    {
        var cipher = brk[0].minus.cm[0].Zero;
        var nttInfos = cipher.NttInfos;
        var n = nttInfos.n;
        var pm = FG.QPoly().Pow(n) + 1;
        var qL = new Rational(nttInfos.w.Mod);
        var x = pm.X;
        var c = n / 2 - 1;
        var f = (2 * c + 1).SeqLazy(-c).Select(j => j * XpowA(j, pm, qL))
            .Aggregate(x.Zero, (acc, v) => acc + v).ResModSigned(pm, qL);

        return BlindRotateNTTgswBGV(ab, f, brk, B, primes);
    }

    public static Vec<NTTCipher> Clone(Vec<NTTCipher> v) => v.Select(c => c.Clone()).ToVec();

    public static RLWECipher BootstrappingNTT(RLWECipher cipher, RLWECipher pk,
        RLWECipher[] skAut,
        ((Vec<NTTCipher> csm, Vec<NTTCipher> cm) plus, (Vec<NTTCipher> csm, Vec<NTTCipher> cm) minus)[] brkntt,
        Rational[] B, Rational[] primes)
    {
        var (pm, t, qL) = pk.PM_T_Q;
        var n = pm.Degree;

        // 1. Extract
        var (ct1, ctprep) = CtPrep(cipher);
        var extract = Extract(-ctprep);

        // 2. BlindRotate
        var ni = (1 - qL) / n;
        var seqBR = extract.Select(ab => ni * BlindRotateNTTgswBGV(ab, brkntt, B, primes)).ToArray();

        // var bag = new ConcurrentDictionary<int, RLWECipher>();
        // var opt = new ParallelOptions() { MaxDegreeOfParallelism = 4 };
        // Parallel.ForEach(extract.Index().ToArray(), opt, (e, _) =>
        // {
        //     var acc = ni * BlindRotateNTTgswBGV(e.Item, brkntt, B, primes);
        //     bag[e.Index] = acc;
        // });
        // var seqBR = bag.OrderBy(e => e.Key).Select(e => e.Value).ToArray();

        // Step 3. Repacking
        var Ni = (1 - t) / (2 * n);
        var ctsm = RepackingBGV(seqBR, skAut);
        return Ni * (ctsm + ct1.CoefsModSigned(qL));
    }
}
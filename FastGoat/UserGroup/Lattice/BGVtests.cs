using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Integers;
using Rq = FastGoat.Structures.VecSpace.KPoly<FastGoat.UserGroup.Integers.Rational>;

namespace FastGoat.UserGroup.Lattice;

public static class BGVtests
{
    static BGVtests()
    {
        Ring.DisplayPolynomial = MonomDisplay.StarCaret;
    }

    public static void FirstBGV()
    {
        var (n, q0) = (16, 424242);
        var t0 = n / 2 - 1;
        var (pm, sk, t, q, pk, rlk) = FHE.KeyGenBGV(n, t0, q0);
        FHE.Show(pm, sk, t, q, pk, rlk);
        
        for (int i = 0; i < 100; i++)
        {
            var mi = FHE.GenUnif(n, t0);
            var ct = FHE.EncryptBGV(mi, pm, t, q, pk);
            var mf = FHE.DecryptBGV(ct, pm, sk, t);
            Console.WriteLine(mi);
            Console.WriteLine(mf);
            Console.WriteLine();
            if (!(mf - mi).IsZero())
                throw new($"step[{i}]");
        }
    }

    public static void HEAddBGV()
    {
        var (n, q0) = (16, 424242);
        var t0 = n / 2 - 1;
        var (pm, sk, t, q, pk, rlk) = FHE.KeyGenBGV(n, t0, q0);
        FHE.Show(pm, sk, t, q, pk, rlk);

        for (int i = 0; i < 100; ++i)
        {
            var m1 = FHE.GenUnif(n, t0);
            var m2 = FHE.GenUnif(n, t0);
            var m1m2 = (m1 + m2).CoefsMod(t);
            var ct1 = FHE.EncryptBGV(m1, pm, t, q, pk);
            var ct2 = FHE.EncryptBGV(m2, pm, t, q, pk);
            var ct = FHE.AddBGV(ct1, ct2, q);

            var decrypt = FHE.DecryptBGV(ct, pm, sk, t);
            Console.WriteLine($"m1 + m2:{m1m2}");
            Console.WriteLine($"decrypt:{decrypt}");
            Console.WriteLine();
            if (!(decrypt - m1m2).IsZero())
                throw new($"step[{i}]");
        }
    }

    public static void HEMulBGV()
    {
        var (n, t0, q0) = (16, 8, 2.Pow(6));
        var (pm, sk, t, q, pk, rlk) = FHE.KeyGenBGV(n, t0, q0);
        FHE.Show(pm, sk, t, q, pk, rlk);

        for (int i = 0; i < 100; ++i)
        {
            var m1 = FHE.GenUnif(n, t0);
            var m2 = FHE.GenUnif(n, t0);
            var m1m2 = (m1 * m2).ResMod(pm, t);
            var ct1 = FHE.EncryptBGV(m1, pm, t, q, pk);
            var ct2 = FHE.EncryptBGV(m2, pm, t, q, pk);
            var ct = FHE.MulRelinBGV(ct1, ct2, pm, q, rlk);

            var decrypt = FHE.DecryptBGV(ct, pm, sk, t);
            Console.WriteLine($"m1 * m2:{m1m2}");
            Console.WriteLine($"decrypt:{decrypt}");
            Console.WriteLine();
            if (!(decrypt - m1m2).IsZero())
                throw new($"step[{i}]");
        }
    }

    public static void HEMul2BGV()
    {
        var l = 4;
        var n = 1 << l;
        var t0 = n;
        var q0 = n * n;
        var (pm, sk, t, q, pk, rlk) = FHE.KeyGenBGV(n, t0, q0);
        FHE.Show(pm, sk, t, q, pk, rlk);

        var (es2, es) = FHE.ESKBGV(pm, sk, t, q, pk);
        es.Show("Encrypt(SK)");
        es2.Show("Encrypt(SK^2)");
        Console.WriteLine();

        for (int k = 0; k < 50; ++k)
        {
            var m1 = FHE.GenUnif(n, t);
            var m2 = FHE.GenUnif(n, t);
            var cm1 = FHE.EncryptBGV(m1, pm, t, q, pk);
            var cm2 = FHE.EncryptBGV(m2, pm, t, q, pk);
            
            var cm1m2 = FHE.MulSwkBGV(cm1, cm2, pm, q, es, es2);
            
            var dm1m2 = FHE.DecryptBGV(cm1m2, pm, sk, t);
            var m1m2 = (m1 * m2).ResMod(pm, t);
            Console.WriteLine($"m1     :{m1}");
            Console.WriteLine($"m2     :{m2}");
            Console.WriteLine($"m1 * m2:{m1m2}");
            Console.WriteLine($"       :{dm1m2}");
            Console.WriteLine();

            if (!dm1m2.Equals(m1m2))
                throw new($"step[{k}]");
        }
    }

    public static void TrackingErrors()
    {
        var depth = 8;
        var size = 1 << depth;
        var n = 7;
        var (t0, q0) = (59, 7 * 59);
        var (pm, sk, t, q, pk, rlk) = FHE.KeyGenBGV(n, t0, q0);
        FHE.Show(pm, sk, t, q, pk, rlk);
        var (es2, es) = FHE.ESKBGV(pm, sk, t, q, pk);

        for (int k = 0; k < 50; ++k)
        {
            Console.WriteLine($"Size {size} Depth:{depth}");
            var seqMsg = size.SeqLazy().Select(_ => FHE.GenUnif(n, t0)).ToArray();
            var seqCipher = seqMsg.Select(m => FHE.EncryptBGV(m, pm, t, q, pk)).ToArray();
            
            var sum = seqMsg.Aggregate((e0, e1) => (e0 + e1).ResMod(pm, t));
            var mul = seqMsg.Aggregate((e0, e1) => (e0 * e1).ResMod(pm, t));

            var qSum = new Queue<RLWECipher>(seqCipher);
            while (qSum.Count > 1) qSum.Enqueue(FHE.AddBGV(qSum.Dequeue(), qSum.Dequeue(), q));
            var dSum = FHE.DecryptBGV(qSum.Dequeue(), pm, sk, t);
            
            var qMul1 = new Queue<RLWECipher>(seqCipher);
            while (qMul1.Count > 1) qMul1.Enqueue(FHE.MulRelinBGV(qMul1.Dequeue(), qMul1.Dequeue(), pm, q, rlk));
            var dMul1 = FHE.DecryptBGV(qMul1.Dequeue(), pm, sk, t);
            
            var qMul2 = new Queue<RLWECipher>(seqCipher);
            while (qMul2.Count > 1) qMul2.Enqueue(FHE.MulSwkBGV(qMul2.Dequeue(), qMul2.Dequeue(), pm, q, es, es2));
            var dMul2 = FHE.DecryptBGV(qMul2.Dequeue(), pm, sk, t);
            
            Console.WriteLine($"Sum  Equal: {dSum.Equals(sum)}");
            Console.WriteLine($"sum :{sum}");
            Console.WriteLine($"    :{dSum}");
            Console.WriteLine($"Mul1 Equal: {dMul1.Equals(mul)}");
            Console.WriteLine($"Mul2 Equal: {dMul2.Equals(mul)}");
            Console.WriteLine($"prod:{mul}");
            Console.WriteLine($"    :{dMul1}");
            Console.WriteLine($"    :{dMul2}");
            Console.WriteLine();
            if (!sum.Equals(dSum) || !mul.Equals(dMul1) || !mul.Equals(dMul2))
                throw new($"step[{k}]");
        }
    }

    public static void HEEvalAuto()
    {
        var (n, t0, q0) = (16, 8, 2.Pow(6));
        var (pm, sk, t, q, pk, rlk) = FHE.KeyGenBGV(n, t0, q0);
        FHE.Show(pm, sk, t, q, pk, rlk);
        var ak = FHE.AKBGV(pm, sk, t, q, pk);
        var x = pm.X;

        foreach (var k in IntExt.Coprimes(2 * n))
        {
            var xk = x.Pow(k);
            for (var i = 0; i < 50; ++i)
            {
                var m = FHE.GenUnif(n, t0);
                var mk = m.Substitute(xk).ResMod(pm, t);
                Console.WriteLine($"m   :{m}");
                Console.WriteLine($"m^{k,2}:{mk}");

                var cm = FHE.EncryptBGV(m, pm, t, q, pk);
                var ck = FHE.AutoMorphBGV(cm, k, pm, t, q, pk, ak);
                var dk = FHE.DecryptBGV(ck, pm, sk, t);
                Console.WriteLine($"    :{dk}");
                if (!dk.Equals(mk))
                    throw new($"k:{k} step[{i}] {dk.Div(mk)}");

                Console.WriteLine();
            }
        }
    }

    public static void SwitchKeys()
    {
        var (n, t0, q0) = (16, 8, 2.Pow(6));
        var (pm, sk1, t, q, pk1, rlk1) = FHE.KeyGenBGV(n, t0, q0);
        FHE.Show(pm, sk1, t, q, pk1, rlk1);

        var (_, sk2, _, _, pk2, rlk2) = FHE.KeyGenBGV(n, t0, q0);
        FHE.Show(pm, sk2, t, q, pk2, rlk2);
        var s2s1 = FHE.SWKBGV(pm, sk1, t, q, pk2);
        var s1s2 = FHE.SWKBGV(pm, sk2, t, q, pk1);

        for (int i = 0; i < 50; ++i)
        {
            {
                var m = FHE.GenUnif(n, t0);
                Console.WriteLine($"m :{m}");
                var cm = FHE.EncryptBGV(m, pm, t, q, pk1);

                var dm = FHE.DecryptBGV(FHE.SwitchKeysBGV(cm, pm, t, q, pk2, s2s1), pm, sk2, t);
                Console.WriteLine($"  :{dm}");

                if (!dm.Equals(m))
                    throw new($"step[{i}]");
            }

            {
                var m = FHE.GenUnif(n, t0);
                Console.WriteLine($"m :{m}");
                var cm = FHE.EncryptBGV(m, pm, t, q, pk2);

                var dm = FHE.DecryptBGV(FHE.SwitchKeysBGV(cm, pm, t, q, pk1, s1s2), pm, sk1, t);
                Console.WriteLine($"  :{dm}");

                if (!dm.Equals(m))
                    throw new($"step[{i}]");
            }

            Console.WriteLine();
        }
    }

    public static void BlindRotate()
    {
        var l = 4;
        var n = 1 << l;
        var t0 = 8;
        var q0 = n * 8;
        var (pm, sk, t, q, pk, rlk) = FHE.KeyGenBGV(n, t0, q0);
        FHE.Show(pm, sk, t, q, pk, rlk);

        // var brk = FHE.BRKBGV(pm, sk, t, q, pk);
        var brk = FHE.BRKgswBGV(pm, sk, t, q, pk);

        var x = pm.X;
        var c = n / 4;
        var f = (2 * c + 1).SeqLazy(-c).Select(j => j * FHE.XpowA(j, pm, t))
            .Aggregate(x.Zero, (acc, v) => acc + v).ResMod(pm, t);
        // f = FHE.GenUnif(n, t);
        var nbDigits = $"{q}".Length;
        var fmt = $"{{0,{nbDigits}}}";
        var set = new HashSet<int>();
        var s = n.SeqLazy().Select(i => sk[i].Mod(q)).ToArray();

        for (int k = 0; k < 50; ++k)
        {
            var ai = new Rational(DistributionExt.Dice(1, t0));
            var bi = FHE.GenUnif(n, t0).CoefsExtended(n - 1);
            // var acc = FHE.BlindRotateBGV((ai, bi), f, pm, q, pk, rlk, brk);
            var acc = FHE.BlindRotateRgswBGV((ai, bi), f, pm, q, pk, brk);
            var actual = FHE.DecryptBGV(acc, pm, sk, t);

            Console.WriteLine($"f           :{f}");
            Console.WriteLine($"f           :[{f.CoefsExtended(n - 1).Glue(", ", fmt)}]");
            Console.WriteLine($"ai          :{ai}");
            Console.WriteLine($"bi          :[{bi.Glue(", ", fmt)}]");
            Console.WriteLine($"blind rotate:{actual}");
            Console.WriteLine($"            :[{actual.CoefsExtended(n - 1).Glue(", ", fmt)}]");

            // Testing result
            var u = (ai - FHE.InnerProd(bi, s, q)).Mod(q);
            var u0 = (int)u.Num;
            var expected = (FHE.XpowA(u0, pm, q) * f).ResMod(pm, t);
            Console.WriteLine($"u= a - <b,s>:{u}");
            Console.WriteLine($"f*X^{u,-4}    :{expected}");
            Console.WriteLine($"f*X^{u,-4}    :[{expected.CoefsExtended(n - 1).Glue(", ", fmt)}]");
            var factor = ((int)q.Num).SeqLazy(1).First(k0 => (k0 * actual).CoefsMod(q).Equals(expected));
            Console.WriteLine($"factor      :{factor}");
            Console.WriteLine();
            set.Add(factor);
        }

        Console.WriteLine($"Factors:{set.Order().Glue(", ", fmt)}");
    }

    public static void Bootstrapping()
    {
        var l = 3;
        var n = 1 << l;
        var t0 = n / 2;
        var q0 = 4 * n * n;
        var (pm, sk, t, q, pk, rlk) = FHE.KeyGenBGV(n, t0, q0);
        FHE.Show(pm, sk, t, q, pk, rlk);

        var Q = 32 * q;
        var fact = 2 * n;
        // var brk = FHE.BRKBGV(pm, sk, t, Q, pk);
        var brk = FHE.BRKgswBGV(pm, sk, t, Q, pk);
        var ak = FHE.AKBGV(pm, sk, t, Q, pk);

        var nbTrial = 50;
        for (int k = 0; k < nbTrial; ++k)
        {
            var m1 = FHE.GenUnif(n, t);
            var cm1 = FHE.EncryptBGV(m1, pm, t, q, pk);
            cm1.Show($"ct m1:{m1}");

            // var ctboot = FHE.Bootstrapping(cm1, pm, q, Q, pk, rlk, brk, ak, fact);
            var ctboot = FHE.BootstrappingRgsw(cm1, pm, q, Q, pk, brk, ak, fact);
            ctboot.Show($"ctboot Q = {Q}");
            var decBoot = FHE.DecryptBGV(ctboot, pm, sk, t);

            Console.WriteLine($"decrypt ctboot:{decBoot}");
            Console.WriteLine($"m1            :{m1}");
            Console.WriteLine();
            if (!decBoot.Equals(m1))
                throw new("decrypt");
        }
    }

    public static void RgswMul()
    {
        var (n, t0, q0) = (16, 8, 2.Pow(6));
        var (pm, sk, t, q, pk, rlk) = FHE.KeyGenBGV(n, t0, q0);
        FHE.Show(pm, sk, t, q, pk, rlk);

        for (int l = 0; l < 50; ++l)
        {
            var m1 = FHE.GenUnif(n, t);
            var cm1 = FHE.EncryptBGV(m1, pm, t, q, pk);

            var m2 = FHE.GenUnif(n, t);
            var (csm2, cm2) = FHE.EncryptRgswBGV(m2, pm, t, q, pk);
            var cm1m2 = FHE.SubBGV(FHE.KMulBGV(cm2, cm1.A, pm, q), FHE.KMulBGV(csm2, cm1.B, pm, q), q);

            var m1m2 = (m1 * m2).ResMod(pm, t);
            var dm1m2 = FHE.DecryptBGV(cm1m2, pm, sk, t);
            Console.WriteLine($"m1  :{m1}");
            Console.WriteLine($"m2  :{m2}");
            Console.WriteLine($"m1m2:{m1m2}");
            Console.WriteLine($"    :{dm1m2}");
            Console.WriteLine();
            if (!dm1m2.Equals(m1m2))
                throw new($"diff:{(m1m2 - dm1m2).CoefsMod(q)}");
        }
    }
    
    public static void TestAll()
    {
        FirstBGV();
        HEAddBGV();
        HEMulBGV();
        HEMul2BGV();
        TrackingErrors();
        HEEvalAuto();
        SwitchKeys();
        RgswMul();
        BlindRotate();
        Bootstrapping();
    }
}
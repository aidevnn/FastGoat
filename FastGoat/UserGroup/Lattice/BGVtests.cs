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
        var (n, t0, q0) = (16, 8, 2.Pow(6));
        var (pm, sk, t, q, pk, rlk) = FHE.KeyGenBGV(n, t0, q0);
        FHE.Show(pm, sk, t, q, pk, rlk);

        // var (esk2, esk) = FHE.EncryptRgswBGV(sk, pm, t, q, pk);
        var esk = FHE.EncryptBGV(sk, pm, t, q, pk);
        var esk2 = FHE.EncryptBGV((sk * sk).ResMod(pm, t), pm, t, q, pk);

        for (int i = 0; i < 100; ++i)
        {
            var m1 = FHE.GenUnif(n, t0);
            var m2 = FHE.GenUnif(n, t0);
            var m1m2 = (m1 * m2).ResMod(pm, t);
            var cm1 = FHE.EncryptBGV(m1, pm, t, q, pk);
            var cm2 = FHE.EncryptBGV(m2, pm, t, q, pk);
            
            var scm2 = FHE.SubBGV(FHE.KMulBGV(esk, cm2.A, pm, q), FHE.KMulBGV(esk2, cm2.B, pm, q), q);
            var cm1m2 = FHE.SubBGV(FHE.KMulBGV(cm2, cm1.A, pm, q), FHE.KMulBGV(scm2, cm1.B, pm, q), q);
            // var scm2 = cm2.A * esk - cm2.B * esk2;
            // var cm1m2 = cm1.A * cm2 - cm1.B * scm2;

            var decrypt = FHE.DecryptBGV(cm1m2, pm, sk, t);
            Console.WriteLine($"m1 * m2:{m1m2}");
            Console.WriteLine($"decrypt:{decrypt}");
            Console.WriteLine();
            if (!(decrypt - m1m2).IsZero())
                throw new($"step[{i}]");
        }
    }

    static void TrackingErrors(int n, int t0, int q0, int size = 100)
    {
        var (pm, sk, t, q, pk, rlk) = FHE.KeyGenBGV(n, t0, q0);
        FHE.Show(pm, sk, t, q, pk, rlk);

        Console.WriteLine($"Size {size}");
        var seqMsg = size.SeqLazy().Select(_ => FHE.GenUnif(n, t0)).ToArray();
        var seqCipher = seqMsg.Select(m => FHE.EncryptBGV(m, pm, t, q, pk)).ToArray();
        var seqAddMsg = new List<Rq>();
        var seqMulMsg = new List<Rq>();
        var seqAddCipher = new List<BGVCipher>();
        var seqMulCipher = new List<BGVCipher>();
        for (int i = 0; i < size; i++)
        {
            if (i == 0)
            {
                seqAddMsg.Add(seqMsg[0]);
                seqMulMsg.Add(seqMsg[0]);
                seqAddCipher.Add(seqCipher[0]);
                seqMulCipher.Add(seqCipher[0]);
                continue;
            }

            seqAddMsg.Add((seqAddMsg.Last() + seqMsg[i]).CoefsMod(t));
            seqMulMsg.Add((seqMulMsg.Last() * seqMsg[i]).ResMod(pm, t));
            seqAddCipher.Add(FHE.AddBGV(seqAddCipher.Last(), seqCipher[i], q));
            seqMulCipher.Add(FHE.MulRelinBGV(seqMulCipher.Last(), seqCipher[i], pm, q, rlk));
        }

        var seqAddDecrypt = seqAddCipher.Select(ct => FHE.DecryptBGV(ct, pm, sk, t)).ToArray();
        var setAdd = seqAddDecrypt.Select((c, i) => (c, i)).Where(e => !e.c.Equals(seqAddMsg[e.i])).SingleOrEmpty();
        var idxAdd = setAdd.Length == 0 ? "Empty" : $"At idx:{setAdd[0].i}";
        Console.WriteLine($"Add Cumulative first Error: {idxAdd}");

        var seqMulDecrypt = seqMulCipher.Select(ct => FHE.DecryptBGV(ct, pm, sk, t)).ToArray();
        var setMul = seqMulDecrypt.Select((c, i) => (c, i)).Where(e => !e.c.Equals(seqMulMsg[e.i])).SingleOrEmpty();
        var idxMul = setMul.Length == 0 ? "Empty" : $"At idx:{setMul[0].i}";
        Console.WriteLine($"Mul Cumulative first Error: {idxMul}");

        Console.WriteLine();
    }

    public static void TestTrackingErrors()
    {
        TrackingErrors(16, 7, 35, 500);
        // TrackingErrors(16, 8, 8, 50);
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
    
    public static void TestSwitchKeys()
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
    
    public static void TestBlindRotate()
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
            Console.WriteLine($"f*X^{u,-2}      :{expected}");
            Console.WriteLine($"f*X^{u,-2}      :[{expected.CoefsExtended(n - 1).Glue(", ", fmt)}]");
            var factor = ((int)q.Num).SeqLazy(1).First(k0 => (k0 * actual).CoefsMod(q).Equals(expected));
            Console.WriteLine($"factor      :{factor}");
            Console.WriteLine();
            set.Add(factor);
        }

        Console.WriteLine($"Factors:{set.Order().Glue(", ", fmt)}");
    }

    public static void TestBootstrapping()
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

    public static void TestRgswMul()
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

    public static void TestHEMul2()
    {
        var l = 4;
        var n = 1 << l;
        var t0 = n;
        var q0 = n * n;
        var (pm, sk, t, q, pk, rlk) = FHE.KeyGenBGV(n, t0, q0);
        FHE.Show(pm, sk, t, q, pk, rlk);
        
        var es = FHE.EncryptBGV(sk, pm, t, q, pk);
        var es2 = FHE.EncryptBGV((sk * sk).ResMod(pm, t), pm, t, q, pk);

        es.Show("Enc(SK)");
        es2.Show("Enc(SK^2)");
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

    public static void TestAll()
    {
        FirstBGV();
        HEAddBGV();
        HEMulBGV();
        TestTrackingErrors();
        HEEvalAuto();
        TestSwitchKeys();
        TestRgswMul();
        TestBlindRotate();
        TestBootstrapping();
    }
}
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
            var m1m2 = (m1 * m2).ResMod(pm).CoefsMod(t);
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
            seqMulMsg.Add((seqMulMsg.Last() * seqMsg[i]).ResMod(pm).CoefsMod(t));
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
                var mk = m.Substitute(xk).ResMod(pm).CoefsMod(t);
                Console.WriteLine($"m   :{m}");
                Console.WriteLine($"m^{k,2}:{mk}");

                var cm = FHE.EncryptBGV(m, pm, t, q, pk);
                var ck = FHE.AutoMorphBGV(cm, k, pm, t, q, pk, rlk, ak);
                var dk = FHE.DecryptBGV(ck, pm, sk, t);
                Console.WriteLine($"    :{dk}");
                if (!dk.Equals(mk))
                    throw new($"k:{k} step[{i}] {dk.Div(mk)}");

                Console.WriteLine();
            }
        }
    }
    
    public static void TestKeysExchange()
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

                var dm = FHE.DecryptBGV(FHE.SwitchKeysBGV(cm, pm, t, q, pk2, rlk2, s2s1), pm, sk2, t);
                Console.WriteLine($"  :{dm}");

                if (!dm.Equals(m))
                    throw new($"step[{i}]");
            }

            {
                var m = FHE.GenUnif(n, t0);
                Console.WriteLine($"m :{m}");
                var cm = FHE.EncryptBGV(m, pm, t, q, pk2);

                var dm = FHE.DecryptBGV(FHE.SwitchKeysBGV(cm, pm, t, q, pk1, rlk1, s1s2), pm, sk1, t);
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
        var t0 = 2 * n;
        var q0 = t0 * 3.Pow(3);
        var (pm, sk, t, q, pk, rlk) = FHE.KeyGenBGV(n, t0, q0);
        FHE.Show(pm, sk, t, q, pk, rlk);

        var Q = 25 * t;
        var delta = q / t;
        Console.WriteLine(new { t, Q, delta });

        var brk = FHE.BRKBGV(pm, sk, t, q, pk);

        var x = pm.X;
        var c = n / 4;
        var f = (2 * c + 1).SeqLazy(-c).Select(j => delta * j * FHE.XpowA(j, pm, t))
            .Aggregate(x.Zero, (acc, v) => acc + v).ResMod(pm).CoefsMod(t);
        var nbDigits = $"{Q}".Length;
        var fmt = $"{{0,{nbDigits}}}";
        var set = new HashSet<int>();
        var s = n.SeqLazy().Select(i => sk[i]).ToArray();

        for (int k = 0; k < 50; ++k)
        {
            var ai = new Rational(IntExt.Rng.Next(1, t0));
            var bi = n.SeqLazy().Select(_ => IntExt.Rng.Next(t0) * ai.One).ToArray();
            var acc = FHE.BlindRotateBGV((ai, bi), f, pm, t, q, pk, rlk, brk);
            var actual = FHE.DecryptBGV(acc, pm, sk, t);

            Console.WriteLine($"f           :{f}");
            Console.WriteLine($"f           :[{f.CoefsExtended(n - 1).Glue(", ", fmt)}]");
            Console.WriteLine($"ai          :{ai}");
            Console.WriteLine($"bi          :[{bi.Glue(", ", fmt)}]");
            Console.WriteLine($"blind rotate:{actual}");
            Console.WriteLine($"            :[{actual.CoefsExtended(n - 1).Glue(", ", fmt)}]");

            // Testing result
            var u = (ai - FHE.InnerProd(bi, s, t)).Mod(t);
            var expected = (x.Pow((int)u.Num) * f).ResMod(pm).CoefsMod(t);
            Console.WriteLine($"u= a - <b,s>:{u}");
            Console.WriteLine($"f*X^{u,-2}      :{expected}");
            Console.WriteLine($"f*X^{u,-2}      :[{expected.CoefsExtended(n - 1).Glue(", ", fmt)}]");
            var factor = t0.SeqLazy(1).First(k0 => (k0 * actual).CoefsMod(t).Equals(expected));
            Console.WriteLine($"factor      :{factor}");
            Console.WriteLine();
            set.Add(factor);
        }

        Console.WriteLine($"Factors:{set.Order().Glue(", ", fmt)}");
    }

    public static void TestAll()
    {
        FirstBGV();
        HEAddBGV();
        HEMulBGV();
        TestTrackingErrors();
        HEEvalAuto();
        TestKeysExchange();
    }
}
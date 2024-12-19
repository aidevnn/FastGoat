using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Integers;

namespace FastGoat.UserGroup.Lattice;

public static class BGVtests
{
    static BGVtests()
    {
        Ring.DisplayPolynomial = MonomDisplay.StarCaret;
    }

    public static void FirstBGV()
    {
        var (n, q) = (16, 424242);
        var p = n / 2 - 1;
        var bgv = new BGV(n, p, q);
        bgv.Show();
        for (int i = 0; i < 100; i++)
        {
            var mi = BGVPublic.GenUnif(n, p);
            var ct = bgv.BGVpublic.Encrypt(mi);
            var mf = bgv.Decrypt(ct);
            Console.WriteLine(mi);
            Console.WriteLine(mf);
            Console.WriteLine();
            if (!(mf - mi).IsZero())
                throw new($"step[{i}]");
        }
    }

    public static void HEAddBGV()
    {
        var (n, p, q) = (16, 8, 2.Pow(6));
        var bgv = new BGV(n, p, q);
        bgv.Show();
        var pub = bgv.BGVpublic;

        for (int i = 0; i < 100; ++i)
        {
            var m1 = BGVPublic.GenUnif(n, p);
            var m2 = BGVPublic.GenUnif(n, p);
            var m1m2 = (m1 + m2).ResMod(pub.PM).CoefsMod(bgv.P);
            var ct1 = pub.Encrypt(m1);
            var ct2 = pub.Encrypt(m2);

            var decrypt = bgv.Decrypt(pub.Add(ct1, ct2));
            Console.WriteLine($"m1 + m2:{m1m2}");
            Console.WriteLine($"decrypt:{decrypt}");
            Console.WriteLine();
            if (!(decrypt - m1m2).IsZero())
                throw new($"step[{i}]");
        }
    }

    public static void HEMulBGV()
    {
        var (n, p, q) = (16, 8, 2.Pow(6));
        var bgv = new BGV(n, p, q);
        bgv.Show();
        var pub = bgv.BGVpublic;

        for (int i = 0; i < 100; ++i)
        {
            var m1 = BGVPublic.GenUnif(n, p);
            var m2 = BGVPublic.GenUnif(n, p);
            var m1m2 = (m1 * m2).ResMod(pub.PM).CoefsMod(bgv.P);
            var ct1 = pub.Encrypt(m1);
            var ct2 = pub.Encrypt(m2);

            var decrypt = bgv.Decrypt(pub.Mul(ct1, ct2));
            Console.WriteLine($"m1 * m2:{m1m2}");
            Console.WriteLine($"decrypt:{decrypt}");
            Console.WriteLine();
            if (!(decrypt - m1m2).IsZero())
                throw new($"step[{i}]");
        }
    }

    static void TrackingErrors(int n, int p, int q, int size = 100)
    {
        var bgv = new BGV(n, p, q);
        bgv.Show(showKeys: false);
        var pub = bgv.BGVpublic;

        Console.WriteLine($"Size {size}");
        var seqMsg = size.SeqLazy().Select(_ => BGVPublic.GenUnif(n, p)).ToArray();
        var seqCipher = seqMsg.Select(pub.Encrypt).ToArray();
        var seqAddMsg = new List<KPoly<Rational>>();
        var seqMulMsg = new List<KPoly<Rational>>();
        var seqAddCipher = new List<(KPoly<Rational> ct0, KPoly<Rational> ct1)>();
        var seqMulCipher = new List<(KPoly<Rational> ct0, KPoly<Rational> ct1)>();
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

            seqAddMsg.Add((seqAddMsg.Last() + seqMsg[i]).CoefsMod(bgv.P));
            seqMulMsg.Add((seqMulMsg.Last() * seqMsg[i]).ResMod(pub.PM).CoefsMod(bgv.P));
            seqAddCipher.Add(pub.Add(seqAddCipher.Last(), seqCipher[i]));
            seqMulCipher.Add(pub.Mul(seqMulCipher.Last(), seqCipher[i]));
        }

        var seqAddDecrypt = seqAddCipher.Select(bgv.Decrypt).ToArray();
        var setAdd = seqAddDecrypt.Select((c, i) => (c, i)).Where(e => !e.c.Equals(seqAddMsg[e.i])).SingleOrEmpty();
        var idxAdd = setAdd.Length == 0 ? "Empty" : $"At idx:{setAdd[0].i}";
        Console.WriteLine($"Add Cumulative first Error: {idxAdd}");

        var seqMulDecrypt = seqMulCipher.Select(bgv.Decrypt).ToArray();
        var setMul = seqMulDecrypt.Select((c, i) => (c, i)).Where(e => !e.c.Equals(seqMulMsg[e.i])).SingleOrEmpty();
        var idxMul = setMul.Length == 0 ? "Empty" : $"At idx:{setMul[0].i}";
        Console.WriteLine($"Mul Cumulative first Error: {idxMul}");

        Console.WriteLine();
    }

    public static void TestTrackingErrors()
    {
        TrackingErrors(16, 7, 35, 50);
        // TrackingErrors(16, 8, 8, 50);
    }

    public static void HEEvalAuto()
    {
        var (n, p, q) = (16, 8, 2.Pow(6));
        var bgv = new BGV(n, p, q);
        var pub = bgv.BGVpublic;
        bgv.Show();
        var (N, P, Q, PM) = pub.Params_NPQ_PM;
        var x = PM.X;

        foreach (var k in IntExt.Coprimes(2 * N))
        {
            var xk = x.Pow(k);
            for (var i = 0; i < 50; ++i)
            {
                var m = BGVPublic.GenUnif(n, p);
                var mk = m.Substitute(xk).ResMod(PM).CoefsMod(P);
                Console.WriteLine($"m   :{m}");
                Console.WriteLine($"m^{k,2}:{mk}");

                var cm = pub.Encrypt(m);
                var ck = pub.EvalAuto(cm, k);
                var dk = bgv.Decrypt(ck);
                Console.WriteLine($"    :{dk}");
                if (!dk.Equals(mk))
                    throw new($"k:{k} step[{i}] {dk.Div(mk)}");

                Console.WriteLine();
            }
        }
    }

    public static void TestRGSW()
    {
        var (n, p, q) = (16, 8, 2.Pow(10));
        var bgv = new BGV(n, p, q);
        bgv.Show();

        var pub = bgv.BGVpublic;
        var expected = bgv.SK.CoefsMod(pub.P);
        var actual = bgv.Decrypt((pub.PM.Zero, -pub.PM.One)); // BGV.Enc(s) = (0, -1) = RGSW0
        Console.WriteLine($"ct0:{pub.PM.Zero}");
        Console.WriteLine($"ct1:{pub.PM.One}");
        Console.WriteLine($"s      :{expected}");
        Console.WriteLine($"decrypt:{actual}");
        Console.WriteLine();
        if (!expected.Equals(actual))
            throw new();

        for (int k = 0; k < 10; ++k)
        {
            var r0 = pub.RGSW(pub.PM.One);
            actual = bgv.Decrypt(r0);

            Console.WriteLine($"ct0:{r0.ct0}");
            Console.WriteLine($"ct1:{r0.ct1}");
            Console.WriteLine($"s      :{expected}");
            Console.WriteLine($"decrypt:{actual}");
            Console.WriteLine();
            if (!expected.Equals(actual))
                throw new();
        }

        for (int k = 0; k < 10; ++k)
        {
            var m = BGVPublic.GenUnif(n, p);
            expected = (bgv.SK * m).ResMod(pub.PM).CoefsMod(pub.P);

            var m0 = pub.RGSW(m);
            actual = bgv.Decrypt(m0); // BGV.Enc(s * m) = BGV.Enc(0) + (0, -m)

            Console.WriteLine($"ct0:{m0.ct0}");
            Console.WriteLine($"ct1:{m0.ct1}");
            Console.WriteLine($"s * m  :{expected}");
            Console.WriteLine($"decrypt:{actual}");
            Console.WriteLine();
            if (!expected.Equals(actual))
                throw new();
        }
    }

    public static void TestAll()
    {
        FirstBGV();
        HEAddBGV();
        HEMulBGV();
        TestTrackingErrors();
        HEEvalAuto();
        TestRGSW();
    }
}
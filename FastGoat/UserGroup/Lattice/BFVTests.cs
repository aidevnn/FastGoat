using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Integers;

namespace FastGoat.UserGroup.Lattice;

public static class BFVTests
{
    public static void firstBFV()
    {
        // IntExt.RngSeed(87456);
        Ring.DisplayPolynomial = MonomDisplay.StarCaret;
        var (n, q) = (34, 424242);
        var bfv = new BFV(n, q);
        bfv.Show();
        for (int i = 0; i < 100; i++)
        {
            var mi = BFVPublic.GenUnif(n, bfv.P);
            var ct = bfv.BFVpublic.Encrypt(mi);
            var mf = bfv.Decrypt(ct);
            Console.WriteLine(mi);
            Console.WriteLine(mf);
            Console.WriteLine();
            if (!(mf - mi).IsZero())
                throw new($"step[{i}]");
        }
    }

    public static void HEAddBFV()
    {
        // IntExt.RngSeed(87456);
        Ring.DisplayPolynomial = MonomDisplay.StarCaret;
        var (n, q) = (16, 424242);
        var bfv = new BFV(n, q);
        bfv.Show();
        var bfvPub = bfv.BFVpublic;

        for (int i = 0; i < 100; ++i)
        {
            var m1 = BFVPublic.GenUnif(n, bfv.P);
            var m2 = BFVPublic.GenUnif(n, bfv.P);
            var m1m2 = BFVPublic.CoefsMod((m1 + m2).Div(bfvPub.PM).rem, bfv.P);
            var ct1 = bfvPub.Encrypt(m1);
            var ct2 = bfvPub.Encrypt(m2);

            var decrypt = bfv.Decrypt(bfvPub.Add(ct1, ct2));
            Console.WriteLine($"m1 + m2:{m1m2}");
            Console.WriteLine($"decrypt:{decrypt}");
            Console.WriteLine();
            if (!(decrypt - m1m2).IsZero())
                throw new($"step[{i}]");
        }
    }

    public static void HEMulBFV()
    {
        // IntExt.RngSeed(87456);
        Ring.DisplayPolynomial = MonomDisplay.StarCaret;
        var (n, p, q) = (16, 7, 165000);
        var bfv = new BFV(n, p, q);
        bfv.Show();
        var bfvPub = bfv.BFVpublic;

        for (int i = 0; i < 1000; i++)
        {
            var m1 = BFVPublic.GenUnif(n, bfv.P);
            var m2 = BFVPublic.GenUnif(n, bfv.P);
            var m1m2 = BFVPublic.CoefsMod((m1 * m2).Div(bfvPub.PM).rem, bfv.P);
            var ct1 = bfvPub.Encrypt(m1);
            var ct2 = bfvPub.Encrypt(m2);

            var decrypt = bfv.Decrypt(bfvPub.Mul(ct1, ct2));
            Console.WriteLine($"m1 * m2:{m1m2}");
            Console.WriteLine($"decrypt:{decrypt}");
            Console.WriteLine();
            if (!(decrypt - m1m2).IsZero())
                throw new($"step[{i}]");
        }
    }

    static void TrackingErrors(int n, int p, int q, int size = 100)
    {
        var bfv = new BFV(n, p, q);
        bfv.Show(showKeys: false);
        var bfvPub = bfv.BFVpublic;

        Console.WriteLine($"Size {size}");
        var seqMsg = size.SeqLazy().Select(_ => BFVPublic.GenUnif(n, bfv.P)).ToArray();
        var seqCipher = seqMsg.Select(bfvPub.Encrypt).ToArray();
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

            seqAddMsg.Add(BFVPublic.CoefsMod(seqAddMsg.Last() + seqMsg[i], p));
            seqMulMsg.Add(BFVPublic.CoefsMod((seqAddMsg.Last() * seqMsg[i]).Div(bfvPub.PM).rem, p));
            seqAddCipher.Add(bfvPub.Add(seqAddCipher.Last(), seqCipher[i]));
            seqMulCipher.Add(bfvPub.Mul(seqAddCipher.Last(), seqCipher[i]));
        }

        var seqAddDecrypt = seqAddCipher.Select(bfv.Decrypt).ToArray();
        var setAdd = seqAddDecrypt.Select((c, i) => (c, i)).Where(e => !e.c.Equals(seqAddMsg[e.i])).SingleOrEmpty();
        var idxAdd = setAdd.Length == 0 ? "Empty" : $"At idx:{setAdd[0].i}";
        Console.WriteLine($"Add Cumulative first Error: {idxAdd}");

        var seqMulDecrypt = seqMulCipher.Select(bfv.Decrypt).ToArray();
        var setMul = seqMulDecrypt.Select((c, i) => (c, i)).Where(e => !e.c.Equals(seqMulMsg[e.i])).SingleOrEmpty();
        var idxMul = setMul.Length == 0 ? "Empty" : $"At idx:{setMul[0].i}";
        Console.WriteLine($"Mul Cumulative first Error: {idxMul}");

        Console.WriteLine();
    }

    public static void TrackingCumulativeErrors()
    {
        Ring.DisplayPolynomial = MonomDisplay.StarCaret;
        TrackingErrors(16, 2, 2.Pow(13));
        TrackingErrors(16, 2, 2.Pow(17));
    }
    // N = 8 Q = 1024 T = 128 P = 2
    // Size 100
    // Add Cumulative first Error: Empty
    // Mul Cumulative first Error: At idx:2
    // 
    // N = 8 Q = 16384 T = 512 P = 2
    // Size 100
    // Add Cumulative first Error: Empty
    // Mul Cumulative first Error: Empty
    // 

}
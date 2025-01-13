using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Lattice;

namespace FastGoat.Examples;

public static class LearningWithErrors
{
    static LearningWithErrors()
    {
        Ring.DisplayPolynomial = MonomDisplay.StarCaret;
        LipsumSentences = Lipsum.BySentences;
        LipsumParagraphes = Lipsum.ByParagraphes;

        Console.WriteLine("################################################");
        Console.WriteLine("#                                              #");
        Console.WriteLine("#           Learning With Errors               #");
        Console.WriteLine("#        for crafting purposes only            #");
        Console.WriteLine("#                                              #");
        Console.WriteLine("################################################");
        Console.WriteLine();
    }

    private static string[] LipsumSentences { get; }
    private static string[] LipsumParagraphes { get; }

    #region Char to Bit
    // 16-bit default C#-char
    const int charBit = 8;

    static IEnumerable<int> String2Bin(string s)
    {
        foreach (var c in s)
        {
            var c0 = Convert.ToInt32(c);
            for (int i = 0; i < charBit; ++i)
            {
                var k = c0 & 1;
                yield return k;
                c0 = c0 >> 1;
            }
        }
    }

    static string Bin2String(IEnumerable<int> bin)
    {
        var s = "";
        foreach (var l in bin.Chunk(charBit))
        {
            var sum = l.Select((c, i) => (c, i)).Sum(e => e.c << e.i);
            s += (char)sum;
        }

        return s;
    }
    
    #endregion

    static void RunLWERegev(Regev lwe, string text, bool showCipher = false, bool showBinary = true)
    {
        Console.WriteLine(text);
        var seq = String2Bin(text).ToArray();
        var seqCiphers = seq.Select(b => (b, cipher: lwe.EncryptBit(b))).ToArray();

        if (showCipher)
            seqCiphers.Println(
                l => $"{l.b} => {l.cipher}",
                $"Cyphers text:{text}"
            );

        var seqDecrypt = seqCiphers.Select(e => lwe.DecryptBit(e.cipher)).ToArray();
        var text2 = Bin2String(seqDecrypt);
        if (showBinary)
        {
            Console.WriteLine($"seqInput  :[{seq.Glue()}]");
            Console.WriteLine($"seqDecrypt:[{seqDecrypt.Glue()}]");
            Console.WriteLine(text2);
        }

        if (string.Equals(text, text2))
            Console.WriteLine("    SUCCESS");
        else
            Console.WriteLine("    FAIL");

        Console.WriteLine();
    }

    static void RunRLWE(RLWE rlwe, string text, bool showBinary = true)
    {
        Console.WriteLine(text);
        
        var seq = String2Bin(text).ToArray();
        var seqCiphers = rlwe.Encrypt(seq);

        var seqDecrypt = rlwe.Decrypt(seqCiphers);
        var text2 = Bin2String(seqDecrypt);
        
        if (showBinary)
        {
            Console.WriteLine($"seqInput  :[{seq.Glue()}]");
            Console.WriteLine($"seqDecrypt:[{seqDecrypt.Glue()}]");
            Console.WriteLine(text2);
        }

        if (string.Equals(text, text2))
            Console.WriteLine("    SUCCESS");
        else
        {
            Console.WriteLine("    FAIL");
            Console.Beep();
        }

        Console.WriteLine();
    }

    public static void Example1Regev()
    {
        GlobalStopWatch.Restart();
        // Warning : Weak parameters
        var reg = new Regev(32);
        reg.Show();

        RunLWERegev(reg, "hello world lwe");
        RunLWERegev(reg, "Hello World LWE");
        RunLWERegev(reg, "AAA+", showCipher: true);

        for (int i = 0; i < 10; i++)
            RunLWERegev(reg, DistributionExt.Dice(LipsumSentences));
        
        // long text
        RunLWERegev(reg, DistributionExt.Dice(LipsumParagraphes), showBinary: false);
        GlobalStopWatch.Show();
    }

    public static void Example2RLWE()
    {
        GlobalStopWatch.Restart();
        // Warning : Weak parameters
        var rlwe = new RLWE(32);
        rlwe.Show();

        RunRLWE(rlwe, "hello world lwe");
        RunRLWE(rlwe, "Hello World LWE");
        RunRLWE(rlwe, "AAA+");

        for (int i = 0; i < 10; i++)
            RunRLWE(rlwe, DistributionExt.Dice(LipsumSentences));
        
        // long text
        RunRLWE(rlwe, DistributionExt.Dice(LipsumParagraphes), showBinary: false);
        GlobalStopWatch.Show();
    }

    public static void Example3AdditionMultiplication()
    {
        // Warning : Weak parameters
        // RLWE N=16=2^4, Φ(N)=8 PM=x^8 + 1 t=97 q=9797=97*101
        var rlwe = new RLWE(16);
        var (n, pm, sk, t, q, pk, rlk) = rlwe;
        rlwe.Show();

        var m1 = FHE.GenUnif(n, t);
        var e1 = FHE.EncryptBGV(m1, pm, t, q, pk);
        var m2 = FHE.GenUnif(n, rlwe.T);
        var e2 = FHE.EncryptBGV(m2, pm, t, q, pk);

        e1.Show($"e1 = Encrypt(m1 = {m1})");
        e2.Show($"e2 = Encrypt(m2 = {m2})");
        Console.WriteLine();

        var add_m1m2 = (m1 + m2).CoefsMod(t);
        var add_e1e2 = (e1 + e2).CoefsMod(q);
        var d_add = FHE.DecryptBGV(add_e1e2, pm, sk, t);
        
        add_e1e2.Show("e1 + e2");
        Console.WriteLine($"m1 + m2          = {add_m1m2}");
        Console.WriteLine($"Decrypt(e1 + e2) = {d_add}");
        Console.WriteLine();
        
        var mul_m1m2 = (m1 * m2).ResMod(pm, t);
        var mul_e1e2 = FHE.MulRelinBGV(e1, e2, pm, q, rlk);
        var d_mul = FHE.DecryptBGV(mul_e1e2, pm, sk, t);

        mul_e1e2.Show("e1 * e2");
        Console.WriteLine($"m1 * m2          = {mul_m1m2}");
        Console.WriteLine($"Decrypt(e1 * e2) = {d_mul}");
    }

    public static void Example4LogicGates()
    {
        // Warning : Weak parameters
        // RLWE N=16=2^4, Φ(N)=8 PM=x^8 + 1 t=97 q=9797=97*101
        var rlwe = new RLWE(16);
        rlwe.Show();

        int[] m0 = [0, 1];
        var e0 = rlwe.Encrypt(m0);
        int[] m1 = [0, 1, 0, 1];
        int[] m2 = [0, 0, 1, 1];
        var e1 = rlwe.Encrypt(m1);
        var e2 = rlwe.Encrypt(m2);

        Console.WriteLine("NOT");
        Console.WriteLine($"   [{m0.Glue()}]");
        Console.WriteLine($" = [{rlwe.Decrypt(rlwe.NOT(e0)).Glue()}]");
        Console.WriteLine();
        
        Console.WriteLine("AND");
        Console.WriteLine($"   [{m1.Glue()}]");
        Console.WriteLine($"   [{m2.Glue()}]");
        Console.WriteLine($" = [{rlwe.Decrypt(rlwe.AND(e1, e2)).Glue()}]");
        Console.WriteLine();
        
        Console.WriteLine("OR");
        Console.WriteLine($"   [{m1.Glue()}]");
        Console.WriteLine($"   [{m2.Glue()}]");
        Console.WriteLine($" = [{rlwe.Decrypt(rlwe.OR(e1, e2)).Glue()}]");
        Console.WriteLine();
    }

    public static void Example5AdditionWithCarry()
    {
        // RLWE N=16=2^4, Φ(N)=8 PM=x^8 + 1 t=97 q=9797=97*101
        var rlwe = new RLWE(16);
        rlwe.Show();

        var bits = 32;
        var fmt = $"{{0,{(int)(bits * double.Log10(2)) + 1}}}";
        string FMT(long a) => string.Format(fmt, a);
        
        var m1 = DistributionExt.DiceSample(bits - 1, [0, 1]).Append(0).ToArray();
        var m2 = DistributionExt.DiceSample(bits - 1, [0, 1]).Append(0).ToArray();
        var a1 = Convert.ToInt64(m1.Reverse().Glue(), 2);
        var a2 = Convert.ToInt64(m2.Reverse().Glue(), 2);
        var e1 = rlwe.Encrypt(m1);
        var e2 = rlwe.Encrypt(m2);

        var sum = a1 + a2;
        var add_m1m2 = Convert.ToString(sum, 2);
        add_m1m2 = Enumerable.Repeat(0, bits - add_m1m2.Length).Glue() + add_m1m2;

        var add_e1e2 = rlwe.ADD(e1, e2);
        var d_add = rlwe.Decrypt(add_e1e2);

        Console.WriteLine($"   0b{m1.Reverse().Glue()} = {FMT(a1)}");
        Console.WriteLine($" + 0b{m2.Reverse().Glue()} = {FMT(a2)}");
        Console.WriteLine($" = 0b{d_add.Reverse().Glue()} = {FMT(sum)}");
        Console.WriteLine($"   0b{add_m1m2}");
    }
}
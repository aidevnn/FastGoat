using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.LWE;

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
        // Weak and invalid parameters
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
        // Weak and invalid parameters
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
        // Weak and invalid parameters
        // RLWE N=16=2^4, Φ(N)=8 PM=x^8 + 1 t=97 q=9797=97*101
        var rlwe = new RLWE(16);
        var (n, pm, sk, t, q, pk, rlk) = rlwe;
        rlwe.Show();

        var m1 = RLWE.GenUnif(n, t);
        var e1 = RLWE.EncryptBGV(m1, pm, t, q, pk);
        var m2 = RLWE.GenUnif(n, rlwe.T);
        var e2 = RLWE.EncryptBGV(m2, pm, t, q, pk);

        e1.Show($"e1 = Encrypt(m1 = {m1})");
        e2.Show($"e2 = Encrypt(m2 = {m2})");
        Console.WriteLine();

        var add_m1m2 = (m1 + m2).CoefsMod(t);
        var add_e1e2 = (e1 + e2).CoefsMod(q);
        var d_add = RLWE.DecryptBGV(add_e1e2, pm, sk, t);

        add_e1e2.Show("e1 + e2");
        Console.WriteLine($"m1 + m2          = {add_m1m2}");
        Console.WriteLine($"Decrypt(e1 + e2) = {d_add}");
        Console.WriteLine();

        var mul_m1m2 = (m1 * m2).ResMod(pm, t);
        var mul_e1e2 = RLWE.MulRelinBGV(e1, e2, pm, q, rlk);
        var d_mul = RLWE.DecryptBGV(mul_e1e2, pm, sk, t);

        mul_e1e2.Show("e1 * e2");
        Console.WriteLine($"m1 * m2          = {mul_m1m2}");
        Console.WriteLine($"Decrypt(e1 * e2) = {d_mul}");
    }

    public static void Example4LogicGates()
    {
        // Weak and invalid parameters
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

    public static void Example5HomomorphicAdditionWithCarry()
    {
        // Weak and invalid parameters
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

        var sumi = a1 + a2;
        var add_m1m2 = Convert.ToString(sumi, 2).PadLeft(bits, '0');

        var add_e1e2 = rlwe.ADD(e1, e2);
        var d_add = rlwe.Decrypt(add_e1e2);

        m1.Zip(e1).Take(6).Println($"Encrypt(m1 = 0b{m1.Reverse().Glue()})");
        Console.WriteLine("...");
        m2.Zip(e2).Take(6).Println($"Encrypt(m2 = 0b{m2.Reverse().Glue()})");
        Console.WriteLine("...");
        d_add.Zip(add_e1e2).Take(6).Println($"Decrypt(e1 + e2) = 0b{d_add.Reverse().Glue()}");
        Console.WriteLine("...");

        Console.WriteLine($"   0b{m1.Reverse().Glue()} = {FMT(a1)}");
        Console.WriteLine($" + 0b{m2.Reverse().Glue()} = {FMT(a2)}");
        Console.WriteLine($" = 0b{d_add.Reverse().Glue()} = {FMT(sumi)}");
        Console.WriteLine($"   0b{add_m1m2}");

        var sumf = Convert.ToInt64(d_add.Reverse().Glue(), 2);
        if (sumi != sumf)
            throw new();
    }

    public static void Example6HomomorphicMultiplicationWithCarry()
    {
        // Weak and invalid parameters
        // RLWE N=8=2^3, Φ(N)=4 PM=x^4 + 1 t=17 q=323=17*19
        var rlwe = new RLWE(8);
        rlwe.Show();

        var bits = 32;
        var fmt = $"{{0,{(int)(2 * bits * double.Log10(2)) + 1}}}";
        string FMT(long a) => string.Format(fmt, a);

        var m1 = DistributionExt.DiceSample(bits - 1, [0, 1]).Append(0).ToArray();
        var m2 = DistributionExt.DiceSample(bits - 1, [0, 1]).Append(0).ToArray();
        var a1 = Convert.ToInt64(m1.Reverse().Glue(), 2);
        var a2 = Convert.ToInt64(m2.Reverse().Glue(), 2);
        var e1 = rlwe.Encrypt(m1);
        var e2 = rlwe.Encrypt(m2);

        var prodi = a1 * a2;
        var mult_m1m2 = Convert.ToString(prodi, 2).PadLeft(bits * 2, '0');

        var mult_e1e2 = rlwe.MULT(e1, e2);
        var d_mult = rlwe.Decrypt(mult_e1e2);

        m1.Zip(e1).Take(6).Println($"Encrypt(m1 = 0b{m1.Reverse().Glue()})");
        Console.WriteLine("...");
        m2.Zip(e2).Take(6).Println($"Encrypt(m2 = 0b{m2.Reverse().Glue()})");
        Console.WriteLine("...");
        d_mult.Zip(mult_e1e2).Take(6).Println($"Decrypt(e1 * e2) = 0b{d_mult.Reverse().Glue()}");
        Console.WriteLine("...");

        var zeros = Enumerable.Repeat(0, bits).ToArray();
        Console.WriteLine($"   0b{m1.Concat(zeros).Reverse().Glue()} = {FMT(a1)}");
        Console.WriteLine($" * 0b{m2.Concat(zeros).Reverse().Glue()} = {FMT(a2)}");
        Console.WriteLine($" = 0b{d_mult.Reverse().Glue()} = {FMT(prodi)}");
        Console.WriteLine($"   0b{mult_m1m2}");
        Console.WriteLine();

        var prodf = Convert.ToInt64(d_mult.Reverse().Glue(), 2);
        if (prodi != prodf)
            throw new();
    }

    public static void Example7Regev2RLWE()
    {
        var (reg, rlwe) = Regev.SetupRLWE(16);
        reg.Show();
        rlwe.Show();

        for (int k = 0; k < 10; ++k)
        {
            var m = DistributionExt.DiceSample(reg.N, [0, 1]).ToArray();

            var c11 = reg.Encrypt(m);
            var c21 = rlwe.FromRegevCipher(c11);
            var d1 = rlwe.Decrypt(c21);

            var c12 = rlwe.Encrypt(m);
            var c22 = rlwe.ToRegevCipher(c12);
            var d2 = reg.Decrypt(c22);

            Console.WriteLine($"m :[0b{m.Reverse().Glue()}]");
            Console.WriteLine($"d1:[0b{d1.Reverse().Glue()}] Regev to RLWE");
            Console.WriteLine($"d2:[0b{d2.Reverse().Glue()}] RLWE  to Regev");
            Console.WriteLine();

            if (!m.SequenceEqual(d1) || !m.SequenceEqual(d2))
                throw new($"step[{k}]");
        }
    }

    public static void Example8TrackingErrors()
    {
        // Weak and invalid parameters
        // RLWE N=16=2^4, Φ(N)=8 PM=x^8 + 1 t=97 q=9797=97*101
        var (reg, rlwe) = Regev.SetupRLWE(16);
        var t = rlwe.T;

        var m1 = DistributionExt.DiceSample(reg.N, [0, 1]).ToArray();
        var m2 = DistributionExt.DiceSample(reg.N, [0, 1]).ToArray();
        var c1a = rlwe.Encrypt(m1);
        var c2a = rlwe.Encrypt(m2);
        var xa = rlwe.XOR(
            rlwe.XOR(c1a.Take(rlwe.n / 2).ToArray(), c1a.Skip(rlwe.n / 2).ToArray()),
            rlwe.XOR(c2a.Take(rlwe.n / 2).ToArray(), c2a.Skip(rlwe.n / 2).ToArray())
        );
        var eia = rlwe.Errors(c1a).Concat(rlwe.Errors(c2a)).Select(e => e[0].Signed(t).Absolute).Max();
        var efa = xa.Select(e => rlwe.Errors(e)[0].Signed(t).Absolute).Max();

        var c1b = rlwe.FromRegevCipher(reg.Encrypt(m1));
        var c2b = rlwe.FromRegevCipher(reg.Encrypt(m2));
        var xb = rlwe.XOR(
            rlwe.XOR(c1b.Take(rlwe.n / 2).ToArray(), c1b.Skip(rlwe.n / 2).ToArray()),
            rlwe.XOR(c2b.Take(rlwe.n / 2).ToArray(), c2b.Skip(rlwe.n / 2).ToArray())
        );
        var eib = rlwe.Errors(c1b).Concat(rlwe.Errors(c2b)).Select(e => e[0].Signed(t).Absolute).Max();
        var efb = xb.Select(e => rlwe.Errors(e)[0].Signed(t).Absolute).Max();

        Console.WriteLine($"{reg.Params}");
        Console.WriteLine($"{rlwe.Params}");
        Console.WriteLine();
        Console.WriteLine($"RLWE only             Errors: {eia,3} --> {efa,3}");
        Console.WriteLine($"RLWE from RegevCipher Errors: {eib,3} --> {efb,3}");
        // Regev N:16   P:577    M:171    A:0.0156 A*P:9.0156 P/M:3.3743
        // RLWE N=32=2^5, Φ(N)=16 PM=x^16 + 1 t=577 q=338699=577*587
        // 
        // RLWE only             Errors:   0 -->   0
        // RLWE from RegevCipher Errors:  17 --> 135
        // 
    }

    public static void Example9WrongParameters()
    {
        // Warning : Weak and invalid parameters
        // RLWE N=16=2^4, Φ(N)=8 PM=x^8 + 1 t=97 q=9797=97*101
        for (int k = 2; k < 8; k++)
        {
            var rlwe = new RLWE(1 << k); // ciphertext modulus q is multiple of plaintext modulus t
            var (n, pm, sk, t, q, pk, rlk) = rlwe;
            var t0 = (int)t.Num;
    
            var pka = pk.A.ToZnPoly(t0);
            var pkb = pk.B.ToZnPoly(t0);
    
            var a = FG.EPoly(pka.X.Pow(n) + 1, 'a');
            Console.WriteLine(rlwe.Params);
            pk.CoefsMod(t).Show($"pk mod {t}");
            var sk1 = sk.ToZnPoly(t0).Substitute(a);
            var sk2 = pka.Substitute(a) / pkb.Substitute(a);
            Console.WriteLine($"sk = {sk}");
            Console.WriteLine($"   = {sk1}");
            Console.WriteLine($"   = {sk2}");
            Console.WriteLine($"Invalid:{sk1.Equals(sk2)}"); // always true
            Console.WriteLine();
        }
    }

}
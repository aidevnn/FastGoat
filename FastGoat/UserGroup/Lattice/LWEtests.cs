using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Integers;
using static FastGoat.Commons.IntExt;

namespace FastGoat.UserGroup.Lattice;

public static class LWEtests
{
    static LWEtests()
    {
        LipsumSentences = Lipsum.BySentences;
        LipsumParagraphes = Lipsum.ByParagraphes;
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

    static void RunLWE(LWE lwe, string text, bool showCipher = false, bool showBinary = true)
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
        {
            Console.WriteLine("    FAIL");
            Console.Beep();
        }

        Console.WriteLine();
    }

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
        {
            Console.WriteLine("    FAIL");
            Console.Beep();
        }

        Console.WriteLine();
    }

    public static void TestEncryptDecryptLWE()
    {
        for (int n = 4; n < 31; ++n)
        {
            var lwe = new LWE(n);
            lwe.Show();

            RunLWE(lwe, "hello world lwe");
            RunLWE(lwe, "Hello World LWE");
            RunLWE(lwe, "AAA+", showCipher: true);

            for (int i = 0; i < 10; i++)
                RunLWE(lwe, DistributionExt.Dice(LipsumSentences));
        
            // long text
            RunLWE(lwe, DistributionExt.Dice(LipsumParagraphes), showBinary: false);
            Console.WriteLine();
        }
    }

    public static void ExpandProduct()
    {
        var n = 4;
        var ind = new Indeterminates<Xi>(MonomOrder.GrLex, n.SeqLazy().Select(i => new Xi($"s{i}")).ToArray());
        ind.ExtendAppend(n.SeqLazy().Select(i => new Xi($"a{i}")).ToArray());
        ind.ExtendAppend(n.SeqLazy().Select(i => new Xi($"b{i}")).ToArray());
        var xis = new Polynomial<Rational, Xi>(ind, Rational.KZero());
        var sk = xis.Variables.Take(n).ToVec();
        var ak = xis.Variables.Skip(n).Take(n).ToVec();
        var bk = xis.Variables.Skip(2 * n).Take(n).ToVec();

        var prod1 = (ak * sk).Sum() * (bk * sk).Sum();

        Console.WriteLine($"ak:{ak}");
        Console.WriteLine($"bk:{bk}");
        Console.WriteLine($"sk:{sk}");
        var akbk = n.Range().Grid2D().Select(e => ak[e.t1] * bk[e.t2]).ToVec();
        var sksk = n.Range().Grid2D().Select(e => sk[e.t1] * sk[e.t2]).ToVec();
        var prod2 = (akbk * sksk).Sum();

        Console.WriteLine(prod1);
        Console.WriteLine(akbk);
        Console.WriteLine(sksk);
        Console.WriteLine(prod2);
        Console.WriteLine(prod1 - prod2);
        Console.WriteLine();
    }

    public static void TestHELogicGates()
    {
        for(int n = 4; n < 31; n += 1)
        {
            var lwe = new LWE(n);
            lwe.Show();

            var nbTrials = 100;
            for (int l = 0; l < nbTrials; l++)
            {
                var bit = Rng.Next(2);
                var enc = lwe.EncryptBit(bit);
                var d = lwe.DecryptBit(enc);
                if (l < 10)
                    Console.WriteLine($"b:{bit} {enc} => {d}");

                if (d != bit)
                    throw new($"step[{l}] N:{n} Q:{lwe.Q}");
            }

            for (int l = 0; l < nbTrials; ++l)
            {
                var table = 2.Range().Grid2D().ToDictionary(e => e, e => (lwe.EncryptBit(e.t1), lwe.EncryptBit(e.t2)));
                var tableNot = 2.Range()
                    .Select(e => (e, 1 - e, lwe.DecryptBit(CipherLWE.Not(lwe.EncryptBit(e))))).ToArray();
                var tableXor = table.Select(e =>
                        (e.Key, e.Key.t1 ^ e.Key.t2, lwe.DecryptBit(CipherLWE.Xor(e.Value.Item1, e.Value.Item2))))
                    .ToArray();

                if (l < 10)
                {
                    tableNot.Println($"Test[{l}] NOT");
                    tableXor.Println($"Test[{l}] XOR");
                }

                if (tableNot.Any(e => e.Item2 != e.Item3))
                    throw new($"step[{l}] N:{n} Q:{lwe.Q}");

                if (tableXor.Any(e => e.Item2 != e.Item3))
                    throw new($"step[{l}] N:{n} Q:{lwe.Q}");
            }

            Console.WriteLine("...");
            Console.WriteLine($"SUCCESS ALL {nbTrials} NOT and XOR gates {lwe.Params}");
            Console.WriteLine();

            for (int l = 0; l < nbTrials; ++l)
            {
                var table = 2.Range().Grid2D().ToDictionary(e => e, e => (lwe.EncryptBit(e.t1), lwe.EncryptBit(e.t2)));
                var tableAnd = table
                    .Select(e => (e.Key, e.Key.t1 & e.Key.t2, 
                        lwe.DecryptBit(CipherLWE.And(e.Value.Item1, e.Value.Item2, lwe.EK))))
                    .ToArray();
                if (l < 10)
                    tableAnd.Println($"Test[{l}] AND {(tableAnd.All(e => e.Item2 == e.Item3) ? "SUCCESS" : "FAIL")}");

                if (tableAnd.Any(e => e.Item2 != e.Item3))
                    throw new($"step[{l}] N:{n} Q:{lwe.Q}");
            }

            Console.WriteLine("...");
            Console.WriteLine($"SUCCESS ALL {nbTrials} AND gates {lwe.Params}");
            Console.WriteLine();

            for (int l = 0; l < nbTrials; ++l)
            {
                var table = 2.Range().Grid2D().ToDictionary(e => e, e => (lwe.EncryptBit(e.t1), lwe.EncryptBit(e.t2)));
                var tableAnd = table
                    .Select(e => (e.Key, 1 - (e.Key.t1 & e.Key.t2),
                        lwe.DecryptBit(CipherLWE.Nand(e.Value.Item1, e.Value.Item2, lwe.EK))))
                    .ToArray();
                if (l < 10)
                    tableAnd.Println($"Test[{l}] NAND {(tableAnd.All(e => e.Item2 == e.Item3) ? "SUCCESS" : "FAIL")}");

                if (tableAnd.Any(e => e.Item2 != e.Item3))
                    throw new($"step[{l}] N:{n} Q:{lwe.Q}");
            }

            Console.WriteLine("...");
            Console.WriteLine($"SUCCESS ALL {nbTrials} NAND gates {lwe.Params}");
            Console.WriteLine();
            
            for (int l = 0; l < nbTrials; ++l)
            {
                var table = 2.Range().Grid2D().ToDictionary(e => e, e => (lwe.EncryptBit(e.t1), lwe.EncryptBit(e.t2)));
                var tableAnd = table
                    .Select(e => (e.Key, e.Key.t1 | e.Key.t2,
                        lwe.DecryptBit(CipherLWE.Or(e.Value.Item1, e.Value.Item2, lwe.EK))))
                    .ToArray();
                if (l < 10)
                    tableAnd.Println($"Test[{l}] OR {(tableAnd.All(e => e.Item2 == e.Item3) ? "SUCCESS" : "FAIL")}");

                if (tableAnd.Any(e => e.Item2 != e.Item3))
                    throw new($"step[{l}] N:{n} Q:{lwe.Q}");
            }

            Console.WriteLine("...");
            Console.WriteLine($"SUCCESS ALL {nbTrials} OR gates {lwe.Params}");
            Console.WriteLine();
        }
    }

    public static void TestLWERegev()
    {
        for (int k = 2; k < 15; ++k)
        {
            var reg = new Regev(5 * k);
            reg.Show();

            RunLWERegev(reg, "hello world lwe");
            RunLWERegev(reg, "Hello World LWE");
            RunLWERegev(reg, "AAA+", showCipher: true);

            for (int i = 0; i < 10; i++)
                RunLWERegev(reg, DistributionExt.Dice(LipsumSentences));
        
            // long text
            RunLWERegev(reg, DistributionExt.Dice(LipsumParagraphes), showBinary: false);
            Console.WriteLine();
        }
    }
    // Regev N:10   P:223    M:94     A:0.0287 A*P:6.3903 P/M:2.3723
    // Private key:(189,  85,  81, 135,  17, 149, 115, 149, 115, 151)
    // 
    // hello world lwe
    // seqInput  :[000101101010011000110110001101101111011000000100111011101111011001001110001101100010011000000100001101101110111010100110]
    // seqDecrypt:[000101101010011000110110001101101111011000000100111011101111011001001110001101100010011000000100001101101110111010100110]
    // hello world lwe
    //     SUCCESS
    // ...
    // Regev N:70   P:5261   M:965    A:0.0032 A*P:16.7379 P/M:5.4518
    // Private key:(2969, 3609,  415, 4982, 1963, 2532,  220,  881,  739, 2204, 4253, 2015, 2534, 1946, 4386, 3323,  309,  194, 3308,  372, 5001, 1121, 2467, 4750,  354,  972, 3468,  978, 3398, 4629,  644, 3663, 2027, 2983, 2264,  954, 4139, 1497,  890, 1080,  652, 2244, 2736, 1761,  897, 3653, 1088, 1866, 2816, 1983,  944, 4171, 2519, 5207, 1559, 1273, 2202, 5012,  561, 2827, 1395, 3800, 2941, 1918, 1493, 3268, 4869, 4740, 4755,  774)
    // 
    // hello world lwe
    // seqInput  :[000101101010011000110110001101101111011000000100111011101111011001001110001101100010011000000100001101101110111010100110]
    // seqDecrypt:[000101101010011000110110001101101111011000000100111011101111011001001110001101100010011000000100001101101110111010100110]
    // hello world lwe
    //     SUCCESS
    // 
}
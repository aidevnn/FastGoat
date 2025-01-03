using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Integers;
using static FastGoat.Commons.IntExt;

namespace FastGoat.UserGroup.Lattice;

public static class LWEtests
{
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

    static string RandString(int length)
    {
        var s = "abcdefghijklmnopqrstuvwxyzabcdefghijklmnopqrstuvwxyzabcdefghijklmnopqrstuvwxyz" +
                @"ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789+-*/!@#$%^&*()_+:<>?|\                   ";
        return length.Range().Select(_ => s[Rng.Next(s.Length)]).Glue();
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
        for (int k = 2; k < 11; ++k)
        {
            var lwe = new LWE(k);
            lwe.Show();

            RunLWE(lwe, "hello world lwe");
            RunLWE(lwe, "Hello World LWE");
            RunLWE(lwe, "AAA+", showCipher: true);

            for (int i = 0; i < 10; i++)
                RunLWE(lwe, RandString(Rng.Next(20, 50)));
        
            // long text
            RunLWE(lwe, RandString(1000), showBinary: false);
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
        for(int n = 2; n < 49; ++n) 
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
                RunLWERegev(reg, RandString(Rng.Next(20, 50)));
        
            // long text
            RunLWERegev(reg, RandString(1000), showBinary: false);
            Console.WriteLine();
        }
    }
    // Regev N:10   Q:223    M:156    a:0.0287 aq:6.3903 Q/M:1.4295
    // Private key:(131,  38, 186, 116, 176,   7, 142,  38, 193,  84)
    // 
    // hello world lwe
    // seqInput  :[000101101010011000110110001101101111011000000100111011101111011001001110001101100010011000000100001101101110111010100110]
    // seqDecrypt:[000101101010011000110110001101101111011000000100111011101111011001001110001101100010011000000100001101101110111010100110]
    // hello world lwe
    //     SUCCESS
    // ...
    // Regev N:70   Q:5261   M:1730   a:0.0032 aq:16.7379 Q/M:3.0410
    // Private key:(2775,  254, 2127, 2519, 1431, 4476, 2533, 2966, 4346, 1747, 3786, 3773, 3490, 1261, 4158, 1031, 4407, 1642,  637, 1090,  146, 5114, 2482, 4022, 3940, 5101,  812, 3243, 3980, 1615, 2785, 3247, 1741,   37, 5096, 1367, 1125, 1408, 3695, 1265, 1928, 2761, 3373, 2105,  316, 4651, 2839, 3338, 2325, 2335, 1904, 3176, 4268, 5077,  803, 5247, 4754, 5122, 4505, 1907, 3919, 2791, 2049, 4080, 1707, 4404,  584,  931, 4764, 4154)
    // 
    // hello world lwe
    // seqInput  :[000101101010011000110110001101101111011000000100111011101111011001001110001101100010011000000100001101101110111010100110]
    // seqDecrypt:[000101101010011000110110001101101111011000000100111011101111011001001110001101100010011000000100001101101110111010100110]
    // hello world lwe
    //     SUCCESS
}
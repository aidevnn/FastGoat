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

    public static void TestEncryptDecryptLWE()
    {
        for (int n = 4; n < 31; ++n)
        {
            var lwe = LWE.Setup(n);
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

    public static void TestHELogicGatesLWE()
    {
        for(int n = 4; n < 31; n += 1)
        {
            var lwe = LWE.Setup(n);
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
                    .Select(e => (e, 1 - e, lwe.DecryptBit(LWECipher.Not(lwe.EncryptBit(e))))).ToArray();
                var tableXor = table.Select(e =>
                        (e.Key, e.Key.t1 ^ e.Key.t2, lwe.DecryptBit(LWECipher.Xor(e.Value.Item1, e.Value.Item2))))
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
                        lwe.DecryptBit(LWECipher.And(e.Value.Item1, e.Value.Item2, lwe.EK))))
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
                var tableNand = table
                    .Select(e => (e.Key, 1 - (e.Key.t1 & e.Key.t2),
                        lwe.DecryptBit(LWECipher.Nand(e.Value.Item1, e.Value.Item2, lwe.EK))))
                    .ToArray();
                if (l < 10)
                    tableNand.Println($"Test[{l}] NAND {(tableNand.All(e => e.Item2 == e.Item3) ? "SUCCESS" : "FAIL")}");

                if (tableNand.Any(e => e.Item2 != e.Item3))
                    throw new($"step[{l}] N:{n} Q:{lwe.Q}");
            }

            Console.WriteLine("...");
            Console.WriteLine($"SUCCESS ALL {nbTrials} NAND gates {lwe.Params}");
            Console.WriteLine();
            
            for (int l = 0; l < nbTrials; ++l)
            {
                var table = 2.Range().Grid2D().ToDictionary(e => e, e => (lwe.EncryptBit(e.t1), lwe.EncryptBit(e.t2)));
                var tableOr = table
                    .Select(e => (e.Key, e.Key.t1 | e.Key.t2,
                        lwe.DecryptBit(LWECipher.Or(e.Value.Item1, e.Value.Item2, lwe.EK))))
                    .ToArray();
                if (l < 10)
                    tableOr.Println($"Test[{l}] OR {(tableOr.All(e => e.Item2 == e.Item3) ? "SUCCESS" : "FAIL")}");

                if (tableOr.Any(e => e.Item2 != e.Item3))
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
    
    public static void TestEncryptDecryptRLWE()
    {
        Ring.DisplayPolynomial = MonomDisplay.StarCaret;
        for (int k = 2; k < 8; ++k)
        {
            var rlwe = new RLWE(1 << k);
            rlwe.Show();

            RunRLWE(rlwe, "hello world lwe");
            RunRLWE(rlwe, "Hello World LWE");
            RunRLWE(rlwe, "AAA+");

            for (int i = 0; i < 10; i++)
                RunRLWE(rlwe, DistributionExt.Dice(LipsumSentences));
        
            // long text
            RunRLWE(rlwe, DistributionExt.Dice(LipsumParagraphes), showBinary: false);
            Console.WriteLine();
        }
    }
    
    public static void TestHELogicGatesRLWE()
    {
        Ring.DisplayPolynomial = MonomDisplay.StarCaret;
        for(int k = 2; k < 6; k += 1)
        {
            var n = 1 << k;
            var rlwe = new RLWE(n);
            rlwe.Show();

            var nbTrials = 100;
            for (int l = 0; l < nbTrials; l++)
            {
                var bit = Rng.Next(2);
                var enc = rlwe.EncryptBit(bit);
                var d = rlwe.DecryptBit(enc);
                if (l < 10)
                    Console.WriteLine($"b:{bit} => {d}");

                if (d != bit)
                    throw new($"step[{l}] N:{n} Q:{rlwe.Q}");
            }

            for (int l = 0; l < nbTrials; ++l)
            {
                var tableNot = 2.Range()
                    .Select(e => (e, 1 - e, rlwe.DecryptBit(rlwe.NOT(rlwe.EncryptBit(e))))).ToArray();

                if (l < 10)
                    tableNot.Println($"Test[{l}] NOT");

                if (tableNot.Any(e => e.Item2 != e.Item3))
                    throw new($"step[{l}] N:{n} Q:{rlwe.Q}");
            }

            Console.WriteLine("...");
            Console.WriteLine($"SUCCESS ALL {nbTrials} NOT gates {rlwe.Params}");
            Console.WriteLine();

            for (int l = 0; l < nbTrials; ++l)
            {
                var table = 2.Range().Grid2D().ToDictionary(e => e, e => (rlwe.EncryptBit(e.t1), rlwe.EncryptBit(e.t2)));
                var tableAnd = table
                    .Select(e => (e.Key, e.Key.t1 & e.Key.t2, rlwe.DecryptBit(rlwe.AND(e.Value.Item1, e.Value.Item2))))
                    .ToArray();
                if (l < 10)
                    tableAnd.Println($"Test[{l}] AND {(tableAnd.All(e => e.Item2 == e.Item3) ? "SUCCESS" : "FAIL")}");

                if (tableAnd.Any(e => e.Item2 != e.Item3))
                    throw new($"step[{l}] N:{n} Q:{rlwe.Q}");
            }

            Console.WriteLine("...");
            Console.WriteLine($"SUCCESS ALL {nbTrials} AND gates {rlwe.Params}");
            Console.WriteLine();

            for (int l = 0; l < nbTrials; ++l)
            {
                var table = 2.Range().Grid2D().ToDictionary(e => e, e => (rlwe.EncryptBit(e.t1), rlwe.EncryptBit(e.t2)));
                var tableNand = table
                    .Select(e => (e.Key, 1 - (e.Key.t1 & e.Key.t2),
                        rlwe.DecryptBit(rlwe.NAND(e.Value.Item1, e.Value.Item2))))
                    .ToArray();
                if (l < 10)
                    tableNand.Println($"Test[{l}] NAND {(tableNand.All(e => e.Item2 == e.Item3) ? "SUCCESS" : "FAIL")}");

                if (tableNand.Any(e => e.Item2 != e.Item3))
                    throw new($"step[{l}] N:{n} Q:{rlwe.Q}");
            }

            Console.WriteLine("...");
            Console.WriteLine($"SUCCESS ALL {nbTrials} NAND gates {rlwe.Params}");
            Console.WriteLine();
            
            for (int l = 0; l < nbTrials; ++l)
            {
                var table = 2.Range().Grid2D().ToDictionary(e => e, e => (rlwe.EncryptBit(e.t1), rlwe.EncryptBit(e.t2)));
                var tableOr = table
                    .Select(e => (e.Key, e.Key.t1 | e.Key.t2, rlwe.DecryptBit(rlwe.OR(e.Value.Item1, e.Value.Item2))))
                    .ToArray();
                if (l < 10)
                    tableOr.Println($"Test[{l}] OR {(tableOr.All(e => e.Item2 == e.Item3) ? "SUCCESS" : "FAIL")}");

                if (tableOr.Any(e => e.Item2 != e.Item3))
                    throw new($"step[{l}] N:{n} Q:{rlwe.Q}");
            }

            Console.WriteLine("...");
            Console.WriteLine($"SUCCESS ALL {nbTrials} OR gates {rlwe.Params}");
            Console.WriteLine();
            
            for (int l = 0; l < nbTrials; ++l)
            {
                var table = 2.Range().Grid2D().ToDictionary(e => e, e => (rlwe.EncryptBit(e.t1), rlwe.EncryptBit(e.t2)));
                var tableXor = table
                    .Select(e => (e.Key, e.Key.t1 ^ e.Key.t2, rlwe.DecryptBit(rlwe.XOR(e.Value.Item1, e.Value.Item2))))
                    .ToArray();
                if (l < 10)
                    tableXor.Println($"Test[{l}] OR {(tableXor.All(e => e.Item2 == e.Item3) ? "SUCCESS" : "FAIL")}");

                if (tableXor.Any(e => e.Item2 != e.Item3))
                    throw new($"step[{l}] N:{n} Q:{rlwe.Q}");
            }

            Console.WriteLine("...");
            Console.WriteLine($"SUCCESS ALL {nbTrials} XOR gates {rlwe.Params}");
            Console.WriteLine();
        }
    }

    #region Boolean arithmetic
    
    static int OpBool(string name, int c1, int c2)
    {
        return name switch
        {
            "and" => c1 & c2,
            "nand" => 1 - (c1 & c2),
            "nor" => 1 - (c1 | c2),
            "or" => c1 | c2,
            _ => c1 ^ c2
        };
    }

    static long Bin2Int(int[] l) => l.Reverse().Aggregate((long)0, (acc, i) => acc * 2 + i);

    static void add(int[] f1, int[] f2)
    {
        var l = f1.Length;
        var a0 = 0;
        for (int i = 0; i < l; i++)
        {
            if (i == 1) a0 = f2[0];
            if (i > 1) (f2[i - 1], a0) = (a0, f2[i - 1]);
            var (e1, e2) = (f1[i], f2[i]);
            (f1[i], f2[i]) = (e1 ^ e2, e1 & e2);
        }
    
        f2[l - 1] = a0;
        f2[0] = 0;
    }
    #endregion
    
    static void RunAddWithCarryRLWE(int N, int bits, int nbTrials = 50)
    {
        var rlwe = new RLWE(N);
        rlwe.Show();
        var mod = (long)1 << bits;

        GlobalStopWatch.Restart();
        for (int k = 0; k < nbTrials; ++k)
        {
            var m1 = DistributionExt.DiceSample(bits - 1, [0, 1]).Append(0).ToArray();
            var e1 = rlwe.Encrypt(m1);
            var m2 = DistributionExt.DiceSample(bits - 1, [0, 1]).Append(0).ToArray();
            var e2 = rlwe.Encrypt(m2);

            var (f1, f2) = (m1.ToArray(), m2.ToArray());
            for (int i = 0; i < bits; i++)
                add(f1, f2);

            var (a1, a2) = (Bin2Int(m1), Bin2Int(m2));
            var sumi = (a1 + a2) % mod;
            var sumf = Bin2Int(f1) % mod;

            var add_c1c2 = rlwe.ADD(e1, e2);
            var dc = rlwe.Decrypt(add_c1c2);

            var fmt = $"{{0,{(int)double.Log10(mod) + 1}}}";
            string FMT(long a) => string.Format(fmt, a);

            Console.WriteLine($"Addition No{k + 1} {rlwe.Params}");
            Console.WriteLine($"   0b{m1.Reverse().Glue()} = {FMT(a1)}");
            Console.WriteLine($" + 0b{m2.Reverse().Glue()} = {FMT(a2)}");
            Console.WriteLine($" = 0b{dc.Reverse().Glue()} = {FMT(sumi)}");
            Console.WriteLine();
            if (sumi != sumf || !dc.SequenceEqual(f1))
                throw new();
        }
        
        GlobalStopWatch.Show($"END {nbTrials} tests of {bits}-bits Addition with carry");
        Console.WriteLine();
    }

    public static void TestAddWithCarryRLWE()
    {
        Ring.DisplayPolynomial = MonomDisplay.StarCaret;
        var (bits, nbTrials) = (32, 10);
        RunAddWithCarryRLWE(N: 4, bits, nbTrials);
        RunAddWithCarryRLWE(N: 8, bits, nbTrials);
        RunAddWithCarryRLWE(N: 16, bits, nbTrials);
        RunAddWithCarryRLWE(N: 32, bits, nbTrials);
        
        RunAddWithCarryRLWE(N: 128, bits, nbTrials);
        // # END 10 tests of 32-bits Addition with carry Time:3m56s
    }
    // Addition No9 RLWE N=32=2^5, Φ(N)=16 t=193 q=38021
    //    0b01011001101001011000100111111110 = 1504020990
    //  + 0b00111110101000000000000010010011 = 1050673299
    //  = 0b10011000010001011000101010010001 = 2554694289
    // 
    // Addition No10 RLWE N=32=2^5, Φ(N)=16 t=193 q=38021
    //    0b00110110010000101010011010011110 =  910337694
    //  + 0b01011001100010011100100101110100 = 1502202228
    //  = 0b10001111110011000111000000010010 = 2412539922
    // 
    // # END 10 tests of 32-bits Addition with carry Time:16.973s
    // 
}
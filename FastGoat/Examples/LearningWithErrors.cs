using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Integers;
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

    static void SymbBGV()
    {
        var ind = new Indeterminates<Xi>(MonomOrder.GrLex, new Xi("k"), new Xi("T"), new Xi("X"), new Xi("sk"));
        ind.ExtendAppend("epk", "pkb", "erlk", "rlkb");
        ind.ExtendAppend("m1", "u1", "e1a", "e1b");
        ind.ExtendAppend("m2", "u2", "e2a", "e2b");
        ind.ExtendAppend("b");

        var z = new Polynomial<Rational, Xi>(ind, Rational.KZero());
        var Vars = z.Variables;
        var (k, t, X, sk) = Vars.Take(4).Deconstruct();
        var (epk, pkb, erlk, rlkb) = Vars.Skip(4).Take(4).Deconstruct();
        var pka = (t * epk + pkb * sk);
        var b = Vars.Last();
        var p = b * t + 1;
        var rlka = (t * erlk + rlkb * sk + p * sk.Pow(2));

        Console.WriteLine("##### Encryption/Decryption #####");
        Console.WriteLine("Decrypt pk: c0pk - sk * c1pk       = {0}", pka - sk * pkb);
        Console.WriteLine("Decrypt pk: c0pk - sk * c1pk mod t = {0}", (pka - sk * pkb).Div(t).rem);
        Console.WriteLine();

        Console.WriteLine("Decrypt rlk: rlka - sk * rlkb       = {0}", rlka - sk * rlkb);
        Console.WriteLine("Decrypt rlk: rlka - sk * rlkb mod t = {0}", (rlka - sk * rlkb).Div(t).rem);
        Console.WriteLine("                             sk^2   = {0}", sk.Pow(2));
        Console.WriteLine();

        var (m1, u1, e1a, e1b) = Vars.Skip(8).Take(4).Deconstruct();
        var m1a = m1 + t * e1a + u1 * pka;
        var m1b = t * e1b + u1 * pkb;
        Console.WriteLine($"m1a = {m1a}");
        Console.WriteLine($"m1b = {m1b}");
        Console.WriteLine($"Decrypt m1: m1a - sk*m1b       = {m1a - sk * m1b}");
        Console.WriteLine($"Decrypt m1: m1a - sk*m1b mod t = {(m1a - sk * m1b).Div(t).rem}");
        Console.WriteLine();

        var (m2, u2, e2a, e2b) = Vars.Skip(12).Take(4).Deconstruct();
        var m2a = m2 + t * e2a + u2 * pka;
        var m2b = t * e2b + u2 * pkb;
        Console.WriteLine($"m2a = {m2a}");
        Console.WriteLine($"m2b = {m2b}");
        Console.WriteLine($"Decrypt m2: m2a - sk*m1b       = {m2a - sk * m2b}");
        Console.WriteLine($"Decrypt m2: m2a - sk*m2b mod t = {(m2a - sk * m2b).Div(t).rem}");
        Console.WriteLine();

        var add_m1m2 = m1 + m2;
        var add_m1m2a = m1a + m2a;
        var add_m1m2b = m1b + m2b;
        Console.WriteLine("##### Homomorphic Addition of Cipher #####");
        Console.WriteLine($"add_m1m2a = {add_m1m2a}");
        Console.WriteLine($"add_m1m2b = {add_m1m2b}");
        Console.WriteLine($"Decrypt add_m1m2: add_m1m2a - sk*add_m1m2b       = {add_m1m2a - sk * add_m1m2b}");
        Console.WriteLine(
            $"Decrypt add_m1m2: add_m1m2a - sk*add_m1m2b mod t = {(add_m1m2a - sk * add_m1m2b).Div(t).rem}");
        Console.WriteLine($"                                  add_m1m2       = {add_m1m2}");
        Console.WriteLine();

        var mul_km1 = k * m1;
        var mul_km1a = k * m1a;
        var mul_km1b = k * m1b;
        Console.WriteLine("##### Homomorphic Multiplication by Scalar #####");
        Console.WriteLine($"mul_km1a = {mul_km1a}");
        Console.WriteLine($"mul_km1b = {mul_km1b}");
        Console.WriteLine($"Decrypt mul_km1: mul_km1a - sk*mul_km1b       = {mul_km1a - sk * mul_km1b}");
        Console.WriteLine($"Decrypt mul_km1: mul_km1a - sk*mul_km1b mod t = {(mul_km1a - sk * mul_km1b).Div(t).rem}");
        Console.WriteLine($"                                mul_km1       = {mul_km1}");
        Console.WriteLine();

        var mul_m1m2 = m1 * m2;

        var c0 = m1a * m2a;
        var c1 = m1a * m2b + m1b * m2a;
        var c2 = m1b * m2b;

        var mul_m1m2a = p * c0 + c2 * rlka;
        var mul_m1m2b = p * c1 + c2 * rlkb;

        Console.WriteLine("##### Homomorphic Multiplication of Ciphers #####");
        Console.WriteLine($"mul_m1m2a = {mul_m1m2a}");
        Console.WriteLine($"mul_m1m2a mod p = {mul_m1m2a.Div(p).rem}");
        Console.WriteLine($"mul_m1m2b = {mul_m1m2b}");
        Console.WriteLine();
        Console.WriteLine($"Decrypt mul_m1m2: mul_m1m2a - sk*mul_m1m2b       = {mul_m1m2a - sk * mul_m1m2b}");
        Console.WriteLine();
        Console.WriteLine(
            $"Decrypt mul_m1m2: mul_m1m2a - sk*mul_m1m2b mod t = {(mul_m1m2a - sk * mul_m1m2b).Div(t).rem}");
        Console.WriteLine($"                                      mul_m1m2   = {mul_m1m2}");
        Console.WriteLine();

        var ska = sk + t * e1a + u1 * pka;
        var skb = t * e1b + u1 * pkb;
        var swk_m2a = m2a - m2b * ska;
        var swk_m2b = -m2b * skb;
        Console.WriteLine("##### Switch key #####");
        Console.WriteLine($"swk_m2a = {swk_m2a}");
        Console.WriteLine($"swk_m2b = {swk_m2b}");
        Console.WriteLine($"Decrypt swk_m2: swk_m2a - sk*swk_m2b       = {swk_m2a - sk * swk_m2b}");
        Console.WriteLine($"Decrypt swk_m2: swk_m2a - sk*swk_m2b mod t = {(swk_m2a - sk * swk_m2b).Div(t).rem}");
        Console.WriteLine($"                                     m2    = {m2}");
        Console.WriteLine();
    }

    public static void RunLeveledBGV(int N, int t0, int level, bool differentPrimes = true)
    {
        var sk = RLWE.SKBGV(N / 2);
        var (pm, _, t, primes, pk, swks) = RLWE.SetupBGV(N, t0, level, sk, differentPrimes);
        Console.WriteLine($"pm = {pm} T = {t} Primes = [{primes.Glue(", ")}]");
        Console.WriteLine($"sk = {sk}");
        Console.WriteLine($"pk => {pk.Params}");

        foreach (var swk in swks)
            Console.WriteLine($"swk[{swk.Key}] => {swk.Value.skPow[0].Params}");

        Console.WriteLine();

        var size = 1 << level;
        var n = pm.Degree;

        var seqMsg = size.SeqLazy().Select(_ => RLWE.GenUnif(n, t)).ToArray();
        var seqCipher = seqMsg.Select(m => RLWE.EncryptBGV(m, pk)).ToArray();
        var mul = seqMsg.Aggregate((xi, xj) => (xi * xj).ResModSigned(pm, t));

        var qMul = new Queue<RLWECipher>(seqCipher);
        while (qMul.Count > 1)
        {
            var c1 = qMul.Dequeue();
            var c2 = qMul.Dequeue();
            var (nextMod, swk) = swks[c1.Q];
            var c1c2 = RLWE.MulRelinBGV(c1, c2, swk[2]).ModSwitch(nextMod);
            qMul.Enqueue(c1c2);
        }

        var cMul = qMul.Dequeue();
        var d_mul = RLWE.DecryptBGV(cMul, sk);

        if (size < 17)
            seqMsg.Println($"level:{level} Size:{size}");
        else
        {
            seqMsg.Take(8).Println($"level:{level} Size:{size}");
            seqMsg.TakeLast(8).Println("    ...");
        }

        Console.WriteLine(" *  {0}", Enumerable.Repeat('-', seqMsg.Max(l => $"{l}".Length)).Glue());
        Console.WriteLine($" =  {mul}");
        Console.WriteLine($"    {d_mul}");
        if (!d_mul.Equals(mul))
            throw new("fail");
    }

    public static void RunLeveledBGV(int N, int t0, int level, bool differentPrimes = true, int nbTests = 500)
    {
        GlobalStopWatch.AddLap();
        for (int i = 0; i < nbTests; i++)
        {
            Console.WriteLine($"Test[{i + 1}]");
            RunLeveledBGV(N, t0, level, differentPrimes);
            Console.WriteLine();
        }

        GlobalStopWatch.Show($"Pass {nbTests} tests. N={N} T={t0} Level={level}.");
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
        // Weak parameters
        // RLWE N=16=2^4, Φ(N)=8 PM=x^8 + 1 t=17 q=1361
        var rlwe = new RLWE(16);
        var (n, pm, sk, t, q, pk, rlk) = rlwe;
        rlwe.Show();

        var m1 = RLWE.GenUnif(n, t);
        var e1 = RLWE.EncryptBGV(m1, pk);
        var m2 = RLWE.GenUnif(n, t);
        var e2 = RLWE.EncryptBGV(m2, pk);
        var k = RLWE.GenUnif(n, t);

        e1.Show($"e1 = Encrypt(m1 = {m1})");
        e2.Show($"e2 = Encrypt(m2 = {m2})");
        var d1 = RLWE.DecryptBGV(e1, sk);
        var d2 = RLWE.DecryptBGV(e2, sk);
        Console.WriteLine($"m1          = {e1}");
        Console.WriteLine($"Decrypt(e1) = {d1}");
        Console.WriteLine($"m2          = {e2}");
        Console.WriteLine($"Decrypt(e2) = {d2}");
        Console.WriteLine();
        if (!d1.Equals(m1) || !d2.Equals(m2))
            throw new("encrypt decrypt");

        var add_m1m2 = (m1 + m2).CoefsModSigned(t);
        var add_e1e2 = e1 + e2;
        var d_add = RLWE.DecryptBGV(add_e1e2, sk);

        add_e1e2.Show("e1 + e2");
        Console.WriteLine($"m1 + m2          = {add_m1m2}");
        Console.WriteLine($"Decrypt(e1 + e2) = {d_add}");
        Console.WriteLine();
        if (!d_add.Equals(add_m1m2))
            throw new("m1 + m2");

        var km1 = (k * m1).ResModSigned(pm, t);
        var ke1 = k * e1;
        var d_k1 = RLWE.DecryptBGV(ke1, sk);

        ke1.Show("k * e1");
        Console.WriteLine($"k               = {k}");
        Console.WriteLine($"k * m1          = {km1}");
        Console.WriteLine($"Decrypt(k * e1) = {d_k1}");
        Console.WriteLine();
        if (!d_k1.Equals(km1))
            throw new("k * m1");

        var mul_m1m2 = (m1 * m2).ResModSigned(pm, t);
        var mul_e1e2 = RLWE.MulRelinBGV(e1, e2, rlk).ModSwitch(q);
        var d_mul = RLWE.DecryptBGV(mul_e1e2, sk);

        mul_e1e2.Show("e1 * e2");
        Console.WriteLine($"m1 * m2          = {mul_m1m2}");
        Console.WriteLine($"Decrypt(e1 * e2) = {d_mul}");
        Console.WriteLine();
        if (!d_mul.Equals(mul_m1m2))
            throw new("m1 * m2");
    }

    public static void Example4LogicGates()
    {
        // Weak parameters
        // RLWE N=16=2^4, Φ(N)=8 PM=x^8 + 1 t=17 q=1361
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

    #region Fail

    public static void Example5HomomorphicAdditionWithCarry()
    {
        // Weak parameters
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
        // Weak parameters
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
        // Weak parameters
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

    #endregion

    public static void Example9WrongParameters()
    {
        for (int k = 2; k < 8; k++)
        {
            var n = 1 << (k - 1);
            var t0 = IntExt.Primes10000.First(t0 => t0 % (2 * n) == 1);
            var _q0 = IntExt.Primes10000.First(t1 => t1 % (2 * n) == 1 && t1 > t0);
            var q0 = new Rational(_q0) * t0; // ciphertext and plaintext not coprimes
            var pm = FG.QPoly().Pow(n) + 1;
            var t = new Rational(t0);
            var sk = 10000.SeqLazy().Select(_ => RLWE.GenTernary(n))
                .First(s => !s[n - 1].IsZero() && s.Coefs.Count(e => e.IsZero()) <= n / 4);

            var epk = RLWE.GenDiscrGauss(n);
            var c1pk = RLWE.GenUnif(n, q0);
            var c0pk = (t * epk + c1pk * sk).ResModSigned(pm, q0);
            var pk = new RLWECipher(c0pk, c1pk, pm, t, q0);

            var pka = pk.A.ToZnPoly(t0);
            var pkb = pk.B.ToZnPoly(t0);

            var a = FG.EPoly(pka.X.Pow(n) + 1, 'a');
            Console.WriteLine(new { n, t, q0, pm });
            var sk1 = sk.ToZnPoly(t0).Substitute(a);
            var sk2 = pka.Substitute(a) / pkb.Substitute(a);
            Console.WriteLine($"sk = {sk}");
            Console.WriteLine($"   = {sk1}");
            Console.WriteLine($"   = {sk2}");
            Console.WriteLine($"Invalid:{sk1.Equals(sk2)}");
            Console.WriteLine();
        }
    }

    public static void Example10LeveledBGV()
    {
        IntExt.RecomputeAllPrimesUpTo(1500000);

        RunLeveledBGV(16, 769, 6, true, 5);
    }
}
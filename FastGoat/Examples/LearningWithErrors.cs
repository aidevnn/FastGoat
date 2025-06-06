using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.LWE;
using Rq = FastGoat.Structures.VecSpace.KPoly<FastGoat.UserGroup.Integers.Rational>;

namespace FastGoat.Examples;

public static class LearningWithErrors
{
    static LearningWithErrors()
    {
        GlobalStopWatch.Restart();
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
            throw new("FAIL");

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
            throw new("FAIL");

        Console.WriteLine();
    }

    static void RunLeveledBGV(int N, int level)
    {
        var (pm, sk, t, primes, sp, pk, rlks) = RLWE.SetupBGV(N, level);
        Console.WriteLine(
            $"N = {N} pm = {pm} T = {t} Primes = [{primes.Glue(", ")}] sp = {sp} Sigma = {RLWE.Sigma(pm.Degree, t):f3}");
        Console.WriteLine($"pk => {pk.Params}");

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
            var (nextMod, rlk) = rlks[c1.Q];
            var c1c2 = RLWE.MulRelinBGV(c1, c2, rlk).ModSwitch(nextMod);
            qMul.Enqueue(c1c2);
        }

        var cMul = qMul.Dequeue();
        var d_mul = RLWE.DecryptBGV(cMul, sk);

        if (n < 17)
        {
            Console.WriteLine($"sk = {sk}");
            foreach (var rlk in rlks)
                Console.WriteLine($"rlk[{rlk.Key}] => {rlk.Value.rlk.Params}");

            Console.WriteLine();

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
        }

        if (!d_mul.Equals(mul))
            throw new("fail");
    }

    static void RunLeveledBGV(int N, int level, int nbTests)
    {
        GlobalStopWatch.AddLap();
        for (int i = 0; i < nbTests; i++)
        {
            Console.WriteLine($"Test[{i + 1}]");
            RunLeveledBGV(N, level);
            Console.WriteLine();
        }

        GlobalStopWatch.Show($"Pass {nbTests} tests. N={N} Level={level}.");
        Console.WriteLine();
    }

    static void CheckBR(Rational[] s, Rq pm, (Rational ai, Rational[] bi) ab, Rq f, Rq actual, Rational q)
    {
        var n = pm.Degree;
        var (ai, bi) = ab;
        Console.WriteLine($"f           :{f}");
        Console.WriteLine($"f           :[{f.CoefsExtended(n - 1).Glue(", ", "{0,4}")}]");
        Console.WriteLine($"ai          :{ai}");
        Console.WriteLine($"bi          :[{bi.Glue(", ", "{0,4}")}]");
        Console.WriteLine($"blind rotate:{actual}");
        Console.WriteLine($"            :[{actual.CoefsExtended(n - 1).Glue(", ", "{0,4}")}]");

        // Testing result
        var u = (ai - bi.Zip(s).Aggregate(q.Zero, (sum, e) => sum + (e.First * e.Second)));
        var expected = (RLWE.XpowA((int)u.Num, pm, q) * f).ResModSigned(pm, q);
        var check = (expected - actual).IsZero();
        Console.WriteLine($"u= a - <b,s>:{u}");
        Console.WriteLine($"f*X^{u,-4}    :{expected}");
        Console.WriteLine($"f*X^{u,-4}    :[{expected.CoefsExtended(n - 1).Glue(", ", "{0,4}")}]");
        Console.WriteLine($"{(check ? "success" : "fail")}");
        Console.WriteLine();
        if (!check)
            throw new();
    }

    public static void Example1Regev()
    {
        GlobalStopWatch.AddLap();
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
        GlobalStopWatch.AddLap();
        // Weak parameters
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

    public static void Example3AdditionMultiplicationBGV()
    {
        // Weak parameters
        var rlwe = new RLWE(32);
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

    public static void Example3SymbBGV()
    {
        var ind = new Indeterminates<Xi>(MonomOrder.GrLex, new Xi("k"), new Xi("T"), new Xi("X"), new Xi("sk"),
            new Xi("sk2"));
        ind.ExtendAppend("epk", "pkb", "erlk", "rlkb");
        ind.ExtendAppend("m1", "u1", "e1a", "e1b");
        ind.ExtendAppend("m2", "u2", "e2a", "e2b");
        ind.ExtendAppend("b");

        var z = new Polynomial<Rational, Xi>(ind, Rational.KZero());
        var Vars = z.AllVariables;
        var (k, t, X, sk, sk2) = Vars.Take(5).Deconstruct();
        var (epk, pkb, erlk, rlkb) = Vars.Skip(5).Take(4).Deconstruct();
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

        var (m1, u1, e1a, e1b) = Vars.Skip(9).Take(4).Deconstruct();
        var m1a = m1 + t * e1a + u1 * pka;
        var m1b = t * e1b + u1 * pkb;
        Console.WriteLine($"m1a = {m1a}");
        Console.WriteLine($"m1b = {m1b}");
        Console.WriteLine($"Decrypt m1: m1a - sk*m1b       = {m1a - sk * m1b}");
        Console.WriteLine($"Decrypt m1: m1a - sk*m1b mod t = {(m1a - sk * m1b).Div(t).rem}");
        Console.WriteLine();

        var (m2, u2, e2a, e2b) = Vars.Skip(13).Take(4).Deconstruct();
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

        var ska = sk2 + t * e1a + u1 * pka;
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

    public static void Example4LogicGates()
    {
        // Weak parameters
        var rlwe = new RLWE(32);
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

        foreach (var opName in new[] { "AND", "NAND", "NOR", "OR", "XOR" })
        {
            Console.WriteLine(opName);
            Console.WriteLine($"   [{m1.Glue()}]");
            Console.WriteLine($"   [{m2.Glue()}]");
            Console.WriteLine($" = [{rlwe.Decrypt(rlwe.OP(opName, e1, e2)).Glue()}]");
            Console.WriteLine();
        }
    }

    public static void Example5Regev2RLWE()
    {
        // Weak parameters
        var (reg, rlwe, swk, exsk) = Regev.SetupRLWE(16);
        reg.Show();
        rlwe.Show();

        for (int k = 0; k < 10; ++k)
        {
            var m = DistributionExt.DiceSample(reg.N, [0, 1]).ToArray();

            var c11 = reg.Encrypt(m);
            var c21 = rlwe.FromRegevCipher(c11, exsk);
            var d1 = rlwe.Decrypt(c21);

            var c12 = rlwe.Encrypt(m);
            var c22 = rlwe.ToRegevCipher(c12, swk);
            var d2 = reg.Decrypt(c22);

            Console.WriteLine($"m :[0b{m.Reverse().Glue()}]");
            Console.WriteLine($"d1:[0b{d1.Reverse().Glue()}] Regev to RLWE");
            Console.WriteLine($"d2:[0b{d2.Reverse().Glue()}] RLWE  to Regev");
            Console.WriteLine();

            if (!m.SequenceEqual(d1) || !m.SequenceEqual(d2))
                throw new($"step[{k}]");
        }
    }

    public static void Example6WrongParameters()
    {
        for (int k = 2; k < 8; k++)
        {
            var n = 1 << (k - 1);
            var t = RLWE.CiphertextModulusBGV(n);
            var q = RLWE.FirstPrimeEqualOneMod(2 * n * t, t) * t;
            var pm = FG.QPoly().Pow(n) + 1;
            Console.WriteLine(new { n, t, q, pm });
            var t0 = (int)t.Num;
            var sk = RLWE.SKBGV(n);

            var epk = RLWE.GenDiscrGauss(n, RLWE.Sigma(n, t));
            var c1pk = RLWE.GenUnif(n, q);
            var c0pk = (t * epk + c1pk * sk).ResModSigned(pm, q);
            var pk = new RLWECipher(c0pk, c1pk, pm, t, q);

            var pka = pk.A.ToZnPoly(t0);
            var pkb = pk.B.ToZnPoly(t0);

            var a = FG.EPoly(pka.X.Pow(n) + 1, 'a');
            var sk1 = sk.ToZnPoly(t0).Substitute(a);
            var sk2 = pka.Substitute(a) / pkb.Substitute(a);
            Console.WriteLine($"sk = {sk}");
            Console.WriteLine($"   = {sk1}");
            Console.WriteLine($"   = {sk2}");
            Console.WriteLine($"Invalid:{sk1.Equals(sk2)}");
            Console.WriteLine();
        }
    }

    public static void Example7LeveledBGV()
    {
        // Weak parameters
        IntExt.RecomputeAllPrimesUpTo(500000);

        var nbTests = 5;
        RunLeveledBGV(N: 16, level: 8, nbTests);
        RunLeveledBGV(N: 32, level: 8, nbTests);
        RunLeveledBGV(N: 64, level: 6, nbTests);
        RunLeveledBGV(N: 128, level: 4, nbTests);
        RunLeveledBGV(N: 256, level: 4, nbTests);
        RunLeveledBGV(N: 512, level: 2, nbTests);
        RunLeveledBGV(N: 1024, level: 2, nbTests);
        RunLeveledBGV(N: 2048, level: 1, nbTests);
    }

    public static void Example8RGSW()
    {
        // Weak parameters

        var N = 32;
        var n = N / 2;
        var level = 2;
        var (pm, sk, t, primes, sp, pk, rlks) = RLWE.SetupBGV(N, level);
        var q = primes[0];
        var qL = pk.Q;
        var B = RLWE.RNSGadgetBase(primes);
        Console.WriteLine($"BGV level = {level}, Primes = [{primes.Glue(", ")}]");
        Console.WriteLine($"pm = {pm} T = {t} q = {q} sp = {sp} qL = {qL} Sigma = {RLWE.Sigma(n, t):f4}");
        Console.WriteLine($"sk = {sk}");
        Console.WriteLine($"pk => {pk.Params}");
        Console.WriteLine();

        for (int k = 0; k < 5; ++k)
        {
            var m1 = RLWE.GenUnif(n, t);
            var cm1 = RLWE.EncryptBGV(m1, pk);
            var m2 = RLWE.GenUnif(n, t);
            var _cm2 = RLWE.EncryptBGV(m2, pk);
            var (csm2, cm2) = RLWE.EncryptRgswBGV(m2, pk, B);
            var m1m2 = (m1 * m2).ResModSigned(pm, t);
            var cm1m2gsw = RLWE.MulRgsw(RLWE.DecompRNS(cm1, primes), cm2, csm2);
            var cm1m2rlk = RLWE.MulRelinBGV(cm1, _cm2, rlks[qL].rlk).ModSwitch(qL);
            var dm1m2gsw = RLWE.DecryptBGV(cm1m2gsw, sk);
            var dm1m2rlk = RLWE.DecryptBGV(cm1m2rlk, sk);
            Console.WriteLine($"m1      = {m1}");
            Console.WriteLine($"m2      = {m2}");
            Console.WriteLine($"m1 * m2 = {m1m2}");
            Console.WriteLine($"        = {dm1m2gsw}");
            Console.WriteLine($"        = {dm1m2rlk}");
            Console.WriteLine($"egsw    = {RLWE.ErrorsBGV(cm1m2gsw, sk).NormInf()}");
            Console.WriteLine();
            if (!dm1m2gsw.Equals(m1m2))
                throw new();
        }
    }

    public static void Example9BlindRotate()
    {
        // Weak parameters

        var N = 32;
        var n = N / 2;
        var level = 2;
        var (pm, sk, t, primes, sp, pk, _) = RLWE.SetupBGV(N, level);
        var q = primes[0];
        var qL = pk.Q;
        var B = RLWE.RNSGadgetBase(primes);
        Console.WriteLine($"BGV level = {level}, Primes = [{primes.Glue(", ")}]");
        Console.WriteLine($"pm = {pm} T = {t} q = {q} sp = {sp} qL = {qL} sigma = {RLWE.Sigma(n, t):f4}");
        Console.WriteLine($"sk = {sk}");
        Console.WriteLine($"pk => {pk.Params}");
        Console.WriteLine();

        var brk = RLWE.BRKgswBGV(sk, pk, B);

        var x = pm.X;
        var c = n / 2 - 1;
        var f = (2 * c + 1).SeqLazy(-c).Select(j => j * RLWE.XpowA(j, pm, t))
            .Aggregate(x.Zero, (acc, v) => acc + v).ResModSigned(pm, t);

        var s = n.SeqLazy().Select(i => sk[i].Signed(t)).ToVec();

        for (int k = 0; k < 50; ++k)
        {
            GlobalStopWatch.AddLap();
            var ai = new Rational(IntExt.Rng.Next(1, n + 1)).Signed(n);
            var bi = RLWE.GenUnif(n, n).CoefsExtended(n - 1);

            var acc = RLWE.BlindRotategswBGV((ai, bi), brk, B, primes);
            var actual = RLWE.DecryptBGV(acc, sk);
            CheckBR(s.ToArray(), pm, (ai, bi), f, actual, n * t.One);
            GlobalStopWatch.Show();
        }
    }

    public static void Example10Repack()
    {
        var (N, level) = (16, 1);
        var (pm, sk, t, primes, sp, pk, rlks) = RLWE.SetupBGV(N, level);
        var n = N / 2;
        var q = primes[0];
        var qL = pk.Q;

        var ak = RLWE.AKBGV(sk, pk, sp.Pow(level));
        var ni = (1 - qL) / n;
        var x = pm.X;
        var c = pm.Degree / 2 - 1;
        var f = (2 * c + 1).SeqLazy(-c).Select(j => j * RLWE.XpowA(j, pm, qL))
            .Aggregate(x.Zero, (acc, v) => acc + v).ResModSigned(pm, qL);

        Console.WriteLine($"BGV level = {level}");
        Console.WriteLine($"pm = {pm} T = {t} q = {q} sp = {sp} qL = {qL}");
        Console.WriteLine($"sk = {sk}");
        Console.WriteLine($"pk => {pk.Params}");
        Console.WriteLine();

        for (int k = 0; k < 10; ++k)
        {
            var uis = n.SeqLazy().Select(_ => IntExt.Rng.Next(-n / 2 + 1, n / 2)).ToArray();
            var u = uis.ToKPoly(Rational.KOne());
            var f_uis = uis.Select(ui => (f * RLWE.XpowA(ui, pm, qL)).ResModSigned(pm, qL)).ToArray();
            var accs = f_uis.Select(fui => new RLWECipher(fui, pm.Zero, pm, t, qL)).ToArray();
            accs.Println("f * X^ui");
            var repack = RLWE.RepackingBGV(accs.Select(e => e * ni).ToArray(), ak);
            Console.WriteLine($"u = {u}");
            repack.Show("repack");
            Console.WriteLine();

            if (!u.Equals(-repack.A))
                throw new();
        }
    }

    public static void Example11Bootstrapping()
    {
        var (N, level) = (16, 4);
        var (pm, sk, t, primes, sp, pk, rlks) = RLWE.SetupBGV(N, level);
        var n = N / 2;
        var q = primes[0];
        var qL = pk.Q;

        var B = RLWE.RNSGadgetBase(primes);
        var tables = RLWE.PrepareNTT(n, t, primes);
        var rlk = rlks[qL].rlk;
        var brk = RLWE.BRKgswBGV(sk, pk, B);
        var brkntt = RLWE.BRKgswNTTBGV(sk, pk, B, tables);
        var ak = RLWE.AKBGV(sk, pk, sp.Pow(level));

        Console.WriteLine($"BGV level = {level}, Primes = [{primes.Glue(", ")}]");
        Console.WriteLine($"pm = {pm} T = {t} q = {q} sp = {sp} qL = {qL}");
        Console.WriteLine($"sk = {sk}");
        Console.WriteLine($"pk => {pk.Params}");
        Console.WriteLine();

        GlobalStopWatch.AddLap();
        for (int k = 0; k < 20; ++k)
        {
            GlobalStopWatch.AddLap();
            var m = RLWE.GenUnif(n, t);
            var m2 = (m * m).ResModSigned(pm, t);
            var cm = RLWE.EncryptBGV(m, pk); // level qL
            var ct = RLWE.MulRelinBGV(cm, cm, rlk).ModSwitch(q); // level q0

            // var ctboot = RLWE.Bootstrapping(ct, pk, ak, brk, B, primes);
            var ctboot = RLWE.BootstrappingNTT(ct, pk, ak, brkntt, B, primes);
            Console.WriteLine($"ct     {ct.Params}");
            Console.WriteLine($"ctboot {ctboot.Params}");

            var m2boot = RLWE.DecryptBGV(ctboot, sk);
            Console.WriteLine($"m      = {m}");
            Console.WriteLine($"m^2    = {m2}");
            Console.WriteLine($"ctboot = {m2boot}");
            Console.WriteLine($"eboot  = {RLWE.ErrorsBGV(ctboot, sk).NormInf()}");
            Console.WriteLine($"emodsw = {RLWE.ErrorsBGV(ct.ModSwitch(qL), sk).NormInf()}");
            Console.WriteLine();
            if (!m2.Equals(m2boot))
                throw new();

            var c2 = RLWE.MulRelinBGV(ctboot, ctboot, rlk).ModSwitch(rlks[qL].nextMod);
            var d2 = RLWE.DecryptBGV(c2, sk);
            var m4 = (m2 * m2).ResModSigned(pm, t);
            Console.WriteLine($"m^4 = {m4}");
            Console.WriteLine($"    = {d2}");

            GlobalStopWatch.Show($"Test[{k + 1}]");
            Console.WriteLine();
            if (!d2.Equals(m4))
                throw new();
        }

        GlobalStopWatch.Show("End");
        Console.WriteLine();
    }

    public static void Example12HomomorphicAdditionWithCarry()
    {
        // Weak parameters

        var bits = 32;
        var rlwe = new RLWE(N: 16, level: 5, bootstrappingMode: true);
        rlwe.Show();

        var fmt = $"{{0,{(int)(bits * double.Log10(2)) + 1}}}";
        string FMT(long a) => string.Format(fmt, a);

        for (int k = 0; k < 10; ++k)
        {
            GlobalStopWatch.AddLap();
            var m1 = DistributionExt.DiceSample(bits, [0, 1]).Append(0).ToArray();
            var m2 = DistributionExt.DiceSample(bits, [0, 1]).Append(0).ToArray();
            var a1 = Convert.ToInt64(m1.Reverse().Glue(), 2);
            var a2 = Convert.ToInt64(m2.Reverse().Glue(), 2);
            var e1 = rlwe.Encrypt(m1);
            var e2 = rlwe.Encrypt(m2);

            var sumi = a1 + a2;
            var add_m1m2 = Convert.ToString(sumi, 2).PadLeft(bits + 1, '0');

            var add_e1e2 = rlwe.ADD(e1, e2);
            var d_add = rlwe.Decrypt(add_e1e2);

            Console.WriteLine($"   0b{m1.Reverse().Glue()} = {FMT(a1)}");
            Console.WriteLine($" + 0b{m2.Reverse().Glue()} = {FMT(a2)}");
            Console.WriteLine($" = 0b{add_m1m2} = {FMT(sumi)}");
            Console.WriteLine($"   0b{d_add.Reverse().Glue()}");

            GlobalStopWatch.Show($"Test[{k + 1}]"); // Time:1.533s
            Console.WriteLine();

            var sumf = Convert.ToInt64(d_add.Reverse().Glue(), 2);
            if (sumi != sumf)
                throw new();
        }
    }
}
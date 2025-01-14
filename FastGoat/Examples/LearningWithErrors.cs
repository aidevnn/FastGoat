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

    public static void Example5HomomorphicAdditionWithCarry()
    {
        // Warning : Weak parameters
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
    // RLWE N=16=2^4, Φ(N)=8 PM=x^8 + 1 t=97 q=9797=97*101
    // Private Key
    // -x^7 - x^6 + x^5 - x^4 + x^2 - x + 1
    // Public Key
    // Q:9797    PM:x^8 + 1
    // A:6088*x^7 + 1227*x^6 + 8782*x^5 + 5579*x^4 + 8969*x^3 + 9703*x^2 + 4517*x + 7252
    // B:5630*x^7 + 2868*x^6 + 8082*x^5 + 4966*x^4 + 5505*x^3 + 2792*x^2 + 340*x + 1703
    // Relinearisation Key
    // Q:9797    PM:x^8 + 1
    // A:2164*x^7 + 4753*x^6 + 5782*x^5 + 212*x^4 + 7543*x^3 + 8311*x^2 + 1794*x + 2306
    // B:5596*x^7 + 8088*x^6 + 4020*x^5 + 1749*x^4 + 5929*x^3 + 2678*x^2 + 1211*x + 4992
    // 
    // Encrypt(m1 = 0b01111100111000101101111001011110)
    //     (0, (A:974*x^7 + 3638*x^6 + 6688*x^5 + 3975*x^4 + 6790*x^3 + 8483*x^2 + 3850*x + 1649, B:5946*x^7 + 2842*x^6 + 2600*x^5 + 3177*x^4 + 7970*x^3 + 4647*x^2 + 4383*x + 5851) mod 9797 mod x^8 + 1)
    //     (1, (A:8107*x^7 + 7240*x^6 + 6153*x^5 + 742*x^4 + 1646*x^3 + 6334*x^2 + 1550*x + 8755, B:771*x^7 + 1529*x^6 + 67*x^5 + 5779*x^4 + 9282*x^3 + 9316*x^2 + 5807*x + 5896) mod 9797 mod x^8 + 1)
    //     (1, (A:2000*x^7 + 5804*x^6 + 4744*x^5 + 2375*x^4 + 7649*x^3 + 6255*x^2 + 2733*x + 1716, B:310*x^7 + 8288*x^6 + 8174*x^5 + 9600*x^4 + 1924*x^3 + 7796*x^2 + 3822*x + 9655) mod 9797 mod x^8 + 1)
    //     (1, (A:2377*x^7 + 4378*x^6 + 6299*x^5 + 1592*x^4 + 9277*x^3 + 571*x^2 + 6656*x + 4430, B:4794*x^7 + 3019*x^6 + 9395*x^5 + 3219*x^4 + 7051*x^3 + 2494*x^2 + 9731*x + 4496) mod 9797 mod x^8 + 1)
    //     (1, (A:1982*x^7 + 9101*x^6 + 3521*x^5 + 5470*x^4 + 1195*x^3 + 9179*x^2 + 2388*x + 6504, B:7978*x^7 + 1622*x^6 + 4811*x^5 + 166*x^4 + 6357*x^3 + 1877*x^2 + 5752*x + 9041) mod 9797 mod x^8 + 1)
    //     (0, (A:2684*x^7 + 6212*x^6 + 8654*x^5 + 8433*x^4 + 9286*x^3 + 3213*x^2 + 7426*x + 9451, B:6799*x^7 + 6249*x^6 + 3611*x^5 + 6789*x^4 + 4918*x^3 + 5963*x^2 + 9372*x + 3141) mod 9797 mod x^8 + 1)
    // ...
    // Encrypt(m2 = 0b01000101001011000010100000110110)
    //     (0, (A:2847*x^7 + 6301*x^6 + 2166*x^5 + 7408*x^4 + 9144*x^3 + 7170*x^2 + 6106*x + 6937, B:9424*x^7 + 5278*x^6 + 1636*x^5 + 8713*x^4 + 9213*x^3 + 9392*x^2 + 3849*x + 8965) mod 9797 mod x^8 + 1)
    //     (1, (A:1229*x^7 + 1609*x^6 + 2837*x^5 + 5988*x^4 + 7812*x^3 + 23*x^2 + 2353*x + 3178, B:889*x^7 + 4016*x^6 + 9456*x^5 + 1196*x^4 + 9716*x^3 + 4410*x^2 + 5554*x + 4716) mod 9797 mod x^8 + 1)
    //     (1, (A:2757*x^7 + 1265*x^6 + 5973*x^5 + 3098*x^4 + 3669*x^3 + 8691*x^2 + 2466*x + 3253, B:3093*x^7 + 8978*x^6 + 8302*x^5 + 8893*x^4 + 8386*x^3 + 4144*x^2 + 8778*x + 5749) mod 9797 mod x^8 + 1)
    //     (0, (A:2750*x^7 + 6495*x^6 + 2166*x^5 + 7117*x^4 + 8853*x^3 + 7655*x^2 + 5815*x + 6743, B:9424*x^7 + 5569*x^6 + 1442*x^5 + 8616*x^4 + 9116*x^3 + 9683*x^2 + 3655*x + 9062) mod 9797 mod x^8 + 1)
    //     (1, (A:4325*x^7 + 1556*x^6 + 870*x^5 + 7013*x^4 + 8962*x^3 + 3688*x^2 + 7646*x + 1786, B:7484*x^7 + 4018*x^6 + 7433*x^5 + 594*x^4 + 3596*x^3 + 4386*x^2 + 8832*x + 7747) mod 9797 mod x^8 + 1)
    //     (1, (A:8038*x^7 + 3485*x^6 + 6035*x^5 + 4297*x^4 + 6890*x^3 + 7520*x^2 + 1356*x + 5973, B:7484*x^7 + 950*x^6 + 5652*x^5 + 1305*x^4 + 1003*x^3 + 7090*x^2 + 8992*x + 1922) mod 9797 mod x^8 + 1)
    // ...
    // Decrypt(e1 + e2) = 0b11000010000011110000011010010100
    //     (0, (A:379*x^7 + 2246*x^6 + 9450*x^5 + 6143*x^4 + 1352*x^3 + 1388*x^2 + 9367*x + 9080, B:966*x^7 + 9729*x^6 + 628*x^5 + 5514*x^4 + 6903*x^3 + 3696*x^2 + 6295*x + 220) mod 9797 mod x^8 + 1)
    //     (0, (A:1559*x^7 + 7804*x^6 + 4999*x^5 + 2071*x^4 + 8291*x^3 + 957*x^2 + 2624*x + 4906, B:8993*x^7 + 17*x^6 + 2849*x^5 + 3069*x^4 + 943*x^3 + 7421*x^2 + 9672*x + 9012) mod 9797 mod x^8 + 1)
    //     (1, (A:2312*x^7 + 9396*x^6 + 9749*x^5 + 5398*x^4 + 7022*x^3 + 4041*x^2 + 9500*x + 2783, B:4623*x^7 + 7283*x^6 + 8823*x^5 + 3559*x^4 + 8462*x^3 + 810*x^2 + 3422*x + 3932) mod 9797 mod x^8 + 1)
    //     (0, (A:1061*x^7 + 613*x^6 + 391*x^5 + 8703*x^4 + 9704*x^3 + 8413*x^2 + 5361*x + 817, B:9138*x^7 + 5430*x^6 + 1416*x^5 + 9144*x^4 + 8298*x^3 + 3062*x^2 + 9018*x + 4747) mod 9797 mod x^8 + 1)
    //     (1, (A:6766*x^7 + 1304*x^6 + 8624*x^5 + 7040*x^4 + 3137*x^3 + 7700*x^2 + 1090*x + 5103, B:3624*x^7 + 637*x^6 + 741*x^5 + 9470*x^4 + 5951*x^3 + 6679*x^2 + 2113*x + 5083) mod 9797 mod x^8 + 1)
    //     (0, (A:9477*x^7 + 3214*x^6 + 531*x^5 + 5959*x^4 + 3403*x^3 + 328*x^2 + 2847*x + 7655, B:7945*x^7 + 3745*x^6 + 279*x^5 + 4879*x^4 + 9559*x^3 + 5914*x^2 + 5144*x + 1830) mod 9797 mod x^8 + 1)
    // ...
    //    0b01111100111000101101111001011110 = 2095242846
    //  + 0b01000101001011000010100000110110 = 1160521782
    //  = 0b11000010000011110000011010010100 = 3255764628
    //    0b11000010000011110000011010010100
    // 

    public static void Example6HomomorphicMultiplicationWithCarry()
    {
        // Warning : Weak parameters
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
        // Warning : Weak parameters
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
        var eia = rlwe.Errors(c1a).Concat(rlwe.Errors(c2a)).Select(e => e.Signed(t).Absolute).Max();
        var efa = xa.Select(e => rlwe.Errors(e).Signed(t).Absolute).Max();
        
        var c1b = rlwe.FromRegevCipher(reg.Encrypt(m1));
        var c2b = rlwe.FromRegevCipher(reg.Encrypt(m2));
        var xb = rlwe.XOR(
            rlwe.XOR(c1b.Take(rlwe.n / 2).ToArray(), c1b.Skip(rlwe.n / 2).ToArray()),
            rlwe.XOR(c2b.Take(rlwe.n / 2).ToArray(), c2b.Skip(rlwe.n / 2).ToArray())
        );
        var eib = rlwe.Errors(c1b).Concat(rlwe.Errors(c2b)).Select(e => e.Signed(t).Absolute).Max();
        var efb = xb.Select(e => rlwe.Errors(e).Signed(t).Absolute).Max();

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
}
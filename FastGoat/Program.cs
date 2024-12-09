using System.Numerics;
using FastGoat.Commons;
using static FastGoat.Commons.IntExt;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");
// 16-bit default C#-char
const int charBit = 8;

IEnumerable<int> String2Bin(string s)
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

string Bin2String(IEnumerable<int> bin)
{
    var s = "";
    foreach (var l in bin.Chunk(charBit))
    {
        var sum = l.Select((c, i) => (c, i)).Sum(e => e.c << e.i);
        s += (char)sum;
    }

    return s;
}

void RunLWE(LWE lwe, string text, bool showCypher = false, bool showBinary = true)
{
    Console.WriteLine(text);
    var seq = String2Bin(text).ToArray();
    var seqCyphers = seq.Select(b => (b, cypher: lwe.EncryptBit(b))).ToArray();

    if (showCypher)
        seqCyphers.Println(
            l => $"{l.b} => [[{l.cypher.Ai.Glue(",", lwe.Fmt)}], {string.Format(lwe.Fmt, l.cypher.Bi)}]",
            $"Cyphers text:{text}"
        );

    var seqDecrypt = seqCyphers.Select(e => lwe.DecryptBit(e.cypher)).ToArray();
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

string RandString(int length)
{
    var s = "abcdefghijklmnopqrstuvwxyzabcdefghijklmnopqrstuvwxyzabcdefghijklmnopqrstuvwxyz" +
            @"ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789+-*/!@#$%^&*()_+:<>?|\                   ";
    return length.Range().Select(_ => s[Rng.Next(s.Length)]).Glue();
}

{
    // RngSeed(95782);
    
    for (int N = 2; N < 21; ++N)
    {
        var lwe = new LWE(N);
        lwe.Show();

        RunLWE(lwe, "hello world lwe");
        RunLWE(lwe, "Hello World LWE");
        RunLWE(lwe, "AAA+", showCypher: true);

        for (int i = 0; i < 100; i++)
            RunLWE(lwe, RandString(Rng.Next(20, 50)));

        // long text
        RunLWE(lwe, RandString(10000), showBinary: false);
    }
}

/// <summary>
/// On Lattices, Learning with Errors,
/// Random Linear Codes, and Cryptography
/// Oded Regev
/// January 9, 2024
/// https://arxiv.org/abs/2401.03703
/// 
/// #5 Public Key Cryptosystem
/// </summary>
public class LWE
{
    public (int[] Ai, int Bi)[] PK { get; }
    private int[] SK { get; }
    public int P { get; }
    public int M { get; }
    public int N { get; }
    public int Err { get; }
    public double A { get; }
    public string Fmt { get; }

    static int MulModP(long a, long b, int p) => (int)((a * b) % p);

    static int AModP(int a, int p)
    {
        int r = a % p;
        return r < 0 ? r + p : r;
    }

    public LWE(int n)
    {
        if (n < 2 || n > 20)
            throw new($"N = {n} must be > 1 and < 21");

        M = (int)(2 * n * double.Log2(n));
        Err = (int)(M - M / (2 * double.Sqrt(n)));
        var p = Primes10000.First(p => p > 4 * Err * Err); // Magik prime
        var a = 1 / (double.Log2(n) * double.Sqrt(n));
        (P, N, A) = (p, n, a);
        SK = N.Range().Select(_ => Rng.Next(1, p)).ToArray();
        var PKai = M.Range().Select(_ => n.Range().Select(_ => Rng.Next(1, p)).ToArray()).ToArray();
        PK = PKai.Select(
                Ai => (Ai, Bi: AModP(Ai.Zip(SK).Sum(f => MulModP(f.First, f.Second, p)) + Rng.Next(1, Err), p)))
            .ToArray();

        var nbDigits = $"{P + 1}".Length;
        Fmt = $"{{0,{nbDigits}}}";
    }

    public int DecryptBit((int[] Ai, int Bi) cypher)
    {
        var c = AModP(cypher.Bi - AModP(cypher.Ai.Zip(SK).Sum(e => MulModP(e.First, e.Second, P)), P), P);
        return c < P / 2 ? 0 : 1;
    }

    public (int[] Ai, int Bi) EncryptBit(int b)
    {
        var m = PK.Length;
        var n = PK[0].Ai.Length;
        var r = Rng.Next(2, m);
        var set = r.Range().Select(_ => PK[Rng.Next(m)]).ToArray();
        var sumBi = AModP(set.Sum(e => e.Bi), P);
        var sumAi = n.Range().Select(i => AModP(set.Sum(e => e.Ai[i]), P)).ToArray();
        if (b != 0)
            sumBi = AModP(sumBi + P / 2, P);

        return (sumAi, sumBi);
    }

    public void Show()
    {
        Console.WriteLine($"P:{P} N:{N} M:{M} Err:{Err} A:{A:G6}");
        Console.WriteLine("Private Key");
        Console.WriteLine($"[{SK.Glue(", ", Fmt)}]");
        PK.Println(e => $"[[{e.Ai.Glue(", ", Fmt)}], {string.Format(Fmt, e.Bi)}]", "Public Key");
        Console.WriteLine();
    }
}
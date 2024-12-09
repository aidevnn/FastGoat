using FastGoat.Commons;
using static FastGoat.Commons.IntExt;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");
// 16-bit default C#-char
const int charBit = 16;

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

void RunLWE(LWE lwe, string text, bool showCypher = false)
{
    var seq = String2Bin(text).ToArray();
    var seqCyphers = seq.Select(b => (b, cypher: lwe.EncryptBit(b))).ToArray();

    if (showCypher)
        seqCyphers.Println(l => $"{l.b} => [[{l.cypher.Ai.Glue(",", "{0,3}")}], {l.cypher.Bi,3}]", $"Cyphers text:{text}");
    
    var seqDecrypt = seqCyphers.Select(e => lwe.DecryptBit(e.cypher)).ToArray();
    Console.WriteLine(text);
    Console.WriteLine($"seqInput  :[{seq.Glue()}]");
    Console.WriteLine($"seqDecrypt:[{seqDecrypt.Glue()}]");
    var text2 = Bin2String(seqDecrypt);
    Console.WriteLine(text2);
    if(string.Equals(text,text2))
        Console.WriteLine("    SUCCESS");
    else
    {
        Console.WriteLine("    FAIL");
        Console.Beep();
    }
    
    Console.WriteLine();
}

{
    // RngSeed(95782);
    
    var lwe = new LWE(613, 5);
    lwe.Show();
    
    RunLWE(lwe, "hello world lwe");
    RunLWE(lwe, "Hello World LWE");
    RunLWE(lwe, "AAA+", showCypher: true);
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
    public int N { get; }
    public int Err { get; }
    public double A { get; }

    public LWE(int p, int n)
    {
        if (n < 3)
            throw new($"N = {n} must be > 2");

        var a = 1 / (double.Log2(n) * double.Sqrt(n));
        (P, N, A) = (p, n, a);
        var M = (int)(2 * n * double.Log2(n));
        Err = (int)(M - M / (2 * double.Sqrt(n)));
        SK = N.Range().Select(_ => Rng.Next(1, p)).ToArray();
        var PKai = M.Range().Select(_ => n.Range().Select(_ => Rng.Next(1, p)).ToArray()).ToArray();
        PK = PKai.Select(Ai => (Ai, Bi: AmodP(Ai.Zip(SK).Sum(f => f.First * f.Second) + Rng.Next(1, Err), p))).ToArray();
    }

    public int DecryptBit((int[] Ai, int Bi) cypher)
    {
        var c = AmodP(cypher.Bi - AmodP(cypher.Ai.Zip(SK).Sum(e => e.First * e.Second), P), P);
        return c < P / 2 ? 0 : 1;
    }

    public (int[] Ai, int Bi) EncryptBit(int b)
    {
        var m = PK.Length;
        var n = PK[0].Ai.Length;
        var r = Rng.Next(2, m);
        var set = r.Range().Select(_ => PK[Rng.Next(m)]).ToArray();
        var sumBi = AmodP(set.Sum(e => e.Bi), P);
        var sumAi = n.Range().Select(i => AmodP(set.Sum(e => e.Ai[i]), P)).ToArray();
        if (b != 0)
            sumBi = AmodP(sumBi + P / 2, P);

        return (sumAi, sumBi);
    }

    public void Show()
    {
        Console.WriteLine($"P:{P} N:{N} Err:{Err} A:{A:G6}");
        var nbDigits = $"{P + 1}".Length;
        var fmt = $"{{0,{nbDigits}}}";
        Console.WriteLine("Private Key");
        Console.WriteLine($"[{SK.Glue(", ", fmt)}]");
        PK.Println(e => $"[[{e.Ai.Glue(", ", fmt)}], {string.Format(fmt, e.Bi)}]", "Public Key");
        Console.WriteLine();
    }
}
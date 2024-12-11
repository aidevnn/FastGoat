using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Lattice;
using static FastGoat.Commons.IntExt;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");
// 16-bit default C#-char
int charBit = 8;

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
        s += Convert.ToChar(sum);
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
            l => $"{l.b} => [[{l.cypher.Ai.Glue(",", lwe.Fmt)}], {string.Format(lwe.Fmt, l.cypher.B)}]",
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

void RunRLWE(RLWE rlwe, string text, bool showCypher = false, bool showBinary = true)
{
    Console.WriteLine(text);
    var seq0 = String2Bin(text).ToArray();
    var rem = seq0.Length % rlwe.N;
    var seq1 = rem == 0 ? seq0 : seq0.Concat(Enumerable.Repeat(0, rlwe.N - rem)).ToArray();
    var diff = seq1.Length - seq0.Length;
    var seq = seq1.Chunk(rlwe.N).ToArray();
    var seqCyphers = seq.Select(b => (b, encode: rlwe.Encode(b)))
        .Select(e => (e.b, e.encode, cypher: rlwe.Encrypt(e.encode)))
        .ToArray();

    if (showCypher)
        seqCyphers.Println(
            l => $"{l.b.Glue()} => [[{l.encode}], {l.cypher}]",
            $"Cyphers text:{text}"
        );

    var seqDecrypt = seqCyphers.Select(e => rlwe.Decode(rlwe.Decrypt(e.cypher))).ToArray();
    var text2 = Bin2String(seqDecrypt.SelectMany(l => l).SkipLast(diff));
    if (showBinary)
    {
        Console.WriteLine($"seqInput  :[{seq.SelectMany(e => e).SkipLast(diff).Glue()}]");
        Console.WriteLine($"seqDecrypt:[{seqDecrypt.SelectMany(e => e).SkipLast(diff).Glue()}]");
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

{
    Ring.DisplayPolynomial = MonomDisplay.StarCaret;
    for (var k = 3; k < 7; ++k)
    {
        var N = 1 << k; // 2^3, 2^4, 2^5, 2^6
        var Q = Primes10000.First(p => p > N * N & p % (2 * N) == 1);
        var rlwe = new RLWE(N, Q, sigma: double.Sqrt(N) / k);
        rlwe.Show();

        RunRLWE(rlwe, "hello world lwe");
        RunRLWE(rlwe, "Hello World LWE");
        RunRLWE(rlwe, "AAA+", showCypher: true);

        for (int i = 0; i < 100; i++)
            RunRLWE(rlwe, RandString(Rng.Next(20, 50)));

        // long text
        RunRLWE(rlwe, RandString(10000), showBinary: false);
    }
}
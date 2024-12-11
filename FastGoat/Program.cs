using FastGoat.Commons;
using FastGoat.UserGroup.Lattice;
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

void TestLWE()
{
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

{
    TestLWE();

    var lwe = new LWE(m: 5, n: 20); // PK sample size 5 and dimension 20
    lwe.Show();
    
    RunLWE(lwe, "hello world lwe");
    RunLWE(lwe, "Hello World LWE");
    RunLWE(lwe, "AAA+", showCypher: true);
    
    for (int i = 0; i < 100; i++)
        RunLWE(lwe, RandString(Rng.Next(20, 50)));
    
    // long text
    RunLWE(lwe, RandString(10000), showBinary: false);
}
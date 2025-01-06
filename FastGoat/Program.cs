using System.Numerics;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Lattice;
using Rq = FastGoat.Structures.VecSpace.KPoly<FastGoat.UserGroup.Integers.Rational>;
using static FastGoat.Commons.IntExt;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");
// RngSeed(57419);

void regev128()
{
    var lweRegev = new Regev(128);
    Console.WriteLine(lweRegev.Params);
    var bits = DistributionExt.BernouilliSample(100, 0.5).ToArray();
    var encrypted = bits.Select(e => lweRegev.EncryptBit(e)).ToArray();
    var decrypted = encrypted.Select(e => lweRegev.DecryptBit(e)).ToArray();
    Console.WriteLine($"bits     :[{bits.Glue(", ")}]");
    Console.WriteLine($"decrypted:[{decrypted.Glue(", ")}]");
    Console.WriteLine($"Are equal:{bits.SequenceEqual(decrypted)}");
    lweRegev.Err.GroupBy(e => e.Signed).Select(e => new { e = e.Key, Nb = e.Count() })
        .OrderBy(e => e.e).Println("Err sample");
}

int[] BinarySum(int[] a, int[] b)
{
    var carry = 0;
    var sum = new List<int>();
    var length = int.Max(a.Length, b.Length);
    var A = length.SeqLazy().Select(i => i < a.Length ? a[i] : 0);
    var B = length.SeqLazy().Select(i => i < b.Length ? b[i] : 0);
    
    var trackCarry = new List<int>();
    foreach (var (ai, bi) in A.Zip(B))
    {
        sum.Add(ai ^ bi ^ carry);
        carry &= ai ^ bi;
        carry ^= ai & bi;
        trackCarry.Add(carry);
    }

    Console.WriteLine($"c: {trackCarry.Glue()}");
    return sum.ToArray();
}

CipherLWE[] BinarySumLWE(LWE lwe, CipherLWE[] a, CipherLWE[] b)
{
    var carry = lwe.Zero;
    var sum = new List<CipherLWE>();
    var length = int.Max(a.Length, b.Length);
    var A = length.SeqLazy().Select(i => i < a.Length ? a[i] : lwe.Zero);
    var B = length.SeqLazy().Select(i => i < b.Length ? b[i] : lwe.Zero);
    var ek = lwe.EK;
    
    var trackCarry = new List<CipherLWE>();
    foreach (var (ai, bi) in A.Zip(B))
    {
        var ai_xor_bi = CipherLWE.Xor(ai, bi);
        sum.Add(CipherLWE.Xor(carry, ai_xor_bi));
        
        var ai_and_bi = CipherLWE.And(ai, bi, ek);
        carry = CipherLWE.Xor(CipherLWE.And(carry, ai_xor_bi, ek), ai_and_bi);
        // carry = lwe.BrutalRefresh(carry);
        trackCarry.Add(carry);
    }

    Console.WriteLine($"c: {trackCarry.Select(lwe.DecryptBit).Glue()} lwe carry");
    Console.WriteLine($"Err:{trackCarry.Select(lwe.Error).Glue(", ", "{0,3}")}");
    return sum.ToArray();
}

int Bin2Int(int[] l) => l.Reverse().Aggregate(0, (acc, i) => acc * 2 + i);

void TestAddition(int n, int bits, int nbTrials = 100, bool noiseMode = true)
{
    var lwe = new LWE(n, q: bits * n.Pow(3), m: n, noiseMode);
    lwe.Show();

    var set = new List<int>();
    for (int k = 0; k < nbTrials; ++k)
    {
        Console.WriteLine($"Step:{k + 1}");
        var a = DistributionExt.DiceSample(bits, [0, 1]).ToArray();
        var b = DistributionExt.DiceSample(bits, [0, 1]).ToArray();
        var s = BinarySum(a, b);
        Console.WriteLine($"a:{a.Glue()}");
        Console.WriteLine($"b:{b.Glue()}");
        Console.WriteLine($"+:{s.Glue()}");

        var an = Bin2Int(a);
        var bn = Bin2Int(b);
        var sn = Bin2Int(s);
        var anbn = (an + bn) % (1 << bits);

        var la = a.Select(ai => lwe.EncryptBit(ai)).ToArray();
        var lb = b.Select(bi => lwe.EncryptBit(bi)).ToArray();
        var ls = BinarySumLWE(lwe, la, lb);
        var ds = ls.Select(lwe.DecryptBit).ToArray();
        var dsn = Bin2Int(ds);
        
        Console.WriteLine($"a:{an} + b:{bn} = {anbn} s:{sn} ds:{dsn}");
        if (sn != anbn)
            throw new(); // never reached
        if (sn != dsn)
        {
            Console.WriteLine($" :{ds.Glue()} Err step[{k}]");
            Console.WriteLine($"err a:{la.Select(lwe.Error).Glue(", ", "{0,3}")}");
            Console.WriteLine($"err b:{lb.Select(lwe.Error).Glue(", ", "{0,3}")}");
            Console.WriteLine(lwe.Params);
            set.Add(k + 1);
        }

        Console.WriteLine();
    }

    Console.WriteLine($"Success {nbTrials - set.Count}/{nbTrials} additions with carry and {lwe.Params}");
    if (set.Any())
        Console.WriteLine($"Errors at Steps:{set.Glue(", ")}");
    
    Console.WriteLine();
}

{
    // TestAddition(n: 8, bits: 8, nbTrials: 1000, noiseMode: false);
    // TestAddition(n: 8, bits: 16, nbTrials: 1000, noiseMode: false);
    // TestAddition(n: 8, bits: 8, nbTrials: 1000);
    TestAddition(n: 8, bits: 16, nbTrials: 1000);
    
    // Step:1000
    // c: 1001101110111110
    // a:1011101010111010
    // b:1001101110110100
    // +:0110110011010001
    // c: 1001101110111110 lwe carry
    // Err:  4,  -1,   4,  -2,   4,  -1,  -3,   6,  14, -28,  13,  -3,  -2,   2,   4,   0
    // a:23901 + b:11737 = 35638 s:35638 ds:35638
    // 
    // Success 999/1000 additions with carry and LWE Params N:8 M:8 Q:8192
    // Errors at Steps:481
    // 
}
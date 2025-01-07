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

long Bin2Int(int[] l) => l.Reverse().Aggregate((long)0, (acc, i) => acc * 2 + i);

void TestAddition(int n, long q, int m, int bits, int nbTrials = 100, bool noiseMode = true)
{
    var lwe = new LWE(n, q, m, noiseMode);
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
        var anbn = (an + bn) % ((long)1 << bits);

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
    var n = 8; // for faster test
    var bits = 32;
    var m = n; // more errors conduct to bigger prime q
    var q = Primes10000.First(q => q > (1 << 13));
    var nbTrials = 1000;
    TestAddition(n, q, m, bits, nbTrials);
    
    // Step:1000
    // c: 01111111111111011111100000000000
    // a:01001100100010011111100000100101
    // b:11111011011111011010000111001000
    // +:10001000000010101010010111101101
    // c: 01111111111111011111100000000000 lwe carry
    // Err: -3,  -4,  -4,   7,  -1,   6,   0,  -3, -12,  -6,  -8,   0,  -3,  -1,   0,  13,   3,  -1,  -5,   2,   0,   3,  -9, -19,  55, -54,  55, -108, 108, 108,  -1,   8
    // a:2753532210 + b:327532255 = 3081064465 s:3081064465 ds:3081064465
    // 
    // Success 996/1000 additions with carry and LWE Params N:8 M:8 Q:8209
    // Errors at Steps:187, 325, 533, 729
    // 
}
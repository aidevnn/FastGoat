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

    if (Logger.Level != LogLevel.Off)
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

    if (Logger.Level != LogLevel.Off)
    {
        Console.WriteLine($"c: {trackCarry.Select(lwe.DecryptBit).Glue()} lwe carry");
        Console.WriteLine($"Err:{trackCarry.Select(lwe.Error).Glue(", ", "{0,3}")}");
    }

    return sum.ToArray();
}

long Bin2Int(int[] l) => l.Reverse().Aggregate((long)0, (acc, i) => acc * 2 + i);

void TestAddition(int n, int l, int bits, int nbTrials = 100)
{
    var set = new List<int>();
    var lweparams = "";
    for (int k = 0; k < nbTrials; ++k)
    {
        // regenerating all keys for each test
        var lwe = LWE.Setup(n, l);

        if (k == 0)
            lweparams = lwe.Params;

        var a = DistributionExt.DiceSample(bits, [0, 1]).ToArray();
        var b = DistributionExt.DiceSample(bits, [0, 1]).ToArray();
        var s = BinarySum(a, b);

        var an = Bin2Int(a);
        var bn = Bin2Int(b);
        var sn = Bin2Int(s);
        var anbn = (an + bn) % ((long)1 << bits);

        var la = a.Select(ai => lwe.EncryptBit(ai)).ToArray();
        var lb = b.Select(bi => lwe.EncryptBit(bi)).ToArray();
        var ls = BinarySumLWE(lwe, la, lb);
        var ds = ls.Select(lwe.DecryptBit).ToArray();
        var dsn = Bin2Int(ds);

        if (Logger.Level != LogLevel.Off)
        {
            Console.WriteLine($"Step:{k + 1}");
            Console.WriteLine($"a:{a.Glue()}");
            Console.WriteLine($"b:{b.Glue()}");
            Console.WriteLine($"+:{s.Glue()}");
            Console.WriteLine($"a:{an} + b:{bn} = {anbn} s:{sn} ds:{dsn}");
        }

        if (sn != anbn)
            throw new(); // never reached

        if (sn != dsn)
        {
            if (Logger.Level != LogLevel.Off)
            {
                Console.WriteLine($" :{ds.Glue()} Err step[{k}]");
                Console.WriteLine($"err a:{la.Select(lwe.Error).Glue(", ", "{0,3}")}");
                Console.WriteLine($"err b:{lb.Select(lwe.Error).Glue(", ", "{0,3}")}");
                Console.WriteLine(lwe.Params);
            }

            set.Add(k + 1);
        }

        if (Logger.Level != LogLevel.Off)
            Console.WriteLine();
    }

    Console.WriteLine(
        $"Tests {nbTrials - set.Count,4}/{nbTrials} for {bits,2}-bits additions with carry and {lweparams}");

    if (set.Any() && Logger.Level != LogLevel.Off)
        Console.WriteLine($"Errors at Steps:{set.Glue(", ")}");
}

{
    Logger.SetOff();
    GlobalStopWatch.Restart();
    GlobalStopWatch.AddLap();
    RecomputeAllPrimesUpTo(5000000);
    GlobalStopWatch.Show($"Primes {Primes10000.Count}");
    
    var bitsLevel = 4.SeqLazy(1).Select(i => 2 << i).ToArray();
    var nLevel = 4.SeqLazy(2).Select(i => 1 << i).ToArray();
    var qLevel = 2.Range(start: 2, step: 8);
    
    // qLevel.Grid3D(bitsLevel, nLevel).ToList().ForEach(e => TestAddition(e.t3, e.t1, e.t2, nbTrials));

    var nbTrials = 100;
    GlobalStopWatch.AddLap();
    foreach (var n in nLevel)
    {
        GlobalStopWatch.AddLap();
        foreach (var bits in bitsLevel)
        {
            foreach (var ql in qLevel)
                TestAddition(n, ql, bits, nbTrials);

            Console.WriteLine();
        }
        
        GlobalStopWatch.Show($"End N:{n}");
        Console.WriteLine();
    }

    GlobalStopWatch.Show();
    Console.Beep();

    // Tests  100/100 for  4-bits additions with carry and LWE Params N:4 M:38 Q:131
    // Tests  100/100 for  4-bits additions with carry and LWE Params N:4 M:82 Q:32771
    // 
    // Tests   89/100 for  8-bits additions with carry and LWE Params N:4 M:38 Q:131
    // Tests  100/100 for  8-bits additions with carry and LWE Params N:4 M:82 Q:32771
    // 
    // Tests   64/100 for 16-bits additions with carry and LWE Params N:4 M:38 Q:131
    // Tests  100/100 for 16-bits additions with carry and LWE Params N:4 M:82 Q:32771
    // 
    // Tests   46/100 for 32-bits additions with carry and LWE Params N:4 M:38 Q:131
    // Tests  100/100 for 32-bits additions with carry and LWE Params N:4 M:82 Q:32771
    // 
    // # End N:4 Time:2.274s
    // 
    // Tests   98/100 for  4-bits additions with carry and LWE Params N:8 M:90 Q:577
    // Tests  100/100 for  4-bits additions with carry and LWE Params N:8 M:169 Q:147457
    // 
    // Tests   90/100 for  8-bits additions with carry and LWE Params N:8 M:90 Q:577
    // Tests  100/100 for  8-bits additions with carry and LWE Params N:8 M:169 Q:147457
    // 
    // Tests   70/100 for 16-bits additions with carry and LWE Params N:8 M:90 Q:577
    // Tests  100/100 for 16-bits additions with carry and LWE Params N:8 M:169 Q:147457
    // 
    // Tests   41/100 for 32-bits additions with carry and LWE Params N:8 M:90 Q:577
    // Tests   98/100 for 32-bits additions with carry and LWE Params N:8 M:169 Q:147457
    // 
    // # End N:8 Time:4.744s
    // 
    // Tests   99/100 for  4-bits additions with carry and LWE Params N:16 M:205 Q:2053
    // Tests  100/100 for  4-bits additions with carry and LWE Params N:16 M:355 Q:524309
    // 
    // Tests   76/100 for  8-bits additions with carry and LWE Params N:16 M:205 Q:2053
    // Tests  100/100 for  8-bits additions with carry and LWE Params N:16 M:355 Q:524309
    // 
    // Tests   45/100 for 16-bits additions with carry and LWE Params N:16 M:205 Q:2053
    // Tests   98/100 for 16-bits additions with carry and LWE Params N:16 M:355 Q:524309
    // 
    // Tests   20/100 for 32-bits additions with carry and LWE Params N:16 M:205 Q:2053
    // Tests   87/100 for 32-bits additions with carry and LWE Params N:16 M:355 Q:524309
    // 
    // # End N:16 Time:15.356s
    // 
    // Tests   99/100 for  4-bits additions with carry and LWE Params N:32 M:459 Q:6421
    // Tests  100/100 for  4-bits additions with carry and LWE Params N:32 M:749 Q:1638431
    // 
    // Tests   66/100 for  8-bits additions with carry and LWE Params N:32 M:459 Q:6421
    // Tests  100/100 for  8-bits additions with carry and LWE Params N:32 M:749 Q:1638431
    // 
    // Tests   19/100 for 16-bits additions with carry and LWE Params N:32 M:459 Q:6421
    // Tests   79/100 for 16-bits additions with carry and LWE Params N:32 M:749 Q:1638431
    // 
    // Tests    2/100 for 32-bits additions with carry and LWE Params N:32 M:459 Q:6421
    // Tests   51/100 for 32-bits additions with carry and LWE Params N:32 M:749 Q:1638431
    // 
    // # End N:32 Time:1m17s
    // 
}
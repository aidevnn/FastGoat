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
Ring.DisplayPolynomial = MonomDisplay.StarCaret;

int OpBool(string name, int c1, int c2)
{
    if (name == "and")
        return c1 & c2;
    if (name == "nand")
        return 1 - (c1 & c2);
    if (name == "nor")
        return 1 - (c1 | c2);
    if (name == "or")
        return c1 | c2;
    else
        return c1 ^ c2;
}

long Bin2Int(int[] l) => l.Reverse().Aggregate((long)0, (acc, i) => acc * 2 + i);

void add(int[] f1, int[] f2)
{
    var l = f1.Length;
    var a0 = 0;
    for (int i = 0; i < l; i++)
    {
        if (i == 1) a0 = f2[0];
        if (i > 1) (f2[i - 1], a0) = (a0, f2[i - 1]);
        var (e1, e2) = (f1[i], f2[i]);
        (f1[i], f2[i]) = (e1 ^ e2, e1 & e2);
    }
    
    f2[l - 1] = a0;
    f2[0] = 0;
}

void TestEncryptDecrypt(int n, int nbTrials = 50)
{
    var rlwe = new RLWE(n);
    rlwe.Show();

    for (int k = 0; k < nbTrials; ++k)
    {
        var mi = DistributionExt.DiceSample(n, [0, 1]).ToArray();
        var ec = rlwe.Encrypt(mi);
        var mf = rlwe.Decrypt(ec);
        Console.WriteLine($"mi:[{mi.Glue(", ")}]");
        Console.WriteLine($"  :[{mf.Glue(", ")}]");
        Console.WriteLine();
        if (!mi.SequenceEqual(mf))
            throw new();
    }
}

void RLWENot(int n, int nbTrials = 50)
{
    var rlwe = new RLWE(n);
    rlwe.Show();

    for (int k = 0; k < nbTrials; ++k)
    {
        var m1 = DistributionExt.DiceSample(n, [0, 1]).ToArray();
        var c1 = rlwe.Encrypt(m1);
        
        var nmi = m1.Select(e => 1 - e).ToArray();
        var nc1 = rlwe.NOT(c1);
        var nmf = rlwe.Decrypt(nc1);
        
        Console.WriteLine($"m   :[{m1.Glue()}]");
        Console.WriteLine($"!m  :[{nmi.Glue()}]");
        Console.WriteLine($"    :[{nmf.Glue()}]");
        if (!nmi.SequenceEqual(nmf))
            throw new($"step[{k}]");
        
        Console.WriteLine();
    }
}

void RLWEGates(int n, string name, int nbTrials = 50)
{
    var rlwe = new RLWE(n);
    rlwe.Show();

    for (int k = 0; k < nbTrials; ++k)
    {
        var m1 = DistributionExt.DiceSample(n, [0, 1]).ToArray();
        var c1 = rlwe.Encrypt(m1);
        var m2 = DistributionExt.DiceSample(n, [0, 1]).ToArray();
        var c2 = rlwe.Encrypt(m2);

        var m1m2i = m1.Zip(m2).Select(e => OpBool(name, e.First, e.Second)).ToArray();
        var c1c2 = rlwe.OP(name, c1, c2);
        var m1m2f = rlwe.Decrypt(c1c2);

        var opm1m2 = $"{name}(m1,m2)";
        Console.WriteLine($"m1         :[{m1.Glue()}]");
        Console.WriteLine($"m2         :[{m2.Glue()}]");
        Console.WriteLine($"{opm1m2,-11}:[{m1m2i.Glue()}]");
        Console.WriteLine($"           :[{m1m2f.Glue()}]");
        if (!m1m2i.SequenceEqual(m1m2f))
            throw new($"step[{k}]");

        Console.WriteLine();
    }
}

void TestAllRLWEGates()
{
    var n = 32;
    RLWENot(n);
    RLWEGates(n, "and");
    RLWEGates(n, "nand");
    RLWEGates(n, "nor");
    RLWEGates(n, "or");
    RLWEGates(n, "xor");
}

void TestADD(int n, int bits, int nbTrials = 50)
{
    var rlwe = new RLWE(n);
    rlwe.Show();
    var mod = (long)1 << bits;

    GlobalStopWatch.Restart();
    for (int k = 0; k < nbTrials; ++k)
    {
        var m1 = DistributionExt.DiceSample(bits - 1, [0, 1]).Append(0).ToArray();
        var e1 = rlwe.Encrypt(m1);
        var m2 = DistributionExt.DiceSample(bits - 1, [0, 1]).Append(0).ToArray();
        var e2 = rlwe.Encrypt(m2);

        var (f1, f2) = (m1.ToArray(), m2.ToArray());
        for (int i = 0; i < bits; i++)
            add(f1, f2);

        var (a1, a2) = (Bin2Int(m1), Bin2Int(m2));
        var sumi = (a1 + a2) % mod;
        var sumf = Bin2Int(f1) % mod;

        var add_c1c2 = rlwe.ADD(e1, e2);
        var dec = rlwe.Decrypt(add_c1c2);

        var fmt = $"{{0,{(int)double.Log10(mod) + 1}}}";
        string FMT(long a) => string.Format(fmt, a);
        
        Console.WriteLine($"m1   :[{m1.Glue()}] = {FMT(a1)}");
        Console.WriteLine($"m2   :[{m2.Glue()}] = {FMT(a2)}");
        Console.WriteLine($"m1+m2:[{f1.Glue()}] = {FMT(sumi)}");
        Console.WriteLine($"     :[{dec.Glue()}]");
        Console.WriteLine();
        if (sumi != sumf || !dec.SequenceEqual(f1))
            throw new();
    }

    GlobalStopWatch.Show($"END nb test:{nbTrials}");
}

{
    // TestEncryptDecrypt(n: 16);
    
    // TestAllRLWEGates();
    
    TestADD(n: 32, bits: 32, nbTrials: 10);
}
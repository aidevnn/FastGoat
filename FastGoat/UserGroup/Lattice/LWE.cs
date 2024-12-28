using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Integers;

namespace FastGoat.UserGroup.Lattice;

public struct LWE
{
    public int N { get; }
    public int M { get; }
    public int Err { get; }
    public int Q { get; }
    public ZnInt64 One { get; }
    public Vec<Vec<ZnInt64>> PK { get; }
    public Vec<ZnInt64> SK { get; }
    public Vec<Vec<ZnInt64>> EK { get; }

    public LWE(int n) : this(n, n)
    {
    }

    public LWE(int n, int m)
    {
        var q = IntExt.Primes10000.First(q => q > n * n * 16 && q % (2 * n) == 1);
        (N, M, Q) = (n, m, q);
        One = new ZnInt64(q, 1);

        Err = (int)IntExt.SqrtBigInt(q);
        var s1 = DistributionExt.DiceSample(1000, -Err, Err).Take(n - 1).Select(k => new ZnInt64(q, k)).ToKMatrix();
        var e1 = DistributionExt.DiceSample(m, [-1, 1]).Select(k => new ZnInt64(q, k)).ToKMatrix(); // TODO: more noise
        var A = DistributionExt.DiceSample(m * (n - 1), 1, q - 1).Select(k => new ZnInt64(q, k)).ToKMatrix(n - 1);

        var a = s1 * A + e1;
        var pk = KMatrix<ZnInt64>.MergeSameCols(A, a);
        var sk = KMatrix<ZnInt64>.MergeSameRows(s1, new[] { -s1.KOne }.ToKMatrix());
        PK = pk.Rows.Select(r => r.ToVec()).ToVec();
        SK = sk.ToVec();

        var o = new ZnInt64(q * q, 1);
        var sks = SK.Select(Signed).ToArray();
        var sksk = sks.Grid2D(sks).Select(t => -t.t1 * t.t2 * q * o).ToKMatrix(n * n);
        var B = DistributionExt.DiceSample(n * n * (n - 1), [-1, 1]).Select(k => k * q * o).ToKMatrix(n * n);
        
        var e2 = DistributionExt.DiceSample(n * n, [-1, 1]).Select(k => new ZnInt64(q * q, k)).ToKMatrix(n * n);
        var diff = B * s1.T + sksk; //  TODO: more noise
        var noise = (N * N).SeqLazy().Select(_ => ScalarSKZero(n, q * q, sks)).ToArray();
        var ek = KMatrix<ZnInt64>.MergeSameRows(B, diff) + KMatrix<ZnInt64>.MergeSameCols(noise);
        EK = ek.Cols.Select(c => c.ToVec()).ToVec();
    }

    static KMatrix<ZnInt64> ScalarSKZero(int n, int q, int[] sks)
    {
        var o = new ZnInt64(q, 1);
        var l = 1000.SeqLazy().Select(_ => IntExt.Rng.Next(o.P / 10, o.P) * o).Where(k => !k.IsZero()).Take(n - 1)
            .ToArray();
        var prod = sks.Zip(l).Select(t => t.First * t.Second).Aggregate((a0, a1) => a0 + a1);
        return l.Append(prod).ToKMatrix();
    }

    public CipherLWE EncryptBit(int bit)
    {
        var n = N;
        var qh = Q / 2 * One;
        var b0 = bit == 0 ? 0 : 1;
        var b = N.SeqLazy().Select(i => i < n - 1 ? qh.Zero : qh * b0).ToVec();
        var r = DistributionExt.DiceSample(M, [-1, 1]).Select(k => k * qh.One).ToVec();
        return new((PK * r).Select(e => e.Sum()).ToVec() + b);
    }

    public int DecryptBit(CipherLWE cipher)
    {
        var d = (SK * cipher.Vec).Sum();
        var d0 = d.K > Q / 2 ? Q - d.K : d.K;
        return d0 < Q / 4 ? 0 : 1;
    }

    public Vec<ZnInt64> Not(Vec<ZnInt64> a)
    {
        var n = N;
        var qh = Q / 2 * One;
        return a.Select((e, i) => i == n - 1 ? e + qh : e).ToVec();
    }
    
    public Vec<ZnInt64> Xor(Vec<ZnInt64> a, Vec<ZnInt64> b) => a + b;
    
    public Vec<ZnInt64> And(Vec<ZnInt64> a, Vec<ZnInt64> b)
    {
        var oq = a.KOne;
        var oq2 = EK[0].KOne;
        var q = Q;
        var cstar = a.Grid2D(b).Select(e => Signed(e.t1) * Signed(e.t2) * oq2).ToVec();
        return (cstar * EK).Select(e => (Signed(e.Sum()) * 2 / q) * oq).ToVec();
    }
    
    public Vec<ZnInt64> Nand(Vec<ZnInt64> a, Vec<ZnInt64> b) => Not(And(a, b));
    public string Params => $"LWE Params N:{N} M:{M} Q:{Q} Q^2:{Q * Q} Err:{Err}";

    public void Show()
    {
        Console.WriteLine(Params);
        Console.WriteLine($"SK: {SK}");
        PK.Println("PK");
        if (N <= 8)
            EK.Transpose().Println("EK");
        
        Console.WriteLine();
    }

    public static int Signed(ZnInt64 a) => a.K * 2 > a.P ? (int)a.K - a.P : (int)a.K;
}
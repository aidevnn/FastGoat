using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Integers;

namespace FastGoat.UserGroup.Lattice;

/// <summary>
/// Lattices and Homomorphic Encryption, Spring 2013
/// Instructors: Shai Halevi, Tal Malkin
/// 
/// LWE-based Homomorphic Encryption
/// 
/// April 12-16, 2013
/// Scribe: Kina Winoto, Clément Canonne
/// <see href="https://shaih.github.io/columbia6261/lecture12-LWE-HE.pdf"/>
/// </summary>
public struct LWE
{
    public static int S { get; set; } = 3;
    public int N { get; }
    public int M { get; }
    public long Q { get; }
    public ZnInt64 KOne { get; }
    public Vec<Vec<ZnInt64>> PK { get; }
    public Vec<ZnInt64> SK { get; }
    public Vec<Vec<ZnInt64>> EK { get; }

    public LWECipher Zero { get; }
    public LWECipher One { get; }

    public LWE(int n, long q, int m, bool noiseMode = true)
    {
        (N, M, Q) = (n, m, q);
        KOne = new ZnInt64(q, 1);

        var s1 = Regev.Ternary(n - 1, q).ToKMatrix();
        var e1 = Regev.DiscGauss(m, q, S).ToKMatrix(); // e1 ~ -1, 0, 1
        var A = Regev.Unif(m * (n - 1), q).ToKMatrix(n - 1);
        if (!noiseMode)
            e1 = e1.Zero;

        var a = s1 * A + e1;
        var pk = KMatrix<ZnInt64>.MergeSameCols(A, a);
        var sk = s1.Append(-s1.KOne).ToKMatrix();
        PK = pk.Rows.Select(r => r.ToVec()).ToVec();
        SK = sk.ToVec();

        long lq = q;
        var s1q2 = s1.Select(e => new ZnInt64(0, e.Signed)).ToKMatrix();
        var sks = SK.Select(e => e.Signed).ToArray();
        var sksk = sks.Grid2D(sks).Select(t => new ZnInt64(0, -t.t1 * t.t2 * q)).ToKMatrix(n * n);
        var B = Regev.Unif(n * n * (n - 1), q).Select(k => new ZnInt64(0,k.Signed * q)).ToKMatrix(n * n);

        var e2 = Regev.DiscGauss(n * n, lq * lq, S).Select(e => new ZnInt64(0, e.K)).ToKMatrix(n * n); // e2 ~ -1, 0, 1
        if (!noiseMode)
            e2 = e2.Zero;

        var b = B * s1q2.T + sksk + e2;
        var noise = (N * N).SeqLazy().Select(_ => BlankNoise(n, lq * lq, sks)).ToArray();
        var ek = KMatrix<ZnInt64>.MergeSameRows(B, b) + KMatrix<ZnInt64>.MergeSameCols(noise);
        EK = ek.Cols.Select(c => c.Select(e => new ZnInt64(lq * lq, e.K)).ToVec()).ToVec();

        Zero = new(N.SeqLazy().Select(_ => new ZnInt64(q, 0)).ToVec());
        One = LWECipher.Not(Zero);
    }

    static KMatrix<ZnInt64> BlankNoise(int n, long q, long[] sks)
    {
        var o = new ZnInt64(0, 1);
        var l = 10000.SeqLazy().Select(_ => IntExt.Rng.NextInt64(q / 10, q) * o).Where(k => k.K % q != 0).Take(n - 1)
            .ToArray();
        var prod = sks.Zip(l).Select(t => t.First * t.Second).Aggregate((a0, a1) => a0 + a1);
        return l.Append(prod).ToKMatrix();
    }

    public LWECipher EncryptBit(int bit)
    {
        var n = N;
        var qh = Q / 2 * KOne;
        var b0 = bit == 0 ? 0 : 1;
        var b = N.SeqLazy().Select(i => i < n - 1 ? qh.Zero : qh * b0).ToVec();
        var (m, q) = (M, Q);
        var r = Regev.DiscGauss(m, q, S);
        return new((PK * r).Select(e => e.Sum()).ToVec() + b);
    }

    public int DecryptBit(LWECipher lweCipher)
    {
        var d = (SK * lweCipher.Vec).Sum();
        return long.Abs(d.Signed) < Q / 4 ? 0 : 1;
    }

    public LWECipher BrutalRefresh(LWECipher lweCipher)
    {
        return EncryptBit(DecryptBit(lweCipher));
    }

    public int Error(LWECipher lweCipher)
    {
        var d = (SK * lweCipher.Vec).Sum();
        var m = long.Abs(d.Signed) < Q / 4 ? 0 : 1;
        return (int)(d - m * Q / 2).Signed;
    }

    public string Params => $"LWE Params N:{N} M:{M} Q:{Q}";

    public void Show()
    {
        Console.WriteLine(Params);
        Console.WriteLine($"SK: {SK}");
        // PK.Println("PK");
        // EK.Transpose().Println("EK");

        Console.WriteLine();
    }

    public static LWE Setup(int n, int k, bool noiseMode = true, bool pow2 = false)
    {
        var (q, m) = Regev.Setup(n, k, pow2);
        return new LWE(n, q, m, noiseMode);
    }

    public static LWE Setup(int n, bool noiseMode = true, bool pow2 = false) => Setup(n, 2, noiseMode, pow2);
}
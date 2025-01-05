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
    public int Q { get; }
    public ZnInt64 One { get; }
    public Vec<Vec<ZnInt64>> PK { get; }
    public Vec<ZnInt64> SK { get; }
    public Vec<Vec<ZnInt64>> EK { get; }

    public LWE(int n)
    {
        var alpha = 1.0 / (4 * double.Log2(n) * double.Log2(n) * double.Sqrt(n));
        var q = IntExt.Primes10000.First(q => alpha * q > 2 * double.Sqrt(n));
        var m = (int)(1.1 * n * double.Log2(q));
        (N, M, Q) = (n, m, q);
        One = new ZnInt64(q, 1);
        
        var s1 = Regev.Ternary(n - 1, q).ToKMatrix();
        var e1 = Regev.DiscGauss(m, q, S).ToKMatrix(); // e1 ~ -1, 0, 1
        var A = Regev.Unif(m * (n - 1), q).ToKMatrix(n - 1);
        var a = s1 * A + e1;
        var pk = KMatrix<ZnInt64>.MergeSameCols(A, a);
        var sk = s1.Append(-s1.KOne).ToKMatrix();
        PK = pk.Rows.Select(r => r.ToVec()).ToVec();
        SK = sk.ToVec();

        var o = new ZnInt64(q * q, 1);
        var s1q2 = s1.Select(e => e.Signed * o).ToKMatrix();
        var sks = SK.Select(e => e.Signed).ToArray();
        var sksk = sks.Grid2D(sks).Select(t => -t.t1 * t.t2 * o).ToKMatrix(n * n);
        var B = Regev.Unif(n * n * (n - 1), q).Select(k => k.K * q * o).ToKMatrix(n * n);

        var e2 = Regev.DiscGauss(n * n, o.P, S).ToKMatrix(n * n);  // e2 ~ -1, 0, 1
        var b = B * s1q2.T + sksk * q + e2;
        var noise = (N * N).SeqLazy().Select(_ => BlankNoise(n, o.P, sks)).ToArray();
        var ek = KMatrix<ZnInt64>.MergeSameRows(B, b) + KMatrix<ZnInt64>.MergeSameCols(noise);
        EK = ek.Cols.Select(c => c.ToVec()).ToVec();
    }

    static KMatrix<ZnInt64> BlankNoise(int n, int q, long[] sks)
    {
        var o = new ZnInt64(q, 1);
        var l = 1000.SeqLazy().Select(_ => IntExt.Rng.Next(q / 10, q) * o).Where(k => !k.IsZero()).Take(n - 1)
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
        var r = Regev.DiscGauss(M, Q, S).ToVec();
        return new((PK * r).Select(e => e.Sum()).ToVec() + b);
    }

    public int DecryptBit(CipherLWE cipher)
    {
        var d = (SK * cipher.Vec).Sum();
        return long.Abs(d.Signed) < Q / 4 ? 0 : 1;
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
}
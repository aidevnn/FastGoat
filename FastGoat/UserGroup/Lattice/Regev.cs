using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Integers;
using static FastGoat.Commons.IntExt;

namespace FastGoat.UserGroup.Lattice;

/// <summary>
/// On Lattices, Learning with Errors,
/// Random Linear Codes, and Cryptography
/// Oded Regev
/// January 9, 2024
/// https://arxiv.org/abs/2401.03703
/// </summary>
public class Regev
{
    public int N { get; }
    public int Q { get; }
    public int M { get; }
    public double A { get; }
    private Vec<ZnInt64> SK { get; }
    public Vec<ZnInt64> Err { get; }
    public (Vec<ZnInt64>[] A, Vec<ZnInt64> B) PK { get; }

    public Regev(int n)
    {
        // All conditions
        
        // 1) m = (1+ε)*n*log2(n) ε=1
        var m = (int)(2 * n * double.Log2(n));
        
        // 2) lim a(n)*log2(n)*√n = lim 1 / (4*log10(n)) = 0
        var a = 1.0 / (4 * double.Log10(n) * double.Log2(n) * double.Sqrt(n));
        
        // 3) a*q > 2*√n
        var q = Primes10000.First(q => a * q > 2 * double.Sqrt(n));
        (N, Q, M, A) = (n, q, m, a);
        
        // 4) Discrete gaussian sample with s = a*q
        var err = Err = DiscGauss(m, q, s: a * q);

        var sk = SK = Unif(n, q);
        var ai = m.SeqLazy().Select(_ => Unif(n, q)).ToArray();
        var b = ai.Zip(err).Select(e => (e.First * sk).Sum() + e.Second).ToVec();
        PK = (ai, b);
    }

    public RegevCipher EncryptBit(int m)
    {
        // 5) Set S ⊂ [0..M-1] 
        var r = DistributionExt.DiceSample(M, [true, false]).Zip(M.Range())
            .Where(e => e.First).Select(e => e.Second).ToArray();
        var m1 = m == 0 ? 0 : 1;
        var acc0 = (a: PK.A[0].Zero, b: new ZnInt64(Q, m1 * (Q / 2)));
        return r.Select(e => (a: PK.A[e], b: PK.B[e])).Aggregate(acc0, (acc, e) => (acc.a + e.a, acc.b + e.b));
    }

    public int DecryptBit(RegevCipher cipher)
    {
        // 6) b - <s,a> distance to 0 and to Q/2
        var d = cipher.B - (cipher.A * SK).Sum();
        return long.Abs(d.Signed) < Q / 4 ? 0 : 1;
    }

    public string Params =>
        $"Regev N:{N,-4} Q:{Q,-6} M:{M,-6} a:{A:F4} aq:{A * Q:F4} Q/M:{Q / M}";

    public void Show()
    {
        Console.WriteLine(Params);
        Console.WriteLine($"Private key:{SK}");
        // M.SeqLazy().Select(i => (PK.A[i], PK.B[i])).Println("Public Key");
        Console.WriteLine();
    }

    public static Vec<ZnInt64> Unif(int n, int q) =>
        DistributionExt.DiceSample(n, 0, q - 1).Select(i => new ZnInt64(q, i)).ToVec();

    public static Vec<ZnInt64> DiscGauss(int n, int q, double s)
    {
        var sigma = s / double.Sqrt(2 * double.Pi);
        return DistributionExt.DiscreteGaussianSample(n, sigma, tau: 1.0).Select(i => new ZnInt64(q, i)).ToVec();
    }
}

public readonly struct RegevCipher
{
    public Vec<ZnInt64> A { get; }
    public ZnInt64 B { get; }

    public RegevCipher(Vec<ZnInt64> a, ZnInt64 b)
    {
        (A, B) = (a, b);
    }

    public override string ToString() => $"[{A}, {B}]";

    public static implicit operator RegevCipher((Vec<ZnInt64> a, ZnInt64 b) e) => new(e.a, e.b);
}
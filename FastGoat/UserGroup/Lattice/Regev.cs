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
    public int P { get; }
    public int M { get; }
    public double A { get; }
    private Vec<ZnInt64> SK { get; }
    public Vec<ZnInt64> Err { get; }
    public (Vec<ZnInt64>[] A, Vec<ZnInt64> B) PK { get; }

    public Regev(int n)
    {
        // All conditions
        
        // 1) lim a(n)*log2(n)*√n = lim 1/log2(n) = 0
        var a = 1.0 / (double.Log2(n) * double.Log2(n) * double.Sqrt(n));
        
        // 2) a*p > 2*√n
        var p = Primes10000.First(p => a * p > 2 * double.Sqrt(n));
        
        // 3) m = (1+ε)*(n+1)*log2(p) ε=0.1
        var m = (int)(1.1 * (n + 1) * double.Log2(p));
        (N, P, M, A) = (n, p, m, a);
        
        // 4) Discrete gaussian sample with s = a*p
        var err = Err = DiscGauss(m, p, s: a * p);

        var sk = SK = Unif(n, p);
        var ai = m.SeqLazy().Select(_ => Unif(n, p)).ToArray();
        var b = ai.Zip(err).Select(e => (e.First * sk).Sum() + e.Second).ToVec();
        PK = (ai, b);
    }

    public Regev(int n, int m, int p, double sigma)
    {
        var s = double.Sqrt(2 * double.Pi) * sigma;
        (N, P, M, A) = (n, p, m, s / p);
        
        var err = Err = DiscGauss(m, p, s);
        var sk = SK = Unif(n, p);
        var ai = m.SeqLazy().Select(_ => Unif(n, p)).ToArray();
        var b = ai.Zip(err).Select(e => (e.First * sk).Sum() + e.Second).ToVec();
        PK = (ai, b);
    }

    public RegevCipher EncryptBit(int m)
    {
        var _m = m == 0 ? 0 : 1;
        var acc0 = (a: PK.A[0].Zero, b: new ZnInt64(P, _m * (P / 2)));
        
        // 5) Set S ⊂ [0..M-1] 
        return DistributionExt.DiceSample(M, [true, false]).Zip(M.Range())
            .Where(e => e.First)
            .Select(e => (a: PK.A[e.Second], b: PK.B[e.Second]))
            .Aggregate(acc0, (acc, e) => (acc.a + e.a, acc.b + e.b));
    }

    public int DecryptBit(RegevCipher cipher)
    {
        // 6) b - <s,a> distance to 0 and to P/2
        var d = cipher.B - (cipher.A * SK).Sum();
        return long.Abs(d.Signed) < P / 4 ? 0 : 1;
    }

    public string Params => $"Regev N:{N,-4} P:{P,-6} M:{M,-6} A:{A:F4} A*P:{A * P:F4} P/M:{P / (1.0 * M):F4}";

    public void Show()
    {
        Console.WriteLine(Params);
        Console.WriteLine($"Private key:{SK}");
        // M.SeqLazy().Select(i => (PK.A[i], PK.B[i])).Println("Public Key");
        Console.WriteLine();
    }

    public static Vec<ZnInt64> Unif(int n, long q) =>
        DistributionExt.DiceSample(n, 0, q - 1).Select(i => new ZnInt64(q, i)).ToVec();

    public static Vec<ZnInt64> Ternary(int n, long q) =>
        DistributionExt.DiceSample(n, [-1, 0, 1]).Select(i => new ZnInt64(q, i)).ToVec();

    public static Vec<ZnInt64> DiscGauss(int n, long q, double s, double tau= 1.0)
    {
        var sigma = s / double.Sqrt(2 * double.Pi);
        return DistributionExt.DiscreteGaussianSample(n, sigma, tau).Select(i => new ZnInt64(q, i)).ToVec();
    }

    public static (int q, int m) Setup(int n, int k, bool pow2 = false)
    {
        var l = 1 << k;
        var alpha = 1.0 / (l * double.Log2(n) * double.Log2(n) * double.Sqrt(n));
        var seq = pow2 ? 16.SeqLazy(1).Select(i => 1 << i) : IntExt.Primes10000;
        var q = seq.First(q0 => alpha * q0 > 2 * double.Sqrt(n));
        var m = (int)(1.1 * (n + 1) * double.Log2(q));
        return (q, m);
    }
}
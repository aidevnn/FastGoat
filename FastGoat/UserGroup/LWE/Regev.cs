using System.Numerics;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Integers;
using static FastGoat.Commons.IntExt;

namespace FastGoat.UserGroup.LWE;

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
    private Vec<ZnBigInt> SK { get; }
    public Vec<ZnBigInt> Err { get; }
    public RegevCipher[] PK { get; }

    public Regev(int n, bool unifSK = true)
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

        var sk = SK = unifSK ? Unif(n, p) : TernarySK(n, p);
        PK = err.Select(e => (a:Unif(n, p), b:e))
            .Select(e => new RegevCipher(e.a, e.b + e.a.InnerProd(sk)))
            .ToArray();
    }

    public RegevCipher EncryptBit(int bit)
    {
        var m = bit == 0 ? 0 : P / 2;
        
        // 5) Set S ⊂ [0..M-1] 
        return PK.Where(_ => Rng.Next(2) == 1).Aggregate((ci, cj) => ci + cj) + m;
    }

    public int DecryptBit(RegevCipher cipher)
    {
        // 6) b - <s,a> distance to 0 and to P/2
        var d = cipher.B - cipher.A.InnerProd(SK);
        return BigInteger.Abs(d.K) <= P / 4 ? 0 : 1;
    }

    public RegevCipher[] Encrypt(int[] bits) => bits.Select(EncryptBit).ToArray();
    public int[] Decrypt(RegevCipher[] ciphers) => ciphers.Select(DecryptBit).ToArray();

    public string Params => $"Regev N:{N,-4} P:{P,-6} M:{M,-6} A:{A:F3} A*P:{A * P:F4}";

    public string ExportParams
    {
        get
        {
            var d = (int)((N - double.Sqrt(N)) / 2);
            var sigma = A * P / double.Sqrt(2 * double.Pi);
            return $"LWEParameters(n={N}, q={P}, Xs=U({P}), Xe=D({sigma:f4}), m=+oo, tag='fg{N}')";
        }
    }
    public void Show()
    {
        Console.WriteLine(Params);
        Console.WriteLine($"Private key:{SK}");
        // PK.Println("Public Key");
        Console.WriteLine();
    }

    public ZnBigInt Errors(RegevCipher cipher)
    {
        var m = DecryptBit(cipher);
        return cipher.B - (cipher.A * SK).Sum() - m * P / 2;
    }

    public ZnBigInt[] Errors(RegevCipher[] ciphers) => ciphers.Select(Errors).ToArray();

    public void Deconstruct(out int n, out int p, out RegevCipher[] pk, out Vec<ZnBigInt> sk, out Vec<ZnBigInt> err)
    {
        (n, p, pk, sk, err) = (N, P, PK, SK, Err);
    }

    public static Vec<ZnBigInt> Unif(int n, long q) =>
        DistributionExt.DiceSample(n, 0, q - 1).Select(i => new ZnBigInt(q, i)).ToVec();

    public static Vec<ZnBigInt> Ternary(int n, long q) =>
        DistributionExt.DiceSample(n, [-1, 0, 1]).Select(i => new ZnBigInt(q, i)).ToVec();

    public static Vec<ZnBigInt> DiscGauss(int n, long q, double s, double tau= 1.0)
    {
        var sigma = s / double.Sqrt(2 * double.Pi);
        return DistributionExt.DiscreteGaussianSample(n, sigma, tau).Select(i => new ZnBigInt(q, i)).ToVec();
    }

    public static Vec<ZnBigInt> TernarySK(int n, long q)
    {
        var o = new ZnBigInt(q, 1);
        var arr = DistributionExt.DiceSample(n, [-o, o]).ToArray();
        var nb = (int)double.Sqrt(n);
        var zeros = DistributionExt.DiceSample(nb, (n - 1).Range()).ToArray();
        foreach (var i in zeros) arr[i] *= 0;
        return arr.ToVec();
    }

    public static (Regev regev, RLWE rlwe, RLWECipher swk, RLWECipher[] exsk) SetupRLWE(int n)
    {
        var regev = new Regev(n, unifSK: false);
        var rlwe = new RLWE(2 * n, regev.P);
        var sk2 = RLWE.IntVecToRq(regev.SK);
        var swk = RLWE.SWKBGV(rlwe.PM, sk2, rlwe.SK, rlwe.T, rlwe.PK.Q, rlwe.SP);
        var exsk = RLWE.EXSK(sk2, rlwe.PK);
        return (regev, rlwe, swk, exsk);
    }
}
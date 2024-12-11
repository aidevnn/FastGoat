using FastGoat.Commons;
using FastGoat.UserGroup.Integers;

namespace FastGoat.UserGroup.Lattice;

/// <summary>
/// On Lattices, Learning with Errors,
/// Random Linear Codes, and Cryptography
/// Oded Regev
/// January 9, 2024
/// https://arxiv.org/abs/2401.03703
/// 
/// #5 Public Key Cryptosystem
/// </summary>
public class LWE
{
    public CypherLWE[] PK { get; }
    private ZnInt64[] SK { get; }
    private ZnInt64 Zero { get; }
    public int P { get; }
    public int M { get; }
    public int N { get; }
    public int Err { get; }
    public double A { get; }
    public string Fmt { get; }

    static ZnInt64 Sum(IEnumerable<ZnInt64> A) => A.Aggregate((a0, a1) => a0 + a1);

    static ZnInt64 InnerProd(IEnumerable<ZnInt64> A, IEnumerable<ZnInt64> B)
        => Sum(A.Zip(B).Select(e => e.First * e.Second));

    public LWE(int n)
    {
        if (n < 2 || n > 20)
            throw new($"N = {n} must be > 1 and < 21");

        M = (int)(2 * n * double.Log2(n));
        var err = (int)(M - M / (2 * double.Sqrt(n)));
        var p = IntExt.Primes10000.First(p => p > 4 * n * err);
        var a = 1 / (double.Log2(n) * double.Sqrt(n));
        (P, N, A, Err) = (p, n, a, err % n);
        Zero = ZnInt64.ZnZero(P);

        SK = DistributionExt.DiceSample(N, 1, p - 1).Select(e => new ZnInt64(P, e)).ToArray();
        PK = M.Range().Select(_ => DistributionExt.DiceSample(N, 1, p - 1).Select(e => new ZnInt64(P, e)).ToArray())
            .Select(Ai => new CypherLWE(Ai, InnerProd(Ai, SK) + DistributionExt.BinomialEquiProb(Err)))
            .ToArray();

        var nbDigits = $"{P + 1}".Length;
        Fmt = $"{{0,{nbDigits}}}";
    }

    public LWE(int m, int n)
    {
        if (n < 2 || n > 20)
            throw new($"N = {n} must be > 1 and < 21");
        
        M = m;
        var err = (int)(M - M / (2 * double.Sqrt(n)));
        var p = IntExt.Primes10000.First(p => p > 4 * n * err);
        var a = 1 / (double.Log2(n) * double.Sqrt(n));
        (P, N, A, Err) = (p, n, a, err % n);
        Zero = ZnInt64.ZnZero(P);

        SK = DistributionExt.DiceSample(N, 1, p - 1).Select(e => new ZnInt64(P, e)).ToArray();
        PK = M.Range().Select(_ => DistributionExt.DiceSample(N, 1, p - 1).Select(e => new ZnInt64(P, e)).ToArray())
            .Select(Ai => new CypherLWE(Ai, InnerProd(Ai, SK) + DistributionExt.BinomialEquiProb(Err)))
            .ToArray();

        var nbDigits = $"{P + 1}".Length;
        Fmt = $"{{0,{nbDigits}}}";
    }

    public CypherLWE EncryptBit(int b)
    {
        var m = PK.Length;
        var n = PK[0].Ai.Length;
        var r = IntExt.Rng.Next(2, m);
        var set = r.Range().Select(_ => PK[IntExt.Rng.Next(m)]).ToArray();
        var sumBi = Sum(set.Select(f => f.B));
        var sumAi = n.Range().Select(i => Sum(set.Select(f => f.Ai[i]))).ToArray();
        if (b != 0)
            sumBi += P / 2;

        return new(sumAi, sumBi);
    }

    public int DecryptBit(CypherLWE cypherLwe)
    {
        var c = cypherLwe.B - InnerProd(cypherLwe.Ai, SK);
        return c.K < P / 2 ? 0 : 1;
    }

    public void Show(bool showKeys = true)
    {
        Console.WriteLine($"P:{P} N:{N} M:{M} Err:{Err} A:{A:G6}");
        if (showKeys)
        {
            Console.WriteLine("Private Key");
            Console.WriteLine($"[{SK.Glue(", ")}]");
            PK.Println(e => $"[[{e.Ai.Glue(", ")}], {e.B}]", "Public Key");
        }

        Console.WriteLine();
    }
}
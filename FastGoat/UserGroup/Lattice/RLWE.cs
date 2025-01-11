using FastGoat.Commons;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Integers;
using Rq = FastGoat.Structures.VecSpace.KPoly<FastGoat.UserGroup.Integers.Rational>;

namespace FastGoat.UserGroup.Lattice;

public class RLWE
{
    public int N { get; }
    public int Q { get; }
    public double Sigma { get; }

    public RLWE(int n)
    {
        if (!int.IsPow2(n))
            throw new($"N = {n} must be 2^k");

        var q = IntExt.Primes10000.First(q => q % (2 * n) == 1);
        (N, Q, Sigma) = (n, q, 3.0);
    }
    
    public RLWE(int n, int q, double sigma = 3.0)
    {
        if (!int.IsPow2(n))
            throw new($"N = {n} must be 2^k");

        if (q % (2 * n) != 1)
            throw new($"Q = {q} must be = 1 mod 2n");
        
        (N, Q, Sigma) = (n, q, sigma);
    }

    public Rq Encode(int[] seq)
    {
        if (seq.Length != N)
            throw new();

        var qh = Q / 2;
        return seq.Select(e => qh * IntExt.AmodP(e, 2)).ToKPoly(Rational.KOne());
    }

    public Rq Mask() => Encode(Enumerable.Repeat(1, N).ToArray());

    public int[] Decode(Rq a)
    {
        return N.Range().Select(i => a[i]).Select(k => k * 2 > Q ? Q - k : k)
            .Select(i => i.Num * 4 < Q ? 0 : 1).ToArray();
    }

    public void Show(bool showKeys = true)
    {
        var k = int.Log2(N);
        Console.WriteLine($"RLWE N = {N} = 2^{k} Q = {Q} Sigma = {Sigma}");
        Console.WriteLine();
    }
}
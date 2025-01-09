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
    public RLWE(int n, int q, double sigma = 1.0)
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

        return seq.Take(N).Select(e => IntExt.AmodP(e, 2)).ToKPoly(Rational.KOne());
    }

    public int[] Decode(Rq a)
    {
        return N.Range().Reverse().Select(i => a[i]).Select(k => k * 2 > Q ? Q - k : k)
            .Select(i => i.Num * 4 < Q ? 0 : 1).ToArray();
    }

    public void Show(bool showKeys = true)
    {
        var k = int.Log2(N);
        Console.WriteLine($"RLWE N = {N} = 2^{k} Q = {Q} Sigma = {Sigma}");
        Console.WriteLine();
    }
}
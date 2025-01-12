using FastGoat.Commons;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Integers;

namespace FastGoat.UserGroup.Lattice;

public class RLWE
{
    public int N { get; }
    public double Sigma { get; }

    public KPoly<Rational> PM { get; }

    public KPoly<Rational> SK { get; }

    public Rational T { get; }

    public Rational Q { get; }

    public RLWECipher PK { get; }

    public RLWECipher RLK { get; }

    public RLWE(int n)
    {
        if (!int.IsPow2(n))
            throw new($"N = {n} must be 2^k");

        var t0 = IntExt.Primes10000.First(t0 => t0 % (2 * n) == 1);
        var t1 = IntExt.Primes10000.First(t1 => t1 > t0);
        (N, Sigma) = (n, 3.0);
        (PM, SK, T, Q, PK, RLK) = FHE.KeyGenBGV(N, t0, t0 * t1);
    }

    public RLWECipher[] Encrypt(int[] seq)
    {
        var qh = (T * PM.One / 2).FloorPoly();
        return seq.Select(e => FHE.EncryptBGV(qh * IntExt.AmodP(e, 2), PM, T, Q, PK)).ToArray();
    }

    public int[] Decrypt(RLWECipher[] cipher)
    {
        return cipher.Select(e => (FHE.DecryptBGV(e, PM, SK, T) * 2 / T).RoundPoly().IsZero() ? 0 : 1).ToArray();
    }

    public void Show()
    {
        var k = int.Log2(N);
        Console.WriteLine($"RLWE N = {N} = 2^{k}");
        FHE.Show(PM, SK, T, Q, PK, RLK);
        Console.WriteLine($"RLWE N = {N} = 2^{k} t = {T} q = {Q}");

        Console.WriteLine("Private Key");
        Console.WriteLine(SK);
        PK.Show("Public Key");
        RLK.Show("Relinearisation Key");
        Console.WriteLine();
    }

    public RLWECipher[] NOT(RLWECipher[] cipher)
    {
        var t0 = (int)T.Num;
        var th = t0 / 2;
        var thi = t0 - th;
        var z = IntExt.InvModPbez(thi, t0);
        var f = z * th % t0;
        return cipher.Select(c => (f * (c + thi)).CoefsMod(Q)).ToArray();
    }

    public RLWECipher[] AND(RLWECipher[] cipher1, RLWECipher[] cipher2)
    {
        var t0 = (int)T.Num;
        var th = t0 / 2;
        var err = th * th % t0;
        var z = IntExt.InvModPbez(err, t0);
        var f = z * th % t0;
        return cipher1.Zip(cipher2).Select(e => f * FHE.MulRelinBGV(e.First, e.Second, PM, Q, RLK)).ToArray();
    }

    public RLWECipher[] NAND(RLWECipher[] cipher1, RLWECipher[] cipher2) => NOT(AND(cipher1, cipher2));
    public RLWECipher[] NOR(RLWECipher[] cipher1, RLWECipher[] cipher2) => AND(NOT(cipher1), NOT(cipher2));
    public RLWECipher[] OR(RLWECipher[] cipher1, RLWECipher[] cipher2) => NOT(NOR(cipher1, cipher2));
    public RLWECipher[] XOR(RLWECipher[] cipher1, RLWECipher[] cipher2)
    {
        var c1_and_nc2 = AND(cipher1, NOT(cipher2));
        var nc1_and_c2 = AND(NOT(cipher1), cipher2);
        return OR(c1_and_nc2, nc1_and_c2);
    }

    public RLWECipher[] OP(string name, RLWECipher[] cipher1, RLWECipher[] cipher2)
    {
        name = name.ToLower();
        if (name == "and")
            return AND(cipher1, cipher2);
        if (name == "nand")
            return NAND(cipher1, cipher2);
        if (name == "nor")
            return NOR(cipher1, cipher2);
        if (name == "or")
            return OR(cipher1, cipher2);
        if (name == "xor")
            return XOR(cipher1, cipher2);
        else
            throw new($"op:{name} not found");
    }
    
    public RLWECipher[] ADD(RLWECipher[] c1, RLWECipher[] c2)
    {
        var (g1, g2) = (c1.ToArray(), c2.ToArray());
        var z = c1[0].Zero;
        for (int i = 0; i < c1.Length; i++)
        {
            var xor_g1g2 = XOR(g1, g2);
            var and_g1g2 = AND(g1, g2).SkipLast(1).Prepend(z).ToArray();
            (g1, g2) = (xor_g1g2, and_g1g2);
        }

        return g1;
    }

}
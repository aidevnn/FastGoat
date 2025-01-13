using FastGoat.Commons;
using FastGoat.UserGroup.Integers;
using Rq = FastGoat.Structures.VecSpace.KPoly<FastGoat.UserGroup.Integers.Rational>;

namespace FastGoat.UserGroup.Lattice;

public class RLWE
{
    public int n { get; }
    public int N { get; }
    public double Sigma { get; }

    public Rq PM { get; }

    public Rq SK { get; }

    public Rational T { get; }

    public Rational Q { get; }

    public RLWECipher PK { get; }

    public RLWECipher RLK { get; }

    public Rational Thalf { get; }
    public Rational InvThalf { get; }
    public RLWECipher Zero => PK.Zero;

    public RLWE(int N)
    {
        if (!int.IsPow2(N))
            throw new($"N = {N} must be 2^k");

        var t0 = IntExt.Primes10000.First(t0 => t0 % (2 * N) == 1);
        var t1 = IntExt.Primes10000.First(t1 => t1 > t0);
        (this.N, n, Sigma) = (N, IntExt.Phi(N), 3.0);
        
        (PM, SK, T, Q, PK, RLK) = FHE.KeyGenBGV(n, t0, t0 * t1);
        Thalf = new(t0 / 2);
        InvThalf = new Rational(IntExt.InvModPbez(t0 / 2, t0));
    }

    public RLWECipher EncryptBit(int m)
    {
        return FHE.EncryptBGV(Thalf * PM.One * IntExt.AmodP(m, 2), PM, T, Q, PK);
    }

    public RLWECipher[] Encrypt(int[] seq) => seq.Select(EncryptBit).ToArray();

    public int DecryptBit(RLWECipher cipher)
    {
        return (FHE.DecryptBGV(cipher, PM, SK, T) * InvThalf).CoefsMod(T).IsZero() ? 0 : 1;
    }

    public int[] Decrypt(RLWECipher[] ciphers) => ciphers.Select(DecryptBit).ToArray();

    public string Params => $"RLWE N={N}=2^{int.Log2(N)}, Φ(N)={n} PM={PM} t={T} q={Q}={T}*{Q/T}";

    public void Deconstruct(out int _n, out Rq pm, out Rq sk, out Rational t, out Rational q, out RLWECipher pk, out RLWECipher rlk)
    {
        (_n, pm, sk, t, q, pk, rlk) = (n, PM, SK, T, Q, PK, RLK);
    }

    public void Show()
    {
        Console.WriteLine(Params);
        Console.WriteLine("Private Key");
        Console.WriteLine(SK);
        PK.Show("Public Key");
        RLK.Show("Relinearisation Key");
        Console.WriteLine();
    }

    public RLWECipher NOT(RLWECipher cipher) => (Thalf - cipher).CoefsMod(Q);

    public RLWECipher[] NOT(RLWECipher[] ciphers) => ciphers.Select(NOT).ToArray();

    public RLWECipher AND(RLWECipher cipher1, RLWECipher cipher2)
    {
        return (InvThalf * FHE.MulRelinBGV(cipher1, cipher2, PM, Q, RLK)).CoefsMod(Q);
    }

    public RLWECipher NAND(RLWECipher cipher1, RLWECipher cipher2) => NOT(AND(cipher1, cipher2));
    public RLWECipher NOR(RLWECipher cipher1, RLWECipher cipher2) => AND(NOT(cipher1), NOT(cipher2));
    public RLWECipher OR(RLWECipher cipher1, RLWECipher cipher2) => NOT(NOR(cipher1, cipher2));
    public RLWECipher XOR(RLWECipher cipher1, RLWECipher cipher2) =>
        OR(AND(cipher1, NOT(cipher2)), AND(NOT(cipher1), cipher2));

    public RLWECipher OP(string name, RLWECipher cipher1, RLWECipher cipher2)
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
    
    public RLWECipher[] AND(RLWECipher[] ciphers1, RLWECipher[] ciphers2)
    {
        return ciphers1.Zip(ciphers2).Select(e => AND(e.First, e.Second)).ToArray();
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
        return cipher1.Zip(cipher2).Select(e => OP(name, e.First, e.Second)).ToArray();
    }

    public RLWECipher[] ADD(RLWECipher[] c1, RLWECipher[] c2)
    {
        var (g1, g2) = (c1.ToArray(), c2.ToArray());
        for (int i = 0; i < c1.Length; i++)
        {
            var xor_g1g2 = XOR(g1, g2);
            var and_g1g2 = AND(g1, g2).SkipLast(1).Prepend(Zero).ToArray();
            (g1, g2) = (xor_g1g2, and_g1g2);
        }

        return g1;
    }
    
    public RLWECipher[] MULT(RLWECipher[] ciphers1, RLWECipher[] ciphers2)
    {
        var list = new List<RLWECipher[]>();
        var length = ciphers1.Length + ciphers2.Length;
        foreach (var (cipher, pos) in ciphers1.Select((c,i)=>(c,i)))
        {
            var prepend = Enumerable.Repeat(Zero, pos).ToArray();
            var append = Enumerable.Repeat(Zero, length - pos).ToArray();
            var tmpMult = ciphers2.Select(c => AND(c, cipher));
            list.Add(prepend.Concat(tmpMult).Concat(append).ToArray());
        }
    
        var sum = Enumerable.Repeat(Zero, length).ToArray();
        foreach (var ciphers in list)
            sum = ADD(sum, ciphers);

        return sum;
    }
}
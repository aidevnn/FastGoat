using System.Numerics;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Integers;
using Rq = FastGoat.Structures.VecSpace.KPoly<FastGoat.UserGroup.Integers.Rational>;

namespace FastGoat.UserGroup.LWE;

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
        
        (PM, SK, T, Q, PK, RLK) = KeyGenBGV(n, t0, t0 * t1);
        Thalf = new(t0 / 2);
        InvThalf = new Rational(IntExt.InvModPbez(t0 / 2, t0));
    }

    public RLWE(int N, int t, Vec<ZnInt64> sk)
    {
        if (!int.IsPow2(N))
            throw new($"N = {N} must be 2^k");

        var t1 = IntExt.Primes10000.First(t1 => t1 > t);
        (this.N, n, Sigma) = (N, IntExt.Phi(N), 3.0);
        
        (PM, SK, T, Q, PK, RLK) = KeyGenBGV(n, t, t * t1, sk);
        Thalf = new(t / 2);
        InvThalf = new Rational(IntExt.InvModPbez(t / 2, t));
    }

    public RLWECipher EncryptBit(int m)
    {
        return EncryptBGV(Thalf * PM.One * IntExt.AmodP(m, 2), PM, T, Q, PK);
    }

    public RLWECipher[] Encrypt(int[] seq) => seq.Select(EncryptBit).ToArray();

    public int DecryptBit(RLWECipher cipher)
    {
        var d = DecryptBGV(cipher, PM, SK, T)[0].Signed(T).Absolute;
        return d * 4 < T ? 0 : 1;
    }

    public int[] Decrypt(RLWECipher[] ciphers) => ciphers.Select(DecryptBit).ToArray();

    public Rational Errors(RLWECipher cipher)
    {
        var m = DecryptBit(cipher);
        return (cipher.A - cipher.B * SK).ResMod(PM, T)[0] - m * Thalf;
    }

    public Rational[] Errors(RLWECipher[] ciphers) => ciphers.Select(Errors).ToArray();

    public RLWECipher FromRegevCipher(RegevCipher cipher)
    {
        var b = new Rational(cipher.B.K) * PM.One;
        var a = ExtractVec(cipher.A).CoefsMod(T);
        return new(b, a, PM, Q);
    }

    public RLWECipher[] FromRegevCipher(RegevCipher[] ciphers) => ciphers.Select(FromRegevCipher).ToArray();

    public RegevCipher ToRegevCipher(RLWECipher cipher)
    {
        var a = ExtractArr(cipher.B, n).Select(e => new ZnInt64((long)T.Num, (long)e.Num)).ToVec();
        var b = new ZnInt64((long)T.Num, (long)cipher.A[0].Num);
        return new(a, b);
    }

    public RegevCipher[] ToRegevCipher(RLWECipher[] ciphers) => ciphers.Select(ToRegevCipher).ToArray();

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
        return (InvThalf * MulRelinBGV(cipher1, cipher2, PM, Q, RLK)).CoefsMod(Q);
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
    
    public static bool NoiseMode { get; set; } = true;
    public static void NoiseOn() => NoiseMode = true;
    public static void NoiseOff() => NoiseMode = false;
    
    public static Rq GenDiscrGauss(int n, double s = 3.2)
    {
        return DistributionExt.DiscreteGaussianSample(n, s).ToKPoly(Rational.KZero());
    }

    public static Rq GenTernary(int n)
    {
        return DistributionExt.DiceSample(n, [-1, 0, 1]).ToKPoly(Rational.KZero());
    }

    public static Rq GenUnif(int n, int q)
    {
        return DistributionExt.DiceSample(n, 0, q - 1).Select(e => new Rational(e)).ToKPoly();
    }

    public static Rq GenUnif(int n, BigInteger q)
    {
        return DistributionExt.DiceSampleBigInt(n, 0, q - 1).Select(e => new Rational(e)).ToKPoly();
    }

    public static Rq GenUnif(int n, Rational q) => GenUnif(n, q.Num);

    public static (Rq pm, Rq sk, Rational t, Rational q, RLWECipher pk, RLWECipher rlk) KeyGenBGV(int n, int t, int q)
    {
        var pm = FG.QPoly().Pow(n) + 1;
        var sk = 10000.SeqLazy().Select(_ => GenTernary(n))
            .First(s => !s[n - 1].IsZero() && s.Coefs.Count(e => e.IsZero()) <= n / 4);
        
        var (T, Q) = (new Rational(t), new Rational(q));
        
        var epk = GenDiscrGauss(n);
        var c1pk = GenUnif(n, q);
        var c0pk = (t * epk + c1pk * sk).ResMod(pm, Q);
        var pk = new RLWECipher(c0pk, c1pk, pm, Q);
        
        var erlk = GenDiscrGauss(n);
        var c1rlk = GenUnif(n, q);
        var c0rlk = (t * erlk + c1rlk * sk + sk.Pow(2)).ResMod(pm, Q);
        var rlk = new RLWECipher(c0rlk, c1rlk, pm, Q);
        
        return (pm, sk, T, Q, pk, rlk);
    }


    public static (Rq pm, Rq sk, Rational t, Rational q, RLWECipher pk, RLWECipher rlk) KeyGenBGV(int n, int t, int q, Vec<ZnInt64> _sk)
    {
        var pm = FG.QPoly().Pow(n) + 1;
        var sk = _sk.Select(e => new Rational(e.Signed)).ToKPoly();
        
        var (T, Q) = (new Rational(t), new Rational(q));
        
        var epk = GenDiscrGauss(n);
        var c1pk = GenUnif(n, q);
        var c0pk = (t * epk + c1pk * sk).ResMod(pm, Q);
        var pk = new RLWECipher(c0pk, c1pk, pm, Q);
        
        var erlk = GenDiscrGauss(n);
        var c1rlk = GenUnif(n, q);
        var c0rlk = (t * erlk + c1rlk * sk + sk.Pow(2)).ResMod(pm, Q);
        var rlk = new RLWECipher(c0rlk, c1rlk, pm, Q);
        
        return (pm, sk, T, Q, pk, rlk);
    }

    public static RLWECipher EncryptBGV(Rq m, Rq pm, Rational t, Rational q, RLWECipher pk, bool noise = true)
    {
        var n = pm.Degree;
        var ea = GenDiscrGauss(n);
        var eb = GenDiscrGauss(n);
        var u = GenTernary(n);
        if (!NoiseMode || !noise)
        {
            ea = eb = eb.Zero;
            u = u.One;
        }

        var a = (u * pk.A + m + t * ea).ResMod(pm, q);
        var b = (u * pk.B + t * eb).ResMod(pm, q);
        return (a, b, pm, q);
    }

    public static Rq DecryptBGV(RLWECipher cipher, Rq pm, Rq sk, Rational t)
    {
        return (cipher.A - sk * cipher.B).ResMod(pm, t);
    }

    public static RLWECipher MulRelinBGV(RLWECipher cipher0, RLWECipher cipher1, Rq pm, Rational q, RLWECipher rlk)
    {
        var d0 = (cipher0.A * cipher1.A).ResMod(pm, q);
        var d1 = (cipher0.A * cipher1.B + cipher0.B * cipher1.A).ResMod(pm, q);
        var d2 = (cipher0.B * cipher1.B).ResMod(pm, q);

        var a = (d0 + d2 * rlk.A).ResMod(pm, q);
        var b = (d1 + d2 * rlk.B).ResMod(pm, q);
        return (a, b, pm, q);
    }
    
    public static Rational[] ExtractArr(Rq poly, int n, int i = 0)
    {
        return n.SeqLazy(start: i, step: -1).Select(j => j >= 0 ? poly[j] : -poly[n + j]).ToArray();
    }

    public static Rq ExtractVec(Vec<ZnInt64> v, int i = 0)
    {
        var x = FG.QPoly();
        var n = v.Length;
        return v.Select(e => new Rational(e.Signed))
            .Select((e, j) => i - j >= 0 ? e * x.Pow(i - j) : -e * x.Pow(n + i - j))
            .Aggregate((xi, xj) => xi + xj);
    }
}
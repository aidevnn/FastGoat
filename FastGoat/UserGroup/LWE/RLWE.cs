using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Integers;
using Rq = FastGoat.Structures.VecSpace.KPoly<FastGoat.UserGroup.Integers.Rational>;

namespace FastGoat.UserGroup.LWE;

public partial class RLWE
{
    public int n { get; }

    public int N { get; }

    public int Level { get; }

    public Rq PM { get; }

    public Rq SK { get; }

    public Rational T { get; }

    public Rational[] Primes { get; }

    public Rational SP { get; }

    public RLWECipher PK { get; }

    public RLWECipher RLK => RLKS[PK.Q].rlk;

    public Dictionary<Rational, (Rational nextMod, RLWECipher rlk)> RLKS { get; }

    public Rational Q => Primes[0];


    public Rational Thalf => (T / 2).Trunc;

    public Rational InvThalf => "-2";

    public bool BootstrappingMode { get; }

    public RLWE(int N)
    {
        (this.N, n, Level) = (N, IntExt.Phi(N), 1);
        (PM, SK, T, Primes, SP, PK, RLKS) = SetupBGV(N, Level);
    }

    public RLWE(int N, int level, bool bootstrappingMode)
    {
        (this.N, n, Level) = (N, IntExt.Phi(N), level);
        (PM, SK, T, Primes, SP, PK, RLKS) = SetupBGV(N, Level);
        BootstrappingMode = bootstrappingMode;
    }

    public RLWE(int N, int t)
    {
        (this.N, n, Level) = (N, IntExt.Phi(N), 1);
        (PM, SK, T, Primes, SP, PK, RLKS) = SetupBGV(N, t, Level);
    }

    public RLWECipher EncryptBit(int bit)
    {
        return EncryptBGV(Thalf * PM.One * IntExt.AmodP(bit, 2), PK);
    }

    public RLWECipher[] Encrypt(int[] bits) => bits.Select(EncryptBit).ToArray();

    public int DecryptBit(RLWECipher cipher)
    {
        var d = DecryptBGV(cipher, SK)[0].Signed(T).Absolute;
        return d <= (T / 4).Trunc ? 0 : 1;
    }

    public int[] Decrypt(RLWECipher[] ciphers) => ciphers.Select(DecryptBit).ToArray();

    public Rq Errors(RLWECipher cipher)
    {
        var m = DecryptBit(cipher);
        return (cipher.A - cipher.B * SK).ResModSigned(PM, cipher.Q) - m * Thalf;
    }

    public Rq[] Errors(RLWECipher[] ciphers) => ciphers.Select(Errors).ToArray();

    public string Params
    {
        get
        {
            var cycloPoly = $"N={N}=2^{int.Log2(N)}, Φ(N)={n}, PM={PM}";
            var primes = $"t={T}, q={Q}, sp={SP}, level={Level} primes=[{Primes.Glue(", ")}]";
            var sigma_gadget = $"Sigma={Sigma(n, T):f3}, GadgetBase = {GadgetBase(T)}";
            return $"RLWE {cycloPoly}, {primes}, {sigma_gadget}";
        }
    }

    public string ExportParams
    {
        get
        {
            var d = (int)((n - double.Sqrt(n)) / 2);
            return $"LWEParameters(n={n}, q={T}, Xs=T({d}, n={n}), Xe=D({Sigma(n, T):f4}), m=+oo, tag='fg{n}')";
        }
    }

    public void Deconstruct(out int _n, out Rq pm, out Rq sk, out Rational t, out Rational q, out RLWECipher pk,
        out RLWECipher rlk)
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

    public RLWECipher Zero => PK.Zero;

    private ((Vec<RLWECipher> csm, Vec<RLWECipher> cm) plus, (Vec<RLWECipher> csm, Vec<RLWECipher> cm) minus)[]
        brk = [];

    public ((Vec<RLWECipher> csm, Vec<RLWECipher> cm) plus, (Vec<RLWECipher> csm, Vec<RLWECipher> cm) minus)[]
        BlindRotateKeys
    {
        get
        {
            if (brk.Length != 0)
                return brk;

            if (!SK.NormInf().IsOne())
                throw new("Only for ternary secret key");

            return brk = BRKgswBGV(SK, PK);
        }
    }

    private RLWECipher[] aks = [];

    public RLWECipher[] AutoMorhKeys
    {
        get
        {
            if (aks.Length != 0)
                return aks;

            return aks = N.SeqLazy()
                .Select(j => SWKBGV(PM, SK, SK.Substitute(PM.X.Pow(j)).ResModSigned(PM, T), T, PK.Q, SP.Pow(Level)))
                .ToArray();
        }
    }

    public RLWECipher FromRegevCipher(RegevCipher cipher, RLWECipher[] exsk)
    {
        var a = new Rational(cipher.B.Signed);
        return a - cipher.A.Zip(exsk).Select(e => new Rational(e.First.Signed) * e.Second).ToVec().Sum();
    }

    public RLWECipher[] FromRegevCipher(RegevCipher[] ciphers, RLWECipher[] exsk) =>
        ciphers.Select(c => FromRegevCipher(c, exsk)).ToArray();

    public RegevCipher ToRegevCipher(RLWECipher cipher, RLWECipher swk)
    {
        var cipher2 = SwitchKeyBGV(cipher, swk).ModSwitch(Q);
        var a = ExtractArr(cipher2.B, n).Select(e => new ZnInt64((long)T.Num, (long)e.Mod(T).Num)).ToVec();
        var b = new ZnInt64((long)T.Num, (long)cipher2.A[0].Mod(T).Num);
        return new(a, b);
    }

    public RegevCipher[] ToRegevCipher(RLWECipher[] ciphers, RLWECipher swk) =>
        ciphers.Select(c => ToRegevCipher(c, swk)).ToArray();

    public RLWECipher Bootstrapping(RLWECipher cipher)
    {
        if (!BootstrappingMode || !RLKS[cipher.Q].nextMod.IsOne())
            return cipher;

        return Bootstrapping(cipher, PK, AutoMorhKeys, BlindRotateKeys).ModSwitch(RLKS[PK.Q].nextMod);
    }

    public RLWECipher[] Bootstrapping(RLWECipher[] ciphers)
    {
        return ciphers.Select(Bootstrapping).ToArray();
    }

    public RLWECipher NOT(RLWECipher cipher) => Thalf - cipher;

    public RLWECipher AND(RLWECipher cipher1, RLWECipher cipher2)
    {
        var (c1, c2) = RLWECipher.AdjustLevel(cipher1, cipher2);
        if (!RLKS.ContainsKey(c1.Q))
            throw new($"Level 0 reached");

        var (nextMod, rlk) = RLKS[c1.Q];
        return InvThalf * MulRelinBGV(c1, c2, rlk).ModSwitch(nextMod);
    }

    public RLWECipher NAND(RLWECipher cipher1, RLWECipher cipher2) => NOT(AND(cipher1, cipher2));

    public RLWECipher NOR(RLWECipher cipher1, RLWECipher cipher2) => AND(NOT(cipher1), NOT(cipher2));

    public RLWECipher OR(RLWECipher cipher1, RLWECipher cipher2) => NOT(NOR(cipher1, cipher2));

    public RLWECipher XOR(RLWECipher cipher1, RLWECipher cipher2) =>
        AND(cipher1, NOT(cipher2)) + AND(NOT(cipher1), cipher2);

    public RLWECipher[] NOT(RLWECipher[] ciphers) => ciphers.Select(NOT).ToArray();

    public RLWECipher[] AND(RLWECipher[] ciphers1, RLWECipher[] ciphers2)
    {
        return ciphers1.Zip(ciphers2).Select(e => AND(e.First, e.Second)).ToArray();
    }

    public RLWECipher[] NAND(RLWECipher[] ciphers1, RLWECipher[] ciphers2) => NOT(AND(ciphers1, ciphers2));

    public RLWECipher[] NOR(RLWECipher[] ciphers1, RLWECipher[] ciphers2) => AND(NOT(ciphers1), NOT(ciphers2));

    public RLWECipher[] OR(RLWECipher[] ciphers1, RLWECipher[] ciphers2) => NOT(NOR(ciphers1, ciphers2));

    public RLWECipher[] XOR(RLWECipher[] ciphers1, RLWECipher[] ciphers2)
    {
        return ciphers1.Zip(ciphers2).Select(e => XOR(e.First, e.Second)).ToArray();
    }

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

    public RLWECipher[] OP(string name, RLWECipher[] ciphers1, RLWECipher[] ciphers2)
    {
        return ciphers1.Zip(ciphers2).Select(e => OP(name, e.First, e.Second)).ToArray();
    }

    public RLWECipher[] ADD(RLWECipher[] ciphers1, RLWECipher[] ciphers2)
    {
        var carry = PK.Zero;
        var sum = new List<RLWECipher>();
        for (int i = 0; i < ciphers1.Length; i++)
        {
            var (xor, and, or) = (XOR(ciphers1[i], ciphers2[i]), AND(ciphers1[i], ciphers2[i]),
                OR(ciphers1[i], ciphers2[i]));
            
            carry = Bootstrapping(carry);
            sum.Add(XOR(xor, carry));
            carry = AND(or, carry);
            
            carry = Bootstrapping(carry);
            carry = OR(and, carry);
        }

        return sum.ToArray();
    }
}
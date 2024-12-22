using System.Numerics;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Integers;

namespace FastGoat.UserGroup.Lattice;

public class BGV
{
    public KPoly<Rational> PM { get; }
    public (KPoly<Rational> pk0, KPoly<Rational> pk1) PK { get; }
    public KPoly<Rational> SK { get; }
    public int N { get; }
    public Rational P { get; }
    public Rational Q { get; }
    public BGVPublic BGVpublic { get; }
    public (KPoly<Rational> ek0, KPoly<Rational> ek1) EK { get; }
    public Dictionary<int, (KPoly<Rational> ct0, KPoly<Rational> ct1)> AK { get; }

    public BGV(int n, int p, BigInteger q) :
        this(n, p, q, DistributionExt.DiceSample(n, [-1, 0, 1]).ToKPoly(Rational.KZero()))
    {
    }

    private BGV(int n, int p, BigInteger q, KPoly<Rational> sk)
    {
        (N, Q, P) = (n, new(q), new(p));
        SK = sk;
        
        var e0 = DistributionExt.DiscreteGaussianSample(N, 0, 3.0).ToKPoly(Rational.KZero());
        var e1 = DistributionExt.DiscreteGaussianSample(N, 0, 3.0).ToKPoly(Rational.KZero());
        if (BGVPublic.SafeMode)
            e0 = e1 = e0.Zero;

        var pk1 = DistributionExt.DiceSampleBigInt(N, 0, Q.Num - 1).Select(e => new Rational(e)).ToKPoly();
        var ek1 = DistributionExt.DiceSampleBigInt(N, 0, Q.Num - 1).Select(e => new Rational(e)).ToKPoly();
        
        var x = FG.QPoly();
        PM = x.Pow(n) + 1;
        var pk0 = (P * e0 + pk1 * SK).ResMod(PM).CoefsMod(Q);
        var ek0 = (P * e1 + ek1 * SK + SK.Pow(2)).ResMod(PM).CoefsMod(Q);
        PK = (pk0, pk1);
        EK = (ek0, ek1);

        AK = FG.UnInt(2 * N).Order()
            .Select(i => (i.K, BGVPublic.Encrypt(SK.Substitute(x.Pow(i.K)).ResMod(PM).CoefsMod(P), P, Q, PM, PK)))
            .ToDictionary(e => e.K, e => e.Item2);
        
        BGVpublic = new(N, P, Q, PM, PK, EK, AK);
    }

    public (KPoly<Rational> ct0, KPoly<Rational> ct1) SwitchKeyGen(BGVPublic bgvPublic) => bgvPublic.Encrypt(SK);

    public ((KPoly<Rational> ct0, KPoly<Rational> ct1) plus, (KPoly<Rational> ct0, KPoly<Rational> ct1) minus)[] 
        BlindRotateKeyGen(Rational Qi, Rational Qf)
    {
        var enc0 = () => BGVPublic.Encrypt(PM.Zero, Qi, Qf, PM, PK); // encrypt zero
        var enc1 = () => BGVPublic.Encrypt(PM.One, Qi, Qf, PM, PK); // encrypt one
        var n = PM.Degree;
        return n.SeqLazy().Select(i => SK[i])
            .Select(c => (plus: ((c - 1).IsZero() ? enc1() : enc0()), minus: ((c + 1).IsZero() ? enc1() : enc0())))
            .ToArray();
    }

    public KPoly<Rational> Decrypt((KPoly<Rational> ct0, KPoly<Rational> ct1) cypher) => Decrypt(cypher, P);

    public KPoly<Rational> Decrypt((KPoly<Rational> ct0, KPoly<Rational> ct1) cypher, Rational mod)
    {
        return (cypher.ct0 - cypher.ct1 * SK).ResMod(PM).CoefsMod(mod);
    }

    public void Show(bool showKeys = true)
    {
        Console.WriteLine($"N = {N} Q = {Q} P = {P}");
        if (showKeys)
        {
            Console.WriteLine("Private Key");
            Console.WriteLine(SK);
            Console.WriteLine("Public Key");
            Console.WriteLine(BGVpublic.PK.pk0);
            Console.WriteLine(BGVpublic.PK.pk1);
            Console.WriteLine("Evaluation Key");
            Console.WriteLine(BGVpublic.EK.ek0);
            Console.WriteLine(BGVpublic.EK.ek1);
        }

        Console.WriteLine();
    }
}
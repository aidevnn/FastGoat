using System.Numerics;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Lattice;
using Rq = FastGoat.Structures.VecSpace.KPoly<FastGoat.UserGroup.Integers.Rational>;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");
IntExt.RngSeed(7532159);
Ring.DisplayPolynomial = MonomDisplay.StarCaret;

// ZnBInt-PadicSequence()
IEnumerable<BigInteger> DecompBase(BigInteger a, int b = 2, int size = -1, bool sizeStrict = false)
{
    var a0 = a;
    var k = 0;
    do
    {
        (a0, var r) = BigInteger.DivRem(a0, b);
        yield return r;
        ++k;
    } while (!a0.IsZero || (sizeStrict && k < size));
}

Vec<Rational> VecBase(int k, int B = 2)
{
    return new(k.SeqLazy().Select(i => new Rational(BigInteger.Pow(B, i))).ToArray());
}

Vec<Rq> DecomposeRq(Rq poly, int degree, int k, int B = 2)
{
    var extr = degree.SeqLazy().Select(i => DecompBase(poly[i].Num, B, k, sizeStrict: true).ToArray()).ToArray();
    // extr.Println(l => $"[{l.Glue(", ")}]", "Decompose);
    return k.SeqLazy().Select(j => degree.SeqLazy().Select(i => new Rational(extr[i][j])).ToKPoly()).ToVec();
}

Vec<BGVCipher> RlwePrime(BGVCipher cipher)
{
    var q = cipher.Q;
    if (!BigInteger.IsPow2(q.Num))
        throw new();

    var d = (int)BigInteger.Log2(q.Num);
    return VecBase(d).MulA<Rq, BGVCipher>(cipher);
}

(Vec<BGVCipher> csm, Vec<BGVCipher> cm) Rgsw(Rq m, Rq pm, Rational t, Rational q, BGVCipher pk)
{
    var (csm, cm) = FHE.EncryptRgswBGV(m, pm, t, q, pk);
    return (RlwePrime(csm), RlwePrime(cm));
}

BGVCipher RqByRlwePrime(Rq poly, Vec<BGVCipher> vCipher)
{
    var pm = vCipher[0].PM;
    var vPoly = DecomposeRq(poly, pm.Degree, vCipher.Length);
    return vPoly.MulA(vCipher).Sum();
}

BGVCipher RlweByRgsw(BGVCipher cipher, (Vec<BGVCipher> csm, Vec<BGVCipher> cm) gCipher)
{
    return RqByRlwePrime(cipher.A, gCipher.cm) - RqByRlwePrime(cipher.B, gCipher.csm);
}

{
    // FHE.NoiseOff();
    var l = 3;
    var n = 1 << l;
    var t0 = n;
    var q0 = n * n;
    var (pm, sk, t, q, pk, rlk) = FHE.KeyGenBGV(n, t0, q0);
    FHE.Show(pm, sk, t, q, pk, rlk);
    
    for(int i = 0; i < 50; ++i)
    {
        var m1 = FHE.GenUnif(n, t);
        var m2 = FHE.GenUnif(n, t);
        var m1m2 = (m1 * m2).ResMod(pm, t);

        var cm1 = FHE.EncryptBGV(m1, pm, t, q, pk);
        var cm2 = FHE.EncryptBGV(m2, pm, t, q, pk);
        var cm1m2 = FHE.EncryptBGV(m1m2, pm, t, q, pk);
        cm1.Show($"ct m1:{m1}");
        cm2.Show($"ct m2:{m2}");

        var cm1p = RlwePrime(cm1);
        cm1p.Println("RlwePrime of m1");
        var gm1 = Rgsw(m1, pm, t, q, pk);
        gm1.csm.Println("Rgsw s * m1");
        gm1.cm.Println("Rgsw m1");

        var m2_cm1p = RqByRlwePrime(m2, cm1p);
        m2_cm1p.Show("m2_cm1p");
        var d_m2_cm1p = FHE.DecryptBGV(m2_cm1p, pm, sk, t);
        (m2 * cm1).Show("m2 * cm1");
        Console.WriteLine($"decrypt:{d_m2_cm1p}");
        Console.WriteLine($"       :{FHE.DecryptBGV(m2 * cm1, pm, sk, t)}");
        Console.WriteLine($"m1 * m2:{m1m2}");

        var cm2_cm1p = RlweByRgsw(cm2, gm1);
        cm1m2.Show("cm1m2");
        cm2_cm1p.Show("cm2_cm1p");
        var d_cm2_cm1p = FHE.DecryptBGV(cm2_cm1p, pm, sk, t);
        Console.WriteLine($"decrypt:{d_cm2_cm1p}");
        Console.WriteLine();
        if (!d_m2_cm1p.Equals(m1m2) || !d_cm2_cm1p.Equals(m1m2))
            throw new($"step[{i}]");
    }
}
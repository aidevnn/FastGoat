using System.Numerics;
using System.Runtime.Intrinsics;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Lattice;
using FastGoat.UserGroup.Polynoms;
using Rq = FastGoat.Structures.VecSpace.KPoly<FastGoat.UserGroup.Integers.Rational>;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");
// IntExt.RngSeed(7532159);
Ring.DisplayPolynomial = MonomDisplay.StarCaret;

// ZnBInt-PadicSequence()
IEnumerable<int> DecompBase(BigInteger a, int b)
{
    var a0 = a;
    do
    {
        (a0, var r) = BigInteger.DivRem(a0, b);
        yield return (int)r;
    } while (!a0.IsZero);
}

Rational[] ExtractArr(int i, Rq P, int N)
{
    return N.SeqLazy(start: i, step: -1).Select(j => j >= 0 ? P[j] : -P[N + j]).ToArray();
}

(Rational ai, Rational[] bi)[] Extract(Rq a, Rq b, int N)
{
    return N.SeqLazy().Select(i => (ai: a[i], bi: ExtractArr(i, b, N).ToArray())).ToArray();
}

Rational InnerProd(Rational[] m0, Rational[] m1, Rational mod)
    => m0.Zip(m1).Aggregate(Rational.KZero(), (acc, e) => (acc + e.First * e.Second).Mod(mod));

(Rq acc0, Rq acc1) BlindRotate(BGVPublic pub, Rq f, (Rational ai, Rational[] bi) ab, Rational Q)
{
    var (N, P, PM) = (new Rational(pub.N), pub.P, pub.PM);

    var n = PM.Degree;
    var x = PM.X;

    var rgsw0 = BGVPublic.Encrypt(PM.One, n, pub.P, Q, PM, pub.PK);
    var (alpha, beta) = ab;
    beta = beta.Select(c => c.Opp()).ToArray();
    var (acc0, acc1) = BGVPublic.Encrypt((f * x.Pow((int)alpha.Mod(N).Num)).ResMod(PM).CoefsMod(Q),
        n, pub.P, Q, PM, pub.PK);
    for (int i = 0; i < n; i++)
    {
        var (rgsw_sp, rgsw_sm) = pub.BRK[i];

        var ai = (int)beta[i].Mod(N).Num;
        var xai = (x.Pow(ai) - 1).ResMod(PM).CoefsMod(P);
        var cxai = pub.Mul(rgsw_sp, BGVPublic.Encrypt(xai, n, pub.P, Q, PM, pub.PK));

        var _ai = (int)beta[i].Opp().Mod(N).Num;
        var x_ai = (x.Pow(_ai) - 1).ResMod(PM).CoefsMod(P);
        var cx_ai = pub.Mul(rgsw_sm, BGVPublic.Encrypt(x_ai, n, pub.P, Q, PM, pub.PK));

        var c0 = (rgsw0.ct0 + cxai.ct0 + cx_ai.ct0).ResMod(PM).CoefsMod(Q);
        var c1 = (rgsw0.ct1 + cxai.ct1 + cx_ai.ct1).ResMod(PM).CoefsMod(Q);
        (acc0, acc1) = pub.Mul((acc0, acc1), (c0, c1));
    }

    return (acc0, acc1);
}

{
    // BGVPublic.SafeMode = true;
    var l = 4;
    var n = 1 << l;
    var bgv = new BGV(n, p: n, q: 2 * n);
    bgv.Show();
    var pubKeys = bgv.BGVpublic;
    var (_, P, Q, PM) = pubKeys.Params_NPQ_PM;

    var x = PM.X;
    var c = n / 4;
    var f = (2 * c + 1).SeqLazy(-c).Select(j => j * x.Pow(IntExt.AmodP(j, 2 * n)))
        .Aggregate(x.Zero, (acc, v) => acc + v).ResMod(PM).CoefsMod(Q);
    for (int k = 0; k < 50; ++k)
    {
        var ai = new Rational(IntExt.Rng.Next(n));
        var bi = n.SeqLazy().Select(_ => IntExt.Rng.Next(n) * ai.One).ToArray();
        var (acc0, acc1) = BlindRotate(pubKeys, f, (ai, bi), Q);
        var actual = bgv.Decrypt((acc0, acc1));

        Console.WriteLine($"f           :{f}");
        Console.WriteLine($"f           :[{f.Coefs.Glue(", ", "{0,2}")}]");
        Console.WriteLine($"ai          :{ai}");
        Console.WriteLine($"bi          :[{bi.Glue(", ", "{0,2}")}]");
        Console.WriteLine($"blind rotate:{actual}");
        
        // Testing result
        var s = n.SeqLazy().Select(i => bgv.SK[i]).ToArray();
        var u = (ai - InnerProd(bi, s, Q)).Mod(P);
        var expected = (x.Pow((int)u.Num) * f).ResMod(PM).CoefsMod(P);
        var factor = (2 * n).SeqLazy(1).First(k0 => ((k0 * actual).CoefsMod(P)).Equals(expected));
        Console.WriteLine($"u= a - <b,s>:{u}");
        Console.WriteLine($"f*X^{u,-2}      :{expected}");
        Console.WriteLine($"factor      :{factor}");
        Console.WriteLine();
    }
}

// f           :x^15 + 2*x^14 + 3*x^13 + 4*x^12 + 4*x^4 + 3*x^3 + 2*x^2 + x
// f           :[ 0,  1,  2,  3,  4,  0,  0,  0,  0,  0,  0,  0,  4,  3,  2,  1]
// ai          :4
// bi          :[ 0,  8, 13,  6,  8,  0,  0, 11, 11,  2,  6,  1,  0,  1,  4,  9]
// blind rotate:3*x^15 + 2*x^14 + x^13 + 15*x^11 + 14*x^10 + 13*x^9 + 12*x^8 + 12
// u= a - <b,s>:12
// f*X^12      :3*x^15 + 2*x^14 + x^13 + 15*x^11 + 14*x^10 + 13*x^9 + 12*x^8 + 12
// factor      :1
// 
// f           :x^15 + 2*x^14 + 3*x^13 + 4*x^12 + 4*x^4 + 3*x^3 + 2*x^2 + x
// f           :[ 0,  1,  2,  3,  4,  0,  0,  0,  0,  0,  0,  0,  4,  3,  2,  1]
// ai          :6
// bi          :[ 9,  2,  7,  0,  4,  1,  6, 13,  1,  8,  8,  6,  9,  4,  3, 15]
// blind rotate:12*x^9 + 13*x^8 + 14*x^7 + 15*x^6 + x^4 + 2*x^3 + 3*x^2 + 4*x
// u= a - <b,s>:5
// f*X^5       :4*x^9 + 3*x^8 + 2*x^7 + x^6 + 15*x^4 + 14*x^3 + 13*x^2 + 12*x
// factor      :15
// 
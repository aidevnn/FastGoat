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

(Rq acc0, Rq acc1) BlindRotate(BGV bgv, int c, (Rational ai, Rational[] bi) ab)
{
    var pub = bgv.BGVpublic;
    var (N, Q, PM) = (new Rational(pub.N), pub.Q, pub.PM);
    var (z, o) = (PM.Zero, PM.One);
    var rgsw0 = (ct1: z, ct2: o);
    var brk = bgv.BRK;
    var Q1 = Q / (2 * N);

    var n = PM.Degree;
    var x = PM.X;

    var (alpha, beta) = ab;
    var f = (2 * c + 1).SeqLazy(-c).Select(j => Q1 * j * x.Pow(IntExt.AmodP(j, 2 * n)))
        .Aggregate(z, (acc, v) => acc + v).ResMod(PM).CoefsMod(Q);
    var acc0 = x.Zero;
    var acc1 = (f * x.Pow((int)alpha.Mod(2 * N).Num)).ResMod(PM).CoefsMod(Q);
    // Console.WriteLine($"f = {f}");
    // Console.WriteLine($"  = [{f.CoefsExtended.Glue(", ")}]");
    for (int i = 0; i < n; i++)
    {
        var ai = beta[i].Mod(2 * N);
        var _ai = (-ai).Mod(2 * N);

        var (rgsw_sp, rgsw_sm) = brk[i];
        var xai = (x.Pow((int)ai.Num) - 1).ResMod(PM).CoefsMod(Q);
        var x_ai = (x.Pow((int)_ai.Num) - 1).ResMod(PM).CoefsMod(Q);
        var c0 = (rgsw0.ct1 + xai * rgsw_sp.ct1 + x_ai * rgsw_sm.ct1).ResMod(PM).CoefsMod(Q);
        var c1 = (rgsw0.ct2 + xai * rgsw_sp.ct2 + x_ai * rgsw_sm.ct2).ResMod(PM).CoefsMod(Q);
        (acc0, acc1) = pub.Mul((acc0, acc1), (c0, c1));
    }

    return (acc0, acc1);
}

// void exp()
// {
//     var (n, p, q) = (16, 4, 2.Pow(10));
//     var bfv = new BFV(n, p, q);
//     bfv.Show();
//     var pub = bfv.BFVpublic;
//     var (N, Q, P, PM) = (new Rational(n), bfv.Q, bfv.P, pub.PM);
//
//     var m = BFVPublic.GenUnif(n, p);
//     var s = bfv.SK;
//
//     var ct = pub.Encrypt(m);
//     var (b, a) = ct;
//     var e = ((b + a * s).ResMod(PM) - pub.FloorQoverP * m).CoefsMod(Q);
//     Console.WriteLine($"m:{m}");
//     Console.WriteLine($"e:{e}");
//     Console.WriteLine();
//     Console.WriteLine($"b + a * s     = {(b + a * s).ResMod(PM).CoefsMod(Q)}");
//     Console.WriteLine($"e + m * [Q/P] = {(e + pub.FloorQoverP * m).CoefsMod(Q)}");
//     Console.WriteLine();
//
//     var Q1 = Q / (2 * N);
//     var ct1 = (ct0: (P * ct.ct0).CoefsMod(Q), ct1: (P * ct.ct1).CoefsMod(Q));
//     var ct2 = (ct0: ct1.ct0.CoefsMod(Q1), ct1: ct1.ct1.CoefsMod(Q1));
//     var ct_prep = (ct0: ((ct1.ct0 - ct2.ct0) / Q1).CoefsMod(2 * N), ct1: ((ct1.ct1 - ct2.ct1) / Q1).CoefsMod(2 * N));
//
//     var (b2, a2) = ct2;
//     var (b_prep, a_prep) = ct_prep;
//
//     var u = (-(b_prep + a_prep * s).ResMod(PM).RoundPoly()).CoefsMod(2 * N);
//     Console.WriteLine($"u:{u}");
//     Console.WriteLine();
//
//     var ue = (P * e + Q1 * u).CoefsMod(Q1);
//     Console.WriteLine($"ct2(s) = b2 + a2 * s    = {(b2 + a2 * s).ResMod(PM).CoefsMod(Q1)}");
//     Console.WriteLine($"       = t * e + Q1 * u = {ue}");
//     Console.WriteLine(((b2 + a2 * s).ResMod(PM) - ue).CoefsMod(Q1));
//     Console.WriteLine();
//
//     var extr = Extract(a_prep, b_prep, n);
//     var si = n.SeqLazy().Select(i => s[i]).ToArray();
//     var checkExtr = n.SeqLazy().Select(i => InnerProd(extr[i].ai, si, Q1) + extr[i].bi)
//         .Select((c, i) => ((-c).Mod(Q1), u[i], (c + u[i]).Mod(Q1)))
//         .ToArray();
//     checkExtr.Println($"Check Extraction b + a * s = -u : {checkExtr.All(c => c.Item3.IsZero())}");
//
//     Console.WriteLine("Blind Rotate");
//     var accs = n.SeqLazy().Select(i => BlindRotate(bfv, extr[i])).ToArray();
//
//     Console.WriteLine("Repacking");
//     var n2 = 2.Pow(int.Log2(n));
//     var x = PM.X;
//     var CT = new (Rq ct1, Rq ct2)[n, n];
//     for (int i = 0; i < n2; i++)
//         CT[i, n] = accs[i];
//
//     for (int k = n; k > 1; k /= 2)
//     {
//         for (int i = 0; i < k / 2; i++)
//         {
//             var (tmp00, tmp01) = CT[i, k];
//             var (tmp10, tmp11) = CT[i + k / 2, k];
//             var xpow = x.Pow(k / 2);
//             CT[i, k / 2] = (tmp00 + xpow * tmp10, tmp01 * xpow * tmp11);
//         }
//     }
// }

(Rq, Rq) Repacking(BGV bgv, (Rq acc0, Rq acc1)[] accs, int n)
{
    var (N, P, Q, PM) = bgv.BGVpublic.Params_NPQ_PM;
    var l = accs.Length;
    // Console.WriteLine(new { n, N, l });
    var x = PM.X;
    var CT = Ring.Matrix((x.Zero, x.Zero), N, N);
    var pub = bgv.BGVpublic;
    for (int i = 0; i < n; i++)
    {
        CT[i, n - 1] = accs[i];
        for (int j = 0; j < N / n; j++)
        {
            // Console.WriteLine($"[i, n - 1]{(i, n - 1)} [i, j * n]{(i, j * n)}");
            var xnj = x.Pow(n * j);
            var c = pub.Mul(pub.Encrypt(xnj), CT[i, j * n]);
            CT[i, n - 1] = pub.Add(c, CT[i, n - 1]);
        }
    }

    for (int k = n; k > 1; k /= 2)
    {
        for (int i = 0; i < k / 2; i++)
        {
            // Console.WriteLine($"[i, k]{(i, k)} [i + k / 2, k]{(i + k / 2, k)} [i, k / 2]{(i, k / 2)}");
            var xpow = x.Pow(k / 2);
            var c = pub.Mul(pub.Encrypt(xpow), CT[i + k / 2, k - 1]);
            CT[i, k / 2] = pub.Add(CT[i, k - 1], c);
            var crot = pub.EvalAuto(pub.Sub(CT[i, k - 1], c), 1 + 2 * N / k);
            CT[i, k / 2] = pub.Add(CT[i, k - 1], crot);
        }
    }

    // Ring.DisplayMatrix(CT);
    return CT[0, 1];
}

Rq Signed(Rq poly, Rational Q)
{
    return poly.CoefsMod(Q).Coefs.Select(c => c * 2 > Q ? c - Q : c).ToKPoly();
}

/***
 * General Bootstrapping Approach for RLWE-based Homomorphic Encryption
 *
 * Andrey Kim1 , Maxim Deryabin1 , Jieun Eom1 , Rakyong Choi1 , Yongwoo Lee1 ,
 * Whan Ghang1 , and Donghoon Yoo1
 */
(Rq ct0, Rq ct1) Bootstrapping(BGV bgv, (Rq ct0, Rq ct1) ct)
{
    var pub = bgv.BGVpublic;
    var (n, P, Q, PM) = bgv.BGVpublic.Params_NPQ_PM;
    Rational N = $"{n}";
    
    // var s = bgv.SK;
    // var m = bgv.Decrypt(ct);
    // var (a, b) = ct;
    // var ev = (((a - b * s).ResMod(PM) - m) / P).CoefsMod(Q);
    // Console.WriteLine($"m:{m}");
    // Console.WriteLine($"e:{ev}");
    // Console.WriteLine();
    // Console.WriteLine($"a - b * s  = {(a - b * s).ResMod(PM).CoefsMod(Q)}");
    // Console.WriteLine($"P * ev + m = {(P * ev + m).CoefsMod(Q)}");
    // Console.WriteLine();

    var Q1 = Q / (2 * N);
    var ct1 = (ct0: ct.ct0.CoefsMod(Q1), ct1: ct.ct1.CoefsMod(Q1));
    var ct_prep = (ct0: ((ct1.ct0 - ct.ct0) / Q1).CoefsMod(2 * N), ct1: ((ct1.ct1 - ct.ct1) / Q1).CoefsMod(2 * N));

    var (a_prep, b_prep) = ct_prep;

    // var (a2, b2) = ct1;
    // var u = (-a_prep + b_prep * s).ResMod(PM).CoefsMod(2 * N);
    // Console.WriteLine($"u:{u}");
    // Console.WriteLine();
    //
    // var v = ((P * ev - Q1 * u).CoefsMod(Q) / Q).CoefsMod(Q);
    // var e = (ev - Q * v / P).CoefsMod(Q);
    // var ue = (m + P * e + Q1 * u).CoefsMod(Q1);
    // Console.WriteLine($"q' = {Q1}");
    // Console.WriteLine($"ct2(s) = a2 - b2 * s    = {(a2 - b2 * s).ResMod(PM).CoefsMod(Q1)}");
    // Console.WriteLine($"       = P * e + Q1 * u = {ue}");
    // Console.WriteLine(((a2 - b2 * s).ResMod(PM) - ue).CoefsMod(Q1));
    // Console.WriteLine();

    // Extraction
    var extr = Extract(a_prep, b_prep, n);
    // var si = n.SeqLazy().Select(i => s[i]).ToArray();
    // var checkExtr = n.SeqLazy().Select(i => extr[i].ai - InnerProd(extr[i].bi, si, Q1))
    //     .Select((c, i) => (c.Mod(Q1), u[i].Opp().Mod(Q1), (c + u[i]).Mod(Q1)))
    //     .ToArray();
    // checkExtr.Println($"Check Extraction a - b * s = -u : {checkExtr.All(c => c.Item3.IsZero())}");
    // Console.WriteLine();

    // var delta = 2 * Rational.Sqrt(N, nbDigits: 4);
    // var gamma = (a - b * s).ResMod(PM).CoefsMod(Q).Coefs.Select(k0 => k0 * 2 >= Q ? k0 - Q : k0)
    //     .MaxBy(k0 => k0.Absolute).Absolute;
    // var c = (int)((delta + 1) / 2 + gamma / Q1).Floor.Num;
    // Console.WriteLine(new { N, delta, gamma, c });
    var c = n / 2;

    // Blind Rotate
    var accs = new List<(Rq acc0, Rq acc1)>();
    foreach (var (i, (ai, bi)) in extr.Select((ab, i) => (i, ab)))
    {
        // Console.WriteLine($"i:{i}");
        var acc = BlindRotate(bgv, c, (ai, bi));
        accs.Add(acc);
        
        // var bix = bi.ToKPoly();
        // var aix = ai * bix.One;
        // var ui = bgv.Decrypt((aix, bix));
        // Console.WriteLine($"[{acc.acc0.CoefsExtended.Glue(", ")}]");
        // Console.WriteLine($"[{acc.acc1.CoefsExtended.Glue(", ")}]");
        // Console.WriteLine($"ui:{ui}");
        // Console.WriteLine($"[{ui.CoefsExtended.Glue(", ")}]");
    }

    // Repacking
    var repack = (ct0: PM.One * Q1, ct1: PM.One * Q1);
    for (int k = 1; k <= int.Log2(bgv.N); k++)
    {
        var tmp = Repacking(bgv, accs.ToArray(), 2.Pow(k));
        repack = pub.Mul(repack, tmp);
    }

    Console.WriteLine("Repack");
    Console.WriteLine(repack);
    Console.WriteLine($"ctsm:{bgv.Decrypt(repack)}");
    // Console.WriteLine($"   :{(-Q1 * u).CoefsMod(P)}");
    var ctboot = pub.Add(ct1, repack);
    return ctboot;
}

{
    var (n, p, q) = (16, 4, 2.Pow(10));
    var bgv = new BGV(n, p, q);
    bgv.Show();
    var pub = bgv.BGVpublic;
    var (N, Q, P, PM) = (new Rational(n), bgv.Q, bgv.P, pub.PM);

    var m = BGVPublic.GenUnif(n, p);
    var ct = pub.Encrypt(m);
    var ctboot = Bootstrapping(bgv, ct);
    Console.WriteLine($"m :{m}");
    Console.WriteLine($"  :{bgv.Decrypt(ctboot)}");

    var err1 = Signed((ct.ct0 - bgv.SK * ct.ct1 - m).ResMod(PM), Q);
    var err2 = Signed((ctboot.ct0 - bgv.SK * ctboot.ct1 - m).ResMod(PM), Q);

    // TODO: check error
    Console.WriteLine($"err1:{err1}");
    Console.WriteLine($"norm:{err1.NormInf()}");
    Console.WriteLine($"err2:{err2}");
    Console.WriteLine($"norm:{err2.NormInf()}");
}
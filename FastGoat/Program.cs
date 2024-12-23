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

Rational InnerProd(Rational[] m0, Rational[] m1, Rational mod)
    => m0.Zip(m1).Aggregate(Rational.KZero(), (acc, e) => (acc + e.First * e.Second).Mod(mod));

void testBR()
{
    var l = 4;
    var n = 1 << l;
    var t0 = 2 * n;
    var q0 = t0 * 3.Pow(3);
    var (pm, sk, t, q, pk, rlk) = FHE.KeyGenBGV(n, t0, q0);
    FHE.Show(pm, sk, t, q, pk, rlk);

    var Q = 25 * t;
    var delta = q / t;
    Console.WriteLine(new { t, Q, delta });

    var brk = FHE.BRKBGV(pm, sk, t, q, pk);

    var x = pm.X;
    var c = n / 4;
    var f = (2 * c + 1).SeqLazy(-c).Select(j => delta * j * FHE.XpowA(j, pm, t))
        .Aggregate(x.Zero, (acc, v) => acc + v).ResMod(pm).CoefsMod(t);
    var nbDigits = $"{Q}".Length;
    var fmt = $"{{0,{nbDigits}}}";
    var set = new HashSet<int>();
    var s = n.SeqLazy().Select(i => sk[i]).ToArray();

    for (int k = 0; k < 50; ++k)
    {
        var ai = new Rational(IntExt.Rng.Next(1, t0));
        var bi = n.SeqLazy().Select(_ => IntExt.Rng.Next(t0) * ai.One).ToArray();
        var acc = FHE.BlindRotateBGV((ai, bi), f, pm, t, q, pk, rlk, brk);
        var actual = FHE.DecryptBGV(acc, pm, sk, t);

        Console.WriteLine($"f           :{f}");
        Console.WriteLine($"f           :[{f.CoefsExtended(n - 1).Glue(", ", fmt)}]");
        Console.WriteLine($"ai          :{ai}");
        Console.WriteLine($"bi          :[{bi.Glue(", ", fmt)}]");
        Console.WriteLine($"blind rotate:{actual}");
        Console.WriteLine($"            :[{actual.CoefsExtended(n - 1).Glue(", ", fmt)}]");

        // Testing result
        var u = (ai - InnerProd(bi, s, t)).Mod(t);
        var expected = (x.Pow((int)u.Num) * f).ResMod(pm).CoefsMod(t);
        Console.WriteLine($"u= a - <b,s>:{u}");
        Console.WriteLine($"f*X^{u,-2}      :{expected}");
        Console.WriteLine($"f*X^{u,-2}      :[{expected.CoefsExtended(n - 1).Glue(", ", fmt)}]");
        var factor = t0.SeqLazy(1).First(k0 => (k0 * actual).CoefsMod(t).Equals(expected));
        Console.WriteLine($"factor      :{factor}");
        Console.WriteLine();
        set.Add(factor);
    }

    Console.WriteLine($"Factors:{set.Order().Glue(", ", fmt)}");
}

BGVCipher RepackingBGV(BGVCipher[] accs, int n, Rq pm, Rational t, Rational q, BGVCipher pk, BGVCipher rlk,
    Dictionary<int, BGVCipher> ak)
{
    var N = pm.Degree;
    var x = pm.X;
    var CT = Ring.Matrix((x.One, x.One), N, N + 1);
    for (int i = 0; i < n; i++)
    {
        CT[i, n] = accs[i];
        for (int j = 1; j < N / n; j++)
        {
            // Console.WriteLine($"[i, n - 1]{(i, n - 1)} [i, j * n]{(i, j * n)}");
            var exnj = FHE.XpowA(n * j, pm, t);
            var c = FHE.KMulBGV(CT[i, i + j * n], exnj, pm, q);
            CT[i, n] = FHE.AddBGV(c, CT[i, n], q);
        }
    }

    for (int k = n; k > 1; k /= 2)
    {
        for (int i = 0; i < k / 2; i++)
        {
            // Console.WriteLine($"[i, k]{(i, k)} [i + k / 2, k]{(i + k / 2, k)} [i, k / 2]{(i, k / 2)}");
            var expowk2 = FHE.XpowA(k / 2, pm, t);
            var c0 = FHE.KMulBGV(CT[i + k / 2, k], expowk2, pm, q);
            CT[i, k / 2] = FHE.AddBGV(CT[i, k], c0, q);

            var c1 = FHE.KMulBGV(CT[i + k / 2, k], expowk2, pm, q);
            var crot = FHE.AutoMorphBGV(FHE.SubBGV(CT[i, k], c1, q), 1 + 2 * N / k, pm, t, q, pk, rlk, ak);
            CT[i, k / 2] = FHE.AddBGV(CT[i, k / 2], crot, q);
        }
    }

    return CT[0, 1];
}

BGVCipher Bootstrapping(BGVCipher ct, Rq pm, Rational t, Rational q, Rational Q, BGVCipher pk, BGVCipher rlk,
    (BGVCipher plus, BGVCipher minus)[] brk, Dictionary<int, BGVCipher> ak, int fact = 2)
{
    var n = pm.Degree;
    var q1 = q / (fact * n);
    if (!q1.IsInteger())
        throw new();
    
    var ct1 = ct.CoefsMod(q1);
    var ctprep = new BGVCipher((ct.A - ct1.A) / q1, (ct.B - ct1.B) / q1).CoefsMod(new(fact * n));
    
    // Step 1. Extraction
    var extract = FHE.Extract(ctprep, n);

    var delta = double.Sqrt(n * 4.0);
    var gamma = q / 4 - q1 / 2 * delta;
    var c = (int)(0.5 * (delta + 1) + gamma * q1.Inv());
    var f = (2 * c + 1).SeqLazy(-c).Where(j => j != 0).Select(j => -q1 * j * FHE.XpowA(j, pm, q))
        .Aggregate((v0, v1) => v0 + v1).ResMod(pm).CoefsMod(q);
    
    // Step 2. Blind Rotate
    var seqBR = new List<BGVCipher>();
    foreach (var (ab, i) in extract.Select((ab, i) => (ab, i)))
        seqBR.Add(FHE.BlindRotateBGV(ab, f, pm, t, Q, pk, rlk, brk));

    // Step 3. Repacking
    seqBR = seqBR.Select(cipher => new BGVCipher(cipher.A * pm.One, cipher.B * pm.One)).ToList();
    var ctsm = RepackingBGV(seqBR.ToArray(), n, pm, q1, Q, pk, rlk, ak);
    ctsm.Show("ctsm");
    
    return FHE.AddBGV(ctsm, ct1, Q);
}

{
    var l = 4;
    var n = 1 << l;
    var t0 = n;
    var q0 = 8 * n * n;
    var (pm, sk, t, q, pk, rlk) = FHE.KeyGenBGV(n, t0, q0);
    FHE.Show(pm, sk, t, q, pk, rlk);

    var Q = t * q;
    var brk = FHE.BRKBGV(pm, sk, t, Q, pk);
    var ak = FHE.AKBGV(pm, sk, t, Q, pk);
    
    for (int k = 0; k < 10; ++k)
    {
        var m1 = FHE.GenUnif(n, t);
        var cm1 = FHE.EncryptBGV(m1, pm, t, q, pk);
        cm1.Show($"ct m1:{m1}");
        
        var ctboot = Bootstrapping(cm1, pm, t, q, Q, pk, rlk, brk, ak);
        var dec = FHE.DecryptBGV(ctboot, pm, sk, t);
        
        ctboot.Show("ctboot");
        Console.WriteLine($"decrypt ctboot:{dec}");
        Console.WriteLine($"m1            :{m1}");
        Console.WriteLine();
        if (!dec.Equals(m1))
            throw new("decrypt");
    }
}

// m:28*x^15 + 8*x^14 + 30*x^13 + 12*x^11 + 9*x^10 + 10*x^9 + x^8 + x^7 + 6*x^6 + 19*x^5 + 17*x^4 + 17*x^3 + 7*x^2 + 29*x + 9
// ct N:16 t:32 q:2048
// A:1435*x^15 + 526*x^14 + 798*x^13 + 668*x^12 + 941*x^11 + 1790*x^10 + 1273*x^9 + 1178*x^8 + 1156*x^7 + 1789*x^6 + 1330*x^5 + 1899*x^4 + 1141*x^3 + 1186*x^2 + 1174*x + 1953
// B:1282*x^15 + 441*x^14 + 1580*x^13 + 1993*x^12 + 1708*x^11 + 1068*x^10 + 617*x^9 + 933*x^8 + 748*x^7 + 1189*x^6 + 166*x^5 + 669*x^4 + 1165*x^3 + 1153*x^2 + 362*x + 127
// 
// ctsm
// A:18713*x^15 + 14031*x^14 + 30948*x^13 + 15197*x^12 + 57846*x^11 + 45173*x^10 + 28700*x^9 + 13708*x^8 + 64966*x^7 + 54165*x^6 + 55504*x^5 + 46870*x^4 + 5075*x^3 + 39555*x^2 + 26418*x + 3277
// B:8695*x^15 + 13307*x^14 + 13379*x^13 + 40175*x^12 + 1182*x^11 + 52247*x^10 + 61497*x^9 + 36063*x^8 + 32775*x^7 + 19096*x^6 + 2926*x^5 + 24517*x^4 + 61409*x^3 + 28164*x^2 + 36019*x + 22713
// 
// ctboot Q=65536
// A:18740*x^15 + 14045*x^14 + 30978*x^13 + 15225*x^12 + 57891*x^11 + 45235*x^10 + 28757*x^9 + 13734*x^8 + 64970*x^7 + 54226*x^6 + 55554*x^5 + 46913*x^4 + 5128*x^3 + 39589*x^2 + 26440*x + 3310
// B:8697*x^15 + 13364*x^14 + 13423*x^13 + 40184*x^12 + 1226*x^11 + 52291*x^10 + 61538*x^9 + 36100*x^8 + 32819*x^7 + 19133*x^6 + 2964*x^5 + 24546*x^4 + 61422*x^3 + 28165*x^2 + 36061*x + 22776
// decrypt ctboot:28*x^15 + 8*x^14 + 30*x^13 + 12*x^11 + 9*x^10 + 10*x^9 + x^8 + x^7 + 6*x^6 + 19*x^5 + 17*x^4 + 17*x^3 + 7*x^2 + 29*x + 9
// m1            :28*x^15 + 8*x^14 + 30*x^13 + 12*x^11 + 9*x^10 + 10*x^9 + x^8 + x^7 + 6*x^6 + 19*x^5 + 17*x^4 + 17*x^3 + 7*x^2 + 29*x + 9
// 
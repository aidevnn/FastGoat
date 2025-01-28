using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.LWE;
using Rq = FastGoat.Structures.VecSpace.KPoly<FastGoat.UserGroup.Integers.Rational>;
using static FastGoat.Commons.IntExt;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");
// RngSeed(2598);
Ring.DisplayPolynomial = MonomDisplay.StarCaret;

void testSetupBGV(int N, int t0, int nbTests = 5)
{
    GlobalStopWatch.AddLap();
    var rlwe = new RLWE(N, t0);
    var (n, pm, sk, t, q0, q1, q2, pk, rlk) = rlwe;
    Console.WriteLine(rlwe.Params);

    var s2 = sk.Pow(2).ResModSigned(pm, t);
    var d_rk = RLWE.DecryptBGV(rlk.ModSwitch(q2, q1), sk);
    if (!(d_rk - s2).IsZero())
        throw new("RLK");

    for (int i = 0; i < nbTests; ++i)
    {
        var m1 = RLWE.GenUnif(n, t);
        var c1 = RLWE.EncryptBGV(m1, pk);
        var m2 = RLWE.GenUnif(n, t);
        var c2 = RLWE.EncryptBGV(m2, pk);
        var d1 = RLWE.DecryptBGV(c1, sk);
        var d2 = RLWE.DecryptBGV(c2, sk);

        if (!d1.Equals(m1) || !d2.Equals(m2))
            throw new($"encrypt decrypt");

        var add_m1m2 = (m1 + m2).ResModSigned(pm, t);
        var add_c1c2 = c1 + c2;
        var d_add = RLWE.DecryptBGV(add_c1c2, sk);
        var add_err = RLWE.ErrorsBGV(add_c1c2, sk);
        
        if (!d_add.Equals(add_m1m2))
            throw new("m1 + m2");

        var k = RLWE.GenUnif(n, t);
        var km1 = (k * m1).ResModSigned(pm, t);
        var kc1 = k * c1;
        var d_k1 = RLWE.DecryptBGV(kc1, sk);
        var k_err = RLWE.ErrorsBGV(kc1, sk);
        if (!d_k1.Equals(km1))
            throw new("k * m1");

        var mul_m1m2 = (m1 * m2).ResModSigned(pm, t);
        var mul_c1c2 = RLWE.MulRelinBGV(c1, c2, rlk);
        var d_mul = RLWE.DecryptBGV(mul_c1c2, sk);
        var mul_err = RLWE.ErrorsBGV(mul_c1c2, sk);

        if (!d_mul.Equals(mul_m1m2))
            throw new("m1 * m2");
        
        if (Logger.Level != LogLevel.Off)
        {
            Console.WriteLine($"m1:{m1}");
            Console.WriteLine($"  :{d1}");
            Console.WriteLine($"m2:{m2}");
            Console.WriteLine($"  :{d2}");
            Console.WriteLine();
            Console.WriteLine($"m1 + m2:{add_m1m2}");
            Console.WriteLine($"       :{d_add}");
            Console.WriteLine($"add_err:{add_err}");
            Console.WriteLine();
            Console.WriteLine($"k     :{k}");
            Console.WriteLine($"k * m1:{km1}");
            Console.WriteLine($"      :{d_k1}");
            Console.WriteLine($"k_err :{k_err}");
            Console.WriteLine();
            Console.WriteLine($"m1 * m2:{mul_m1m2}");
            Console.WriteLine($"       :{d_mul}");
            Console.WriteLine($"mul_err:{mul_err}");
            Console.WriteLine();
        }

    }

    GlobalStopWatch.Show($"Complete {nbTests} Tests (e/d,+,*,x)");
    Console.WriteLine();
}

{
    RecomputeAllPrimesUpTo(5000000);
    testSetupBGV(N: 2048, t0: 83);
    testSetupBGV(N: 1024, t0: 2);
    testSetupBGV(N: 512, t0: 7);
    testSetupBGV(N: 256, t0: 17);
    testSetupBGV(N: 128, t0: 37);
    testSetupBGV(N: 64, t0: 79);
    
    // Logger.Level = LogLevel.Level1;
    testSetupBGV(N: 32, t0: 163);
    testSetupBGV(N: 16, t0: 17);
}
// RLWE N=2048=2^11, Φ(N)=1024 PM=x^1024 + 1 t=83 q=[2210291, 4886120121293=2210291*2210623, 10797314041192422349=2209793*4886120121293]
// # Complete 5 Tests (e/d,+,*,x) Time:1m10s
// 
// RLWE N=1024=2^10, Φ(N)=512 PM=x^512 + 1 t=2 q=[12301, 151585223=12301*12323, 1862830805447=12289*151585223]
// # Complete 5 Tests (e/d,+,*,x) Time:8.801s
// 
// RLWE N=512=2^9, Φ(N)=256 PM=x^256 + 1 t=7 q=[10781, 116833697=10781*10837, 1256312743841=10753*116833697]
// # Complete 5 Tests (e/d,+,*,x) Time:2.255s
// 
// RLWE N=256=2^8, Φ(N)=128 PM=x^128 + 1 t=17 q=[26249, 690794933=26249*26317, 18038728085429=26113*690794933]
// # Complete 5 Tests (e/d,+,*,x) Time:672ms
// 
// RLWE N=128=2^7, Φ(N)=64 PM=x^64 + 1 t=37 q=[9547, 93264643=9547*9769, 883495963139=9473*93264643]
// # Complete 5 Tests (e/d,+,*,x) Time:164ms
// 
// RLWE N=64=2^6, Φ(N)=32 PM=x^32 + 1 t=79 q=[36341, 1343635793=36341*36973, 47555301621649=35393*1343635793]
// # Complete 5 Tests (e/d,+,*,x) Time:69ms
// 
// RLWE N=32=2^5, Φ(N)=16 PM=x^16 + 1 t=163 q=[11411, 152530837=11411*13367, 1591354222421=10433*152530837]
// # Complete 5 Tests (e/d,+,*,x) Time:32ms
// 
// RLWE N=16=2^4, Φ(N)=8 PM=x^8 + 1 t=17 q=[1429, 2187799=1429*1531, 2977594439=1361*2187799]
// # Complete 5 Tests (e/d,+,*,x) Time:7ms
// 
// 
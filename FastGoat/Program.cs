using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Integers;
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

{
    var rlwe = new RLWE(16);
    var (n, pm, sk, t, q0, q1, q2, pk, rlk) = rlwe;
    Console.WriteLine(rlwe.Params);
    
    var s2 = sk.Pow(2).ResMod(pm, t);
    var d_rk = RLWE.DecryptBGV(rlk.ModSwitch(q2, q1), sk);
    if (!(d_rk - s2).IsZero())
        throw new("RLK");

    for(int i = 0; i < 100; ++i)
    {
        var m1 = RLWE.GenUnif(n, t);
        var c1 = RLWE.EncryptBGV(m1, pk);
        var m2 = RLWE.GenUnif(n, t);
        var c2 = RLWE.EncryptBGV(m2, pk);
        var d1 = RLWE.DecryptBGV(c1, sk);
        var d2 = RLWE.DecryptBGV(c2, sk);

        Console.WriteLine($"m1:{m1}");
        Console.WriteLine($"  :{d1}");
        Console.WriteLine($"m2:{m2}");
        Console.WriteLine($"  :{d2}");
        Console.WriteLine();
        if (!d1.Equals(m1) || !d2.Equals(m2))
            throw new($"encrypt decrypt");
        
        var add_m1m2 = (m1 + m2).ResMod(pm, t);
        var add_c1c2 = c1 + c2;
        var d_add = RLWE.DecryptBGV(add_c1c2, sk);
        Console.WriteLine($"m1 + m2:{add_m1m2}");
        Console.WriteLine($"       :{d_add}");
        Console.WriteLine();
        
        if (!d_add.Equals(add_m1m2))
            throw new("m1 + m2");
        
        var k = RLWE.GenUnif(n, t);
        var km1 = (k * m1).ResMod(pm, t);
        var kc1 = k * c1;
        var d_k1 = RLWE.DecryptBGV(kc1, sk);
        var k_err = RLWE.ErrorsBGV(kc1, sk);
        var k_err2 = RLWE.ErrorsBGV(kc1.ModSwitch(q1, q0), sk);
        Console.WriteLine($"k     :{k}");
        Console.WriteLine($"k * m1:{km1}");
        Console.WriteLine($"      :{d_k1}");
        Console.WriteLine($"k_err :{k_err}");
        Console.WriteLine($"      :{k_err2}");
        Console.WriteLine();
        if (!d_k1.Equals(km1))
            throw new("k * m1");
        
        var mul_m1m2 = (m1 * m2).ResMod(pm, t);
        var mul_c1c2 = RLWE.MulRelinBGV(c1, c2, rlk);
        var d_mul = RLWE.DecryptBGV(mul_c1c2, sk);
        var mul_err = RLWE.ErrorsBGV(mul_c1c2, sk);
        var mul_err2 = RLWE.ErrorsBGV(mul_c1c2.ModSwitch(q1, q0), sk);
        Console.WriteLine($"m1 * m2:{mul_m1m2}");
        Console.WriteLine($"       :{d_mul}");
        Console.WriteLine($"mul_err:{mul_err}");
        Console.WriteLine($"       :{mul_err2}");
        Console.WriteLine();
        
        if (!d_mul.Equals(mul_m1m2))
            throw new("m1 * m2");
    }

    Console.WriteLine("Complete !!!");
}
// RLWE N=16=2^4, Φ(N)=8 PM=x^8 + 1 t=17 q=[103, 3372529=103*32743, 1035366403=307*3372529]
// m1:2*x^7 + 6*x^6 + 13*x^5 + 9*x^4 + 10*x^3 + 13*x^2 + 14*x + 7
//   :2*x^7 + 6*x^6 + 13*x^5 + 9*x^4 + 10*x^3 + 13*x^2 + 14*x + 7
// m2:9*x^7 + 4*x^6 + 9*x^5 + 12*x^4 + 6*x^2 + 9*x + 10
//   :9*x^7 + 4*x^6 + 9*x^5 + 12*x^4 + 6*x^2 + 9*x + 10
// 
// m1 + m2:11*x^7 + 10*x^6 + 5*x^5 + 4*x^4 + 10*x^3 + 2*x^2 + 6*x
//        :11*x^7 + 10*x^6 + 5*x^5 + 4*x^4 + 10*x^3 + 2*x^2 + 6*x
// 
// k     :12*x^7 + 8*x^6 + 2*x^5 + 3*x^4 + 5*x^3 + 2*x^2 + 2*x + 13
// k * m1:4*x^7 + 16*x^6 + 2*x^5 + 12*x^4 + 6*x^3 + 16*x^2 + 5*x + 14
//       :4*x^7 + 16*x^6 + 2*x^5 + 12*x^4 + 6*x^3 + 16*x^2 + 5*x + 14
// k_err :-833*x^7 + 5100*x^6 + 3808*x^5 - 3774*x^4 + 782*x^3 + 408*x^2 + 272*x + 1887
//       :-34*x^6 - 34*x^4 - 17*x^3 - 102*x^2 - 102*x + 170
// 
// m1 * m2:15*x^7 + 13*x^6 + 15*x^5 + 13*x^4 + 5*x^3 + x + 14
//        :15*x^7 + 13*x^6 + 15*x^5 + 13*x^4 + 5*x^3 + x + 14
// mul_err:-14773*x^7 - 510068*x^6 - 51459*x^5 + 394451*x^4 + 60605*x^3 - 326485*x^2 + 29988*x - 11305
//        :68*x^7 - 17*x^6 - 34*x^5 - 85*x^3 - 136*x^2 + 187*x - 17
// 
// m1:9*x^7 + 14*x^6 + 3*x^5 + 7*x^3 + 2*x^2 + 9
//   :9*x^7 + 14*x^6 + 3*x^5 + 7*x^3 + 2*x^2 + 9
// m2:10*x^7 + 14*x^6 + 10*x^5 + 7*x^4 + x^3 + 13*x^2 + 14*x + 14
//   :10*x^7 + 14*x^6 + 10*x^5 + 7*x^4 + x^3 + 13*x^2 + 14*x + 14
// 
// m1 + m2:2*x^7 + 11*x^6 + 13*x^5 + 7*x^4 + 8*x^3 + 15*x^2 + 14*x + 6
//        :2*x^7 + 11*x^6 + 13*x^5 + 7*x^4 + 8*x^3 + 15*x^2 + 14*x + 6
// 
// k     :3*x^7 + 3*x^6 + 14*x^5 + 14*x^4 + 3*x^3 + 7*x + 10
// k * m1:x^7 + 6*x^6 + 8*x^5 + 15*x^4 + x^3 + 6*x^2 + 3*x + 16
//       :x^7 + 6*x^6 + 8*x^5 + 15*x^4 + x^3 + 6*x^2 + 3*x + 16
// k_err :255*x^7 + 1054*x^6 - 1428*x^5 - 4964*x^4 - 5508*x^3 - 357*x^2 + 1258*x + 3604
//       :17*x^7 - 34*x^6 + 68*x^5 + 85*x^4 - 17*x^2 - 34
// 
// m1 * m2:10*x^7 + 6*x^6 + 10*x^5 + 7*x^4 + 9*x^3 + 6*x^2 + 9*x + 6
//        :10*x^7 + 6*x^6 + 10*x^5 + 7*x^4 + 9*x^3 + 6*x^2 + 9*x + 6
// mul_err:-102493*x^7 + 514352*x^6 - 380460*x^5 - 28390*x^4 - 691050*x^3 - 319957*x^2 + 87057*x - 112812
//        :-34*x^7 - 17*x^5 + 102*x^4 - 34*x^3 + 85*x^2 + 102*x - 17
// ...
// Complete !!!
// 
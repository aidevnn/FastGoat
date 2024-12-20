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

Rq Signed(Rq poly, Rational Q) => poly.Coefs.Select(c => 2 * c > Q ? c - Q : c).ToKPoly();

{
    var (n, p) = (32, 8); // n = 16, 32, 64
    var q = n * n * p;
    var bgv = new BGV(n, p, q);
    bgv.Show();
    var pub = bgv.BGVpublic;
    var (N, Q, P, PM) = (new Rational(n), bgv.Q, bgv.P, pub.PM);

    var nbFails = 0;
    var nbTests = 100;
    var nbTry = 100;
    var Q1 = P;
    var errNormCan = (min: q * 100.0, max: 0.0);
    var errNormInfty = (min: q * 100.0, max: 0.0);
    for (int i = 0; i < nbTests; ++i)
    {
        var m1 = BGVPublic.GenUnif(n, p);
        var m2 = BGVPublic.GenUnif(n, p);
        var m3 = (m1 * m2).ResMod(PM).CoefsMod(P);
        var ct1 = pub.Encrypt(m1);
        var ct2 = pub.Encrypt(m2);
        var ct = pub.Mul(ct1, ct2);

        var ctmin = ct;
        var min = 100.0 * q;
        for (int k = 0; k < nbTry; ++k)
        {
            ctmin = pub.CoefsMod(ctmin, Q1);
            var c0 = pub.Add(ctmin, pub.Encrypt(m1.Zero));
            var norm = BGVPublic.NormCan(c0, pub.Roots, Q);
            if (norm < min)
            {
                min = norm;
                ctmin = c0;
            }
        }
        
        var d = bgv.Decrypt(ctmin);
        Console.WriteLine($"m  :{m3}");
        Console.WriteLine($"   :{d}");
        if (!d.Equals(m3))
            throw new();

        var diff_init = (ct.ct0 - bgv.SK * ct.ct1 - m3).ResMod(PM).CoefsMod(Q);
        var norm_diff_init = pub.NormCan(diff_init);
        var diff_rand = (ctmin.ct0 - bgv.SK * ctmin.ct1 - m3).ResMod(PM).CoefsMod(Q);
        var norm_diff_rand = pub.NormCan(diff_rand);
        Console.WriteLine($"err init      :{pub.Signed(diff_init)}");
        Console.WriteLine($"err end       :{pub.Signed(diff_rand)}");
        Console.WriteLine($"Q1      :{Q1}");
        Console.WriteLine($"Q/2     :{Q / 2}");
        Console.WriteLine($"norm canonic");
        Console.WriteLine($"||err ct||    :{norm_diff_init}");
        Console.WriteLine($"||err ctmin|| :{norm_diff_rand}");
        Console.WriteLine($"norm infty");
        Console.WriteLine($"||err ct||    :{pub.NormInf(diff_rand)}");
        Console.WriteLine($"||err ctmin|| :{pub.NormInf(diff_rand)}");
        if (norm_diff_rand > 0.5 * q)
        {
            Console.WriteLine($"Q/2       :{Q / 2}");
            Console.WriteLine("Fail");
            ++nbFails;
        }

        errNormCan = (double.Min(errNormCan.min, pub.NormCan(diff_rand)), double.Max(errNormCan.max, pub.NormCan(diff_rand)));
        errNormInfty=(double.Min(errNormInfty.min, pub.NormInf(diff_rand)), double.Max(errNormInfty.max, pub.NormInf(diff_rand)));
        Console.WriteLine();
    }
    
    Console.WriteLine($"Nb Fails:{nbFails}/{nbTests}");
    Console.WriteLine($"norm canonic");
    Console.WriteLine($"diff min:{errNormCan.min}");
    Console.WriteLine($"diff max:{errNormCan.max}");
    Console.WriteLine($"norm infty");
    Console.WriteLine($"diff min:{errNormInfty.min}");
    Console.WriteLine($"diff max:{errNormInfty.max}");
}

// N = 32 Q = 8192 P = 8
// 
// Q/2     :4096
// norm canonic
// ||err ct||    :29316.00947965181
// ||err ctmin|| :293.07944187411
// 
// Nb Fails:0/100
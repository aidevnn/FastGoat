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

{
    var l = 4;
    var n = 1 << l;
    var t0 = n;
    var q0 = 8 * n * n;
    var (pm, sk, t, q, pk, rlk) = FHE.KeyGenBGV(n, t0, q0);
    FHE.Show(pm, sk, t, q, pk, rlk);

    var Q = 2 * q;
    var brk = FHE.BRKBGV(pm, sk, t, Q, pk);
    var ak = FHE.AKBGV(pm, sk, t, Q, pk);
    
    for (int k = 0; k < 10; ++k)
    {
        var m1 = FHE.GenUnif(n, t);
        var cm1 = FHE.EncryptBGV(m1, pm, t, q, pk);
        cm1.Show($"ct m1:{m1}");
        
        var ctboot = FHE.Bootstrapping(cm1, pm, t, q, Q, pk, rlk, brk, ak);
        // var ctbootq = ctboot.CoefsMod(q);
        var dec = FHE.DecryptBGV(ctboot, pm, sk, t);
        
        ctboot.Show($"ctboot Q = {Q}");
        Console.WriteLine($"decrypt ctboot:{dec}");
        Console.WriteLine($"m1            :{m1}");
        Console.WriteLine();
        if (!dec.Equals(m1))
            throw new("decrypt");
    }
}

// N = 16 t = 16 q = 2048
// Private Key
// -x^15 + x^14 + x^12 + x^11 - x^9 + x^8 + x^6 - x^5 + x^4 - x^3 - x^2 + x - 1
// Public Key
// 754*x^15 + 1382*x^14 + 1461*x^13 + 1971*x^12 + 111*x^11 + 216*x^10 + 990*x^9 + 189*x^8 + 520*x^7 + 1130*x^6 + 89*x^5 + 1280*x^4 + 1879*x^3 + 1108*x^2 + 1329*x + 1282
// 169*x^15 + 1850*x^14 + 167*x^13 + 1496*x^12 + 1354*x^11 + 979*x^10 + 220*x^9 + 732*x^8 + 205*x^7 + 1328*x^6 + 123*x^5 + 1572*x^4 + 1125*x^3 + 1934*x^2 + 590*x + 7
// Relinearisation Key
// 37*x^15 + 671*x^14 + 1389*x^13 + 651*x^12 + 1841*x^11 + 1543*x^10 + 610*x^9 + 961*x^8 + 674*x^7 + 1739*x^6 + 919*x^5 + 1357*x^4 + 583*x^3 + 1908*x^2 + 1796*x + 523
// 9*x^15 + 497*x^14 + 1932*x^13 + 347*x^12 + 1002*x^11 + 711*x^10 + 330*x^9 + 569*x^8 + 454*x^7 + 527*x^6 + 816*x^5 + 1034*x^4 + 1027*x^3 + 1065*x^2 + 1091*x + 1906
// 
// ct m1:5*x^15 + 14*x^14 + 13*x^13 + 10*x^12 + 13*x^11 + 12*x^10 + 10*x^9 + 2*x^8 + 3*x^7 + 4*x^6 + 9*x^5 + 2*x^4 + x^3 + 14*x^2 + 8*x + 12
// A:685*x^15 + 256*x^14 + 607*x^13 + 477*x^12 + 1047*x^11 + 757*x^10 + 1110*x^9 + 760*x^8 + 716*x^7 + 423*x^6 + 919*x^5 + 971*x^4 + 1374*x^3 + 786*x^2 + 1834*x + 537
// B:570*x^15 + 1666*x^14 + 83*x^13 + 407*x^12 + 582*x^11 + 1885*x^10 + 1749*x^9 + 563*x^8 + 560*x^7 + 744*x^6 + 1364*x^5 + 1530*x^4 + 1321*x^3 + 1472*x^2 + 605*x + 1596
// ctsm
// A:1735*x^15 + 462*x^14 + 429*x^13 + 1504*x^12 + 1182*x^11 + 1713*x^10 + 753*x^9 + 994*x^8 + 3387*x^7 + 1647*x^6 + 473*x^5 + 1892*x^4 + 3162*x^3 + 133*x^2 + 1817*x + 2179
// B:2677*x^15 + 3122*x^14 + 860*x^13 + 921*x^12 + 1673*x^11 + 526*x^10 + 2593*x^9 + 2048*x^8 + 3052*x^7 + 1884*x^6 + 3590*x^5 + 59*x^4 + 3298*x^3 + 617*x^2 + 2617*x + 1507
// ctboot Q = 4096
// A:1780*x^15 + 462*x^14 + 460*x^13 + 1533*x^12 + 1205*x^11 + 1766*x^10 + 775*x^9 + 1050*x^8 + 3399*x^7 + 1686*x^6 + 496*x^5 + 1903*x^4 + 3192*x^3 + 151*x^2 + 1859*x + 2204
// B:2735*x^15 + 3124*x^14 + 879*x^13 + 944*x^12 + 1679*x^11 + 555*x^10 + 2614*x^9 + 2099*x^8 + 3100*x^7 + 1924*x^6 + 3610*x^5 + 117*x^4 + 3339*x^3 + 617*x^2 + 2646*x + 1567
// decrypt ctboot:5*x^15 + 14*x^14 + 13*x^13 + 10*x^12 + 13*x^11 + 12*x^10 + 10*x^9 + 2*x^8 + 3*x^7 + 4*x^6 + 9*x^5 + 2*x^4 + x^3 + 14*x^2 + 8*x + 12
// m1            :5*x^15 + 14*x^14 + 13*x^13 + 10*x^12 + 13*x^11 + 12*x^10 + 10*x^9 + 2*x^8 + 3*x^7 + 4*x^6 + 9*x^5 + 2*x^4 + x^3 + 14*x^2 + 8*x + 12
// 
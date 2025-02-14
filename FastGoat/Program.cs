using System.Numerics;
using FastGoat.Commons;
using FastGoat.Examples;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.LWE;
using FastGoat.UserGroup.Polynoms;
using Rq = FastGoat.Structures.VecSpace.KPoly<FastGoat.UserGroup.Integers.Rational>;
using static FastGoat.Commons.IntExt;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");
// RngSeed(259128);
Ring.DisplayPolynomial = MonomDisplay.StarCaret;
// RecomputeAllPrimesUpTo(5000000);

void testBootstrapping(int N, int t0, int level, int nbTest)
{
    Console.WriteLine("####       Start        ####");
    var (pm, sk, t, primes, sp, pk, rlks) = RLWE.SetupBGV(N, t0, level);
    var n = N / 2;
    var q = primes[0];
    var qL = pk.Q;

    var u = (int)BigInteger.Log2(t.Num);
    var B = new Rational(BigInteger.Pow(2, u));
    var rlk = rlks[qL].rlk;
    var brk = RLWE.BRKgswBGV(sk, pk, B);

    var ak = N.SeqLazy()
        .Select(j => RLWE.SWKBGV(pm, sk, sk.Substitute(pm.X.Pow(j)).ResModSigned(pm, t), t, qL, sp.Pow(level)))
        .ToArray();

    Console.WriteLine($"BGV level = {level}, Gadget Base = {B}");
    Console.WriteLine($"pm = {pm} T = {t} q = {q} sp = {sp} qL = {qL}");
    Console.WriteLine($"sk = {sk}");
    Console.WriteLine($"pk => {pk.Params}");
    Console.WriteLine();

    Console.WriteLine("#### Bootstrapping test ####");
    for (int k = 0; k < nbTest; ++k)
    {
        GlobalStopWatch.AddLap();
        var m = RLWE.GenUnif(n, t);
        var m2 = (m * m).ResModSigned(pm, t);
        var cm = RLWE.EncryptBGV(m, pk); // level qL
        var ct = RLWE.MulRelinBGV(cm, cm, rlk).ModSwitch(q); // level q0

        var (ctboot, ctsm) = RLWE.Bootstrapping(ct, B, pk, ak, brk);
        Console.WriteLine($"ct     {ct.Params}");
        Console.WriteLine($"ctboot {ctboot.Params}");

        var m2boot = RLWE.DecryptBGV(ctboot, sk);
        Console.WriteLine($"m       = {m}");
        Console.WriteLine($"m^2     = {m2}");
        Console.WriteLine($"ctboot  = {m2boot}");
        Console.WriteLine($"eboot   = {RLWE.ErrorsBGV(ctboot, sk).NormInf()}");
        Console.WriteLine($"emodsw  = {RLWE.ErrorsBGV(ct.ModSwitch(qL), sk).NormInf()}");
        Console.WriteLine();
        if (!m2.Equals(m2boot))
            throw new($"step[{k + 1}]");

        var c2 = RLWE.MulRelinBGV(ctboot, ctboot, rlk).ModSwitch(rlks[qL].nextMod);
        var d2 = RLWE.DecryptBGV(c2, sk);
        var m4 = (m2 * m2).ResModSigned(pm, t);
        Console.WriteLine($"m^4 = {m4}");
        Console.WriteLine($"    = {d2}");

        GlobalStopWatch.Show();
        Console.WriteLine();
        if (!d2.Equals(m4))
            throw new($"step[{k + 1}]");
    }

    Console.WriteLine("####        End         ####");
    Console.WriteLine($"BGV level = {level}, Gadget Base = {B}");
    Console.WriteLine($"pm = {pm} T = {t} q = {q} sp = {sp} qL = {qL}");
    Console.WriteLine($"sk = {sk}");
    Console.WriteLine($"pk => {pk.Params}");
    Console.WriteLine();
    GlobalStopWatch.Show();
    Console.WriteLine();
}

{
    // RLWE.NoiseOff();
    RecomputeAllPrimesUpTo(1000000);
    GlobalStopWatch.Restart();

    testBootstrapping(N: 8, t0: 41, level: 4, nbTest: 5);
    testBootstrapping(N: 16, t0: 97, level: 4, nbTest: 5);
    testBootstrapping(N: 32, t0: 97, level: 3, nbTest: 5);
    testBootstrapping(N: 64, t0: 257, level: 3, nbTest: 5);
    testBootstrapping(N: 128, t0: 257, level: 3, nbTest: 1);

    // BGV level = 3, Gadget Base = 64
    // pm = x^16 + 1 T = 97 q = 43457 sp = 90599 qL = 12815095448073162241
    // sk = x^15 - x^14 - x^12 + x^10 - x^9 - x^8 + x^7 - x^5 + x^4 - x^3 - x^2 - 1
    // pk => RLWECipher Q:12815095448073162241    T:97    PM:x^16 + 1
    // 
    // ct     RLWECipher Q:43457    T:97    PM:x^16 + 1
    // ctboot RLWECipher Q:12815095448073162241    T:97    PM:x^16 + 1
    // m       = 10*x^15 + 38*x^14 - 36*x^13 - 21*x^12 - x^11 - 17*x^10 - 21*x^9 - 44*x^8 - 11*x^7 - 36*x^6 - 36*x^5 - 27*x^4 + 34*x^3 + 9*x^2 - 8*x - 16
    // m^2     = -22*x^15 - 26*x^14 - 8*x^13 - 22*x^12 - 22*x^11 - 10*x^10 - 34*x^9 + 28*x^8 + 42*x^7 - 48*x^6 - 41*x^5 + 13*x^4 + 31*x^3 + 32*x^2 - 13*x - 29
    // ctboot  = -22*x^15 - 26*x^14 - 8*x^13 - 22*x^12 - 22*x^11 - 10*x^10 - 34*x^9 + 28*x^8 + 42*x^7 - 48*x^6 - 41*x^5 + 13*x^4 + 31*x^3 + 32*x^2 - 13*x - 29
    // eboot   = 223715756
    // emodsw  = 78146220257711971
    // 
    // m^4 = 23*x^15 - 33*x^14 - 36*x^13 - 12*x^12 - 13*x^11 + 43*x^10 - 9*x^9 + 45*x^8 - 14*x^7 + 10*x^6 - 11*x^5 + 31*x^4 - 36*x^3 + x^2 - 3*x - 19
    //     = 23*x^15 - 33*x^14 - 36*x^13 - 12*x^12 - 13*x^11 + 43*x^10 - 9*x^9 + 45*x^8 - 14*x^7 + 10*x^6 - 11*x^5 + 31*x^4 - 36*x^3 + x^2 - 3*x - 19
    // #  Time:5.799s
    // 
    // ...
    // 
    // BGV level = 3, Gadget Base = 256
    // pm = x^32 + 1 T = 257 q = 82241 sp = 336157 qL = 614815153288590527297
    // sk = x^31 - x^30 + x^29 + x^28 + x^27 - x^26 - x^24 + x^23 - x^22 - x^21 - x^20 + x^19 - x^17 - x^16 - x^15 - x^14 + x^12 + x^11 - x^10 + x^9 - x^7 + x^6 + x^5 - x^4 - x^3 + x^2 - 1
    // pk => RLWECipher Q:614815153288590527297    T:257    PM:x^32 + 1
    // 
    // ct     RLWECipher Q:82241    T:257    PM:x^32 + 1
    // ctboot RLWECipher Q:614815153288590527297    T:257    PM:x^32 + 1
    // m       = 111*x^31 - 116*x^30 - 99*x^29 + 107*x^28 + x^27 + 25*x^26 - 41*x^25 - 3*x^24 - 126*x^23 + 8*x^22 - 32*x^21 - 99*x^20 + 113*x^19 + 45*x^18 + 3*x^17 - 92*x^16 + 9*x^15 - 105*x^14 - 24*x^13 + 59*x^12 + 72*x^11 + 103*x^10 + 69*x^9 + 101*x^8 - 107*x^7 + 113*x^6 + 107*x^5 + 50*x^4 + 86*x^3 + 36*x^2 - 8*x + 58
    // m^2     = -110*x^31 + 9*x^30 - 93*x^29 + 37*x^28 + 29*x^27 - 29*x^26 - 49*x^25 - 97*x^24 - 99*x^23 - 32*x^22 - 23*x^21 - 69*x^20 - 21*x^19 - 84*x^18 - 24*x^17 + 36*x^16 + 21*x^15 + 81*x^14 + 20*x^13 - 27*x^12 + 128*x^11 + 57*x^10 - 46*x^9 + 115*x^8 + 87*x^7 + 67*x^6 + 51*x^5 - 20*x^4 + 116*x^3 - 99*x^2 + 44*x - 30
    // ctboot  = -110*x^31 + 9*x^30 - 93*x^29 + 37*x^28 + 29*x^27 - 29*x^26 - 49*x^25 - 97*x^24 - 99*x^23 - 32*x^22 - 23*x^21 - 69*x^20 - 21*x^19 - 84*x^18 - 24*x^17 + 36*x^16 + 21*x^15 + 81*x^14 + 20*x^13 - 27*x^12 + 128*x^11 + 57*x^10 - 46*x^9 + 115*x^8 + 87*x^7 + 67*x^6 + 51*x^5 - 20*x^4 + 116*x^3 - 99*x^2 + 44*x - 30
    // eboot   = 2352513750
    // emodsw  = 7924320746171689988
    // 
    // m^4 = -107*x^31 - 128*x^30 + 108*x^29 - 37*x^28 + 92*x^27 + 111*x^26 + 87*x^25 + 34*x^24 + 54*x^23 + 26*x^22 + 62*x^21 - 64*x^20 + 7*x^19 + 79*x^18 + 89*x^17 + 23*x^16 - 45*x^15 - 13*x^14 - 48*x^13 - 69*x^12 - 80*x^11 - 105*x^10 - 16*x^9 + 56*x^8 + 36*x^7 + 31*x^6 - 26*x^5 + 79*x^4 - 16*x^3 + 118*x^2 - 71*x - 84
    //     = -107*x^31 - 128*x^30 + 108*x^29 - 37*x^28 + 92*x^27 + 111*x^26 + 87*x^25 + 34*x^24 + 54*x^23 + 26*x^22 + 62*x^21 - 64*x^20 + 7*x^19 + 79*x^18 + 89*x^17 + 23*x^16 - 45*x^15 - 13*x^14 - 48*x^13 - 69*x^12 - 80*x^11 - 105*x^10 - 16*x^9 + 56*x^8 + 36*x^7 + 31*x^6 - 26*x^5 + 79*x^4 - 16*x^3 + 118*x^2 - 71*x - 84
    // #  Time:1m2s
    // 
    // BGV level = 3, Gadget Base = 256
    // pm = x^64 + 1 T = 257 q = 98689 sp = 592643 qL = 4426622753877260980993
    // sk = -x^63 + x^62 - x^61 + x^60 - x^58 - x^57 + x^56 + x^55 - x^54 + x^53 - x^52 - x^51 - x^50 + x^49 + x^48 - x^47 - x^46 + x^45 + x^43 + x^42 - x^41 - x^39 + x^38 - x^37 + x^36 + x^33 + x^32 + x^30 - x^29 + x^28 + x^27 + x^26 - x^25 - x^24 + x^23 + x^22 + x^21 - x^19 + x^18 - x^17 + x^16 - x^15 - x^14 + x^12 - x^11 - x^10 + x^9 + x^8 - x^7 + x^6 - x^5 + x^4 - x^3 + x^2 - x - 1
    // pk => RLWECipher Q:4426622753877260980993    T:257    PM:x^64 + 1
    // 
    // ct     RLWECipher Q:98689    T:257    PM:x^64 + 1
    // ctboot RLWECipher Q:4426622753877260980993    T:257    PM:x^64 + 1
    // m       = -95*x^63 + 102*x^62 + 114*x^61 + 52*x^60 - 92*x^59 - 79*x^58 - 21*x^57 - 66*x^56 + 2*x^55 - 37*x^54 - 26*x^53 - 65*x^52 - 19*x^51 - 95*x^50 + 72*x^49 + 76*x^48 + 59*x^47 - 27*x^46 - 93*x^45 - 78*x^44 + 60*x^43 - 20*x^42 + 125*x^41 + 117*x^40 - 33*x^39 - 72*x^38 + 62*x^37 + 74*x^36 - 56*x^35 + 25*x^34 - 85*x^33 - 36*x^32 - 47*x^31 + 71*x^30 - 93*x^29 + 66*x^28 - 12*x^27 + 114*x^26 + 58*x^25 + 16*x^24 - 9*x^23 + 99*x^22 - 45*x^21 + 88*x^20 + 50*x^19 + 71*x^18 - 94*x^17 + 109*x^16 - 108*x^15 + 119*x^14 + 124*x^13 + 111*x^12 + 22*x^11 + 12*x^10 + 128*x^9 - 102*x^8 + 92*x^7 - 4*x^6 + 9*x^5 + 51*x^4 - 34*x^3 + 71*x^2 - 12*x - 28
    // m^2     = -110*x^63 + 109*x^62 + 80*x^61 - 23*x^60 + 76*x^59 + 120*x^58 - 11*x^57 + 99*x^56 + 100*x^55 - 89*x^54 - 48*x^53 + 58*x^52 - 9*x^51 - 127*x^50 + 122*x^49 + 120*x^48 - 49*x^47 - 104*x^46 - 63*x^45 - 97*x^44 + 104*x^43 + 47*x^42 - 68*x^41 - 58*x^40 + 49*x^39 + 43*x^38 - 25*x^37 + 9*x^36 - 25*x^35 - 55*x^34 + 69*x^33 + 35*x^32 + 113*x^31 + 85*x^30 + 89*x^29 - 68*x^28 + 23*x^27 + 73*x^26 - 55*x^25 + 75*x^24 + 50*x^23 - 87*x^22 - 77*x^21 - 22*x^20 - 23*x^19 - 63*x^18 + 5*x^17 - 54*x^16 + 107*x^15 + 40*x^14 + 73*x^13 + 99*x^12 - 12*x^11 - 23*x^10 - 7*x^9 + 83*x^8 + 59*x^7 + 68*x^6 + 103*x^5 - 84*x^4 - 80*x^3 + 94*x^2 + 102*x - 61
    // ctboot  = -110*x^63 + 109*x^62 + 80*x^61 - 23*x^60 + 76*x^59 + 120*x^58 - 11*x^57 + 99*x^56 + 100*x^55 - 89*x^54 - 48*x^53 + 58*x^52 - 9*x^51 - 127*x^50 + 122*x^49 + 120*x^48 - 49*x^47 - 104*x^46 - 63*x^45 - 97*x^44 + 104*x^43 + 47*x^42 - 68*x^41 - 58*x^40 + 49*x^39 + 43*x^38 - 25*x^37 + 9*x^36 - 25*x^35 - 55*x^34 + 69*x^33 + 35*x^32 + 113*x^31 + 85*x^30 + 89*x^29 - 68*x^28 + 23*x^27 + 73*x^26 - 55*x^25 + 75*x^24 + 50*x^23 - 87*x^22 - 77*x^21 - 22*x^20 - 23*x^19 - 63*x^18 + 5*x^17 - 54*x^16 + 107*x^15 + 40*x^14 + 73*x^13 + 99*x^12 - 12*x^11 - 23*x^10 - 7*x^9 + 83*x^8 + 59*x^7 + 68*x^6 + 103*x^5 - 84*x^4 - 80*x^3 + 94*x^2 + 102*x - 61
    // eboot   = 3832415097
    // emodsw  = 57144336131074694149
    // 
    // m^4 = -101*x^63 - 108*x^62 + 107*x^61 + 22*x^60 - 21*x^59 - 80*x^58 - 9*x^57 - 60*x^56 - 21*x^55 + 3*x^54 - 93*x^53 - 24*x^52 - 68*x^51 - 5*x^50 - 13*x^49 + 92*x^48 - 53*x^47 - 111*x^46 - 12*x^45 + 14*x^44 - 65*x^43 - 102*x^42 - 4*x^41 - 110*x^40 - 107*x^39 - 118*x^38 + 125*x^37 + 109*x^36 + 127*x^35 + 67*x^34 + 24*x^33 + 72*x^32 + 71*x^31 - 91*x^30 - 24*x^29 - 72*x^28 + 96*x^27 + 121*x^25 - 98*x^24 + 119*x^23 - 53*x^22 - 26*x^21 + 18*x^20 - 98*x^19 - 32*x^18 - 46*x^17 + 116*x^16 - 100*x^15 - 23*x^14 + 74*x^13 - 97*x^12 + 36*x^11 - 92*x^10 - 101*x^9 - 47*x^8 - 12*x^7 - 89*x^6 - 83*x^5 - 40*x^4 - 51*x^3 + 117*x^2 + 37*x - 33
    //     = -101*x^63 - 108*x^62 + 107*x^61 + 22*x^60 - 21*x^59 - 80*x^58 - 9*x^57 - 60*x^56 - 21*x^55 + 3*x^54 - 93*x^53 - 24*x^52 - 68*x^51 - 5*x^50 - 13*x^49 + 92*x^48 - 53*x^47 - 111*x^46 - 12*x^45 + 14*x^44 - 65*x^43 - 102*x^42 - 4*x^41 - 110*x^40 - 107*x^39 - 118*x^38 + 125*x^37 + 109*x^36 + 127*x^35 + 67*x^34 + 24*x^33 + 72*x^32 + 71*x^31 - 91*x^30 - 24*x^29 - 72*x^28 + 96*x^27 + 121*x^25 - 98*x^24 + 119*x^23 - 53*x^22 - 26*x^21 + 18*x^20 - 98*x^19 - 32*x^18 - 46*x^17 + 116*x^16 - 100*x^15 - 23*x^14 + 74*x^13 - 97*x^12 + 36*x^11 - 92*x^10 - 101*x^9 - 47*x^8 - 12*x^7 - 89*x^6 - 83*x^5 - 40*x^4 - 51*x^3 + 117*x^2 + 37*x - 33
    // #  Time:13m25s
    // 
    // ####        End         ####
    // BGV level = 3, Gadget Base = 256
    // pm = x^64 + 1 T = 257 q = 98689 sp = 592643 qL = 4426622753877260980993
    // sk = -x^63 + x^62 - x^61 + x^60 - x^58 - x^57 + x^56 + x^55 - x^54 + x^53 - x^52 - x^51 - x^50 + x^49 + x^48 - x^47 - x^46 + x^45 + x^43 + x^42 - x^41 - x^39 + x^38 - x^37 + x^36 + x^33 + x^32 + x^30 - x^29 + x^28 + x^27 + x^26 - x^25 - x^24 + x^23 + x^22 + x^21 - x^19 + x^18 - x^17 + x^16 - x^15 - x^14 + x^12 - x^11 - x^10 + x^9 + x^8 - x^7 + x^6 - x^5 + x^4 - x^3 + x^2 - x - 1
    // pk => RLWECipher Q:4426622753877260980993    T:257    PM:x^64 + 1
    // 
    // #  Time:17m4s
    // 
}
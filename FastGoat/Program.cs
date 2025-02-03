using System.Numerics;
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

void SymbBGV()
{
    var ind = new Indeterminates<Xi>(MonomOrder.GrLex, new Xi("k"), new Xi("T"), new Xi("X"), new Xi("sk"));
    ind.ExtendAppend("epk", "pkb", "erlk", "rlkb");
    ind.ExtendAppend("m1", "u1", "e1a", "e1b");
    ind.ExtendAppend("m2", "u2", "e2a", "e2b");
    ind.ExtendAppend("b");

    var z = new Polynomial<Rational, Xi>(ind, Rational.KZero());
    var Vars = z.Variables;
    var (k, t, X, sk) = Vars.Take(4).Deconstruct();
    var (epk, pkb, erlk, rlkb) = Vars.Skip(4).Take(4).Deconstruct();
    var pka = (t * epk + pkb * sk);
    var b = Vars.Last();
    var p = b * t + 1;
    var rlka = (t * erlk + rlkb * sk + p * sk.Pow(2));

    Console.WriteLine("##### Encryption/Decryption #####");
    Console.WriteLine("Decrypt pk: c0pk - sk * c1pk       = {0}", pka - sk * pkb);
    Console.WriteLine("Decrypt pk: c0pk - sk * c1pk mod t = {0}", (pka - sk * pkb).Div(t).rem);
    Console.WriteLine();

    Console.WriteLine("Decrypt rlk: rlka - sk * rlkb       = {0}", rlka - sk * rlkb);
    Console.WriteLine("Decrypt rlk: rlka - sk * rlkb mod t = {0}", (rlka - sk * rlkb).Div(t).rem);
    Console.WriteLine("                             sk^2   = {0}", sk.Pow(2));
    Console.WriteLine();

    var (m1, u1, e1a, e1b) = Vars.Skip(8).Take(4).Deconstruct();
    var m1a = m1 + t * e1a + u1 * pka;
    var m1b = t * e1b + u1 * pkb;
    Console.WriteLine($"m1a = {m1a}");
    Console.WriteLine($"m1b = {m1b}");
    Console.WriteLine($"Decrypt m1: m1a - sk*m1b       = {m1a - sk * m1b}");
    Console.WriteLine($"Decrypt m1: m1a - sk*m1b mod t = {(m1a - sk * m1b).Div(t).rem}");
    Console.WriteLine();

    var (m2, u2, e2a, e2b) = Vars.Skip(12).Take(4).Deconstruct();
    var m2a = m2 + t * e2a + u2 * pka;
    var m2b = t * e2b + u2 * pkb;
    Console.WriteLine($"m2a = {m2a}");
    Console.WriteLine($"m2b = {m2b}");
    Console.WriteLine($"Decrypt m2: m2a - sk*m1b       = {m2a - sk * m2b}");
    Console.WriteLine($"Decrypt m2: m2a - sk*m2b mod t = {(m2a - sk * m2b).Div(t).rem}");
    Console.WriteLine();

    var add_m1m2 = m1 + m2;
    var add_m1m2a = m1a + m2a;
    var add_m1m2b = m1b + m2b;
    Console.WriteLine("##### Homomorphic Addition of Cipher #####");
    Console.WriteLine($"add_m1m2a = {add_m1m2a}");
    Console.WriteLine($"add_m1m2b = {add_m1m2b}");
    Console.WriteLine($"Decrypt add_m1m2: add_m1m2a - sk*add_m1m2b       = {add_m1m2a - sk * add_m1m2b}");
    Console.WriteLine($"Decrypt add_m1m2: add_m1m2a - sk*add_m1m2b mod t = {(add_m1m2a - sk * add_m1m2b).Div(t).rem}");
    Console.WriteLine($"                                  add_m1m2       = {add_m1m2}");
    Console.WriteLine();

    var mul_km1 = k * m1;
    var mul_km1a = k * m1a;
    var mul_km1b = k * m1b;
    Console.WriteLine("##### Homomorphic Multiplication by Scalar #####");
    Console.WriteLine($"mul_km1a = {mul_km1a}");
    Console.WriteLine($"mul_km1b = {mul_km1b}");
    Console.WriteLine($"Decrypt mul_km1: mul_km1a - sk*mul_km1b       = {mul_km1a - sk * mul_km1b}");
    Console.WriteLine($"Decrypt mul_km1: mul_km1a - sk*mul_km1b mod t = {(mul_km1a - sk * mul_km1b).Div(t).rem}");
    Console.WriteLine($"                                mul_km1       = {mul_km1}");
    Console.WriteLine();

    var mul_m1m2 = m1 * m2;

    var c0 = m1a * m2a;
    var c1 = m1a * m2b + m1b * m2a;
    var c2 = m1b * m2b;

    var mul_m1m2a = p * c0 + c2 * rlka;
    var mul_m1m2b = p * c1 + c2 * rlkb;

    Console.WriteLine("##### Homomorphic Multiplication of Ciphers #####");
    Console.WriteLine($"mul_m1m2a = {mul_m1m2a}");
    Console.WriteLine($"mul_m1m2a mod p = {mul_m1m2a.Div(p).rem}");
    Console.WriteLine($"mul_m1m2b = {mul_m1m2b}");
    Console.WriteLine();
    Console.WriteLine($"Decrypt mul_m1m2: mul_m1m2a - sk*mul_m1m2b       = {mul_m1m2a - sk * mul_m1m2b}");
    Console.WriteLine();
    Console.WriteLine($"Decrypt mul_m1m2: mul_m1m2a - sk*mul_m1m2b mod t = {(mul_m1m2a - sk * mul_m1m2b).Div(t).rem}");
    Console.WriteLine($"                                      mul_m1m2   = {mul_m1m2}");
    Console.WriteLine();
}

void RunLeveledBGV(int N, int level, bool mod2 = false, bool differentPrimes = true)
{
    var t0 = Primes10000.First(t1 => t1 % N == 1);
    if (mod2)
        t0 = 2;

    var sk = RLWE.SKBGV(N / 2);
    var (pm, _, t, primes, pk, rlks) = RLWE.SetupBGV(N, t0, level, sk, differentPrimes);
    Console.WriteLine($"pm = {pm} T = {t} Primes = [{primes.Glue(", ")}]");
    Console.WriteLine($"sk = {sk}");
    Console.WriteLine($"pk => {pk.Params}");

    foreach (var rlk in rlks)
        Console.WriteLine($"rlk[{rlk.Key}] => {rlk.Value.rlk.Params}");

    Console.WriteLine();

    var size = 1 << level;
    var n = pm.Degree;

    var seqMsg = size.SeqLazy().Select(_ => RLWE.GenUnif(n, t)).ToArray();
    var seqCipher = seqMsg.Select(m => RLWE.EncryptBGV(m, pk)).ToArray();
    var mul = seqMsg.Aggregate((xi, xj) => (xi * xj).ResModSigned(pm, t));

    var qMul = new Queue<RLWECipher>(seqCipher);
    while (qMul.Count > 1)
    {
        var c1 = qMul.Dequeue();
        var c2 = qMul.Dequeue();
        var (nextMod, rlk) = rlks[c1.Q];
        var c1c2 = RLWE.MulRelinBGV(c1, c2, rlk).ModSwitch(nextMod);
        qMul.Enqueue(c1c2);
    }

    var cMul = qMul.Dequeue();
    var d_mul = RLWE.DecryptBGV(cMul, sk);

    if (size < 17)
        seqMsg.Println($"level:{level} Size:{size}");
    else
    {
        seqMsg.Take(8).Println($"level:{level} Size:{size}");
        seqMsg.TakeLast(8).Println("    ...");
    }

    Console.WriteLine(" *  {0}", Enumerable.Repeat('-', seqMsg.Max(l => $"{l}".Length)).Glue());
    Console.WriteLine($" =  {mul}");
    Console.WriteLine($"    {d_mul}");
    if (!d_mul.Equals(mul))
        throw new("fail");
}

void testLeveledBGV(int N, int level, bool mod2, bool differentPrimes, int nbTests = 500)
{
    for (int i = 0; i < nbTests; i++)
    {
        Console.WriteLine($"Test[{i + 1}]");
        RunLeveledBGV(N, level, mod2, differentPrimes);
        Console.WriteLine();
    }

    Console.WriteLine($"Pass {nbTests} tests N = {N} Level = {level}");
    Console.WriteLine();
}

void allTests()
{
    RecomputeAllPrimesUpTo(1000000);
    
    testLeveledBGV(N: 2048, level: 1, mod2: true, differentPrimes: true, nbTests: 5);
    testLeveledBGV(N: 2048, level: 1, mod2: true, differentPrimes: false, nbTests: 5);
    
    testLeveledBGV(N: 1024, level: 2, mod2: true, differentPrimes: true, nbTests: 5);
    testLeveledBGV(N: 1024, level: 2, mod2: true, differentPrimes: false, nbTests: 5);
    
    testLeveledBGV(N: 512, level: 3, mod2: true, differentPrimes: true, nbTests: 5);
    testLeveledBGV(N: 512, level: 3, mod2: true, differentPrimes: false, nbTests: 5);
    
    testLeveledBGV(N: 256, level: 5, mod2: true, differentPrimes: true, nbTests: 5);
    testLeveledBGV(N: 256, level: 5, mod2: true, differentPrimes: true, nbTests: 5);
    
    testLeveledBGV(N: 128, level: 5, mod2: true, differentPrimes: true, nbTests: 5);
    testLeveledBGV(N: 128, level: 5, mod2: false, differentPrimes: false, nbTests: 5);
    
    testLeveledBGV(N: 64, level: 5, mod2: false, differentPrimes: true, nbTests: 5);
    testLeveledBGV(N: 64, level: 8, mod2: false, differentPrimes: true, nbTests: 5);
    
    testLeveledBGV(N: 32, level: 4, mod2: true, differentPrimes: true, nbTests: 500);
    testLeveledBGV(N: 32, level: 4, mod2: false, differentPrimes: true, nbTests: 500);
}

{
    allTests();
    SymbBGV();
}
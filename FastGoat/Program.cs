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

void Symb()
{
    var ind = new Indeterminates<Xi>(MonomOrder.GrLex, new Xi("k"), new Xi("T"), new Xi("X"), new Xi("sk"));
    ind.ExtendAppend("epk", "pkb", "erlk", "rlkb");
    ind.ExtendAppend("m1", "u1", "e1a", "e1b");
    ind.ExtendAppend("m2", "u2", "e2a", "e2b");
    ind.ExtendAppend("P");

    var z = new Polynomial<Rational, Xi>(ind, Rational.KZero());
    var Vars = z.Variables;
    var (k, t, X, sk) = Vars.Take(4).Deconstruct();
    var (epk, pkb, erlk, rlkb) = Vars.Skip(4).Take(4).Deconstruct();
    var pka = (t * epk + pkb * sk);
    var p = Vars.Last();
    var rlka = (t * erlk + rlkb * sk + p * sk.Pow(2));

    Console.WriteLine("Decrypt pk: c0pk - sk * c1pk       = {0}", (pka - sk * pkb));
    Console.WriteLine("Decrypt pk: c0pk - sk * c1pk mod t = {0}", (pka - sk * pkb).Div(t).rem);
    Console.WriteLine();

    Console.WriteLine("Decrypt rlk: rlka - sk * rlkb       = {0}", (rlka - sk * rlkb));
    Console.WriteLine("Decrypt rlk: rlka - sk * rlkb mod t = {0}", (rlka - sk * rlkb).Div(t).rem);
    Console.WriteLine("                         p * sk^2   = {0}", (p * sk.Pow(2)));
    Console.WriteLine();

    var (m1, u1, e1a, e1b) = Vars.Skip(8).Take(4).Deconstruct();
    var m1a = (m1 + t * e1a + u1 * pka);
    var m1b = (t * e1b + u1 * pkb);
    Console.WriteLine($"m1a = {m1a}");
    Console.WriteLine($"m1b = {m1b}");
    Console.WriteLine($"Decrypt m1: m1a - sk*m1b       = {(m1a - sk * m1b)}");
    Console.WriteLine($"Decrypt m1: m1a - sk*m1b mod t = {(m1a - sk * m1b).Div(t).rem}");
    Console.WriteLine();

    var (m2, u2, e2a, e2b) = Vars.Skip(12).Take(4).Deconstruct();
    var m2a = (m2 + t * e2a + u2 * pka);
    var m2b = (t * e2b + u2 * pkb);
    Console.WriteLine($"m2a = {m2a}");
    Console.WriteLine($"m2b = {m2b}");
    Console.WriteLine($"Decrypt m2: m2a - sk*m1b       = {(m2a - sk * m2b)}");
    Console.WriteLine($"Decrypt m2: m2a - sk*m2b mod t = {(m2a - sk * m2b).Div(t).rem}");
    Console.WriteLine();

    var add_m1m2 = (m1 + m2);
    var add_m1m2a = m1a + m2a;
    var add_m1m2b = m1b + m2b;
    Console.WriteLine($"add_m1m2a = {add_m1m2a}");
    Console.WriteLine($"add_m1m2b = {add_m1m2b}");
    Console.WriteLine($"Decrypt add_m1m2: add_m1m2a - sk*add_m1m2b       = {(add_m1m2a - sk * add_m1m2b)}");
    Console.WriteLine($"Decrypt add_m1m2: add_m1m2a - sk*add_m1m2b mod t = {(add_m1m2a - sk * add_m1m2b).Div(t).rem}");
    Console.WriteLine($"                                  add_m1m2       = {add_m1m2}");
    Console.WriteLine();

    var mul_km1 = k * m1;
    var mul_km1a = k * m1a;
    var mul_km1b = k * m1b;
    Console.WriteLine($"mul_km1a = {mul_km1a}");
    Console.WriteLine($"mul_km1b = {mul_km1b}");
    Console.WriteLine($"Decrypt mul_km1: mul_km1a - sk*mul_km1b       = {(mul_km1a - sk * mul_km1b)}");
    Console.WriteLine($"Decrypt mul_km1: mul_km1a - sk*mul_km1b mod t = {(mul_km1a - sk * mul_km1b).Div(t).rem}");
    Console.WriteLine($"                                mul_km1       = {mul_km1}");
    Console.WriteLine();

    var mul_m1m2 = (m1 * m2);

    var c0 = m1a * m2a;
    var c1 = m1a * m2b + m1b * m2a;
    var c2 = m1b * m2b;

    var mul_m1m2a = (p * c0 + c2 * rlka);
    var mul_m1m2b = (p * c1 + c2 * rlkb);

    Console.WriteLine($"mul_m1m2a = {mul_m1m2a}");
    Console.WriteLine($"mul_m1m2a mod p = {mul_m1m2a.Div(p).rem}");
    Console.WriteLine($"mul_m1m2b = {mul_m1m2b}");
    Console.WriteLine($"Decrypt mul_m1m2: mul_m1m2a - sk*mul_m1m2b       = {(mul_m1m2a - sk * mul_m1m2b)}");
    Console.WriteLine($"Decrypt mul_m1m2: mul_m1m2a - sk*mul_m1m2b mod t = {(mul_m1m2a - sk * mul_m1m2b).Div(t).rem}");
    Console.WriteLine($"                                  p * mul_m1m2   = {p * mul_m1m2}");
    Console.WriteLine();
}

(Rq pm, Rq sk, Rational t, Rational[] seqMods, Rational sp, RLWECipher pk, Dictionary<Rational, RLWECipher> rlks)
    KeyGenBGV(int n, int t0, int level = 1)
{
    if (!Primes10000.Contains(t0))
        throw new($"T = {t0} must be prime");

    var pm = FG.QPoly().Pow(n) + 1;
    var t = new Rational(t0);
    var arr = DistributionExt.DiceSample(n, [-t.One, t.One]).ToArray();
    var nb = (int)double.Sqrt(n);
    var zeros = DistributionExt.DiceSample(nb, (n - 1).Range()).ToArray();
    foreach (var i in zeros)
        arr[i] *= 0;
    var sk = arr.ToKPoly();

    var sp = new Rational(Primes10000.First(t1 => t1 % (2 * n) == 1 && t1 % t0 == 1));
    var seqMods = level.SeqLazy(1).Select(i => sp.Pow(i)).ToArray();
    var (qL_1, qL) = level == 1 ? (sp, sp) : (seqMods[level - 2], seqMods[level - 1]);

    var epk = RLWE.GenDiscrGauss(n);
    var c1pk = RLWE.GenUnif(n, qL);
    var c0pk = (t * epk + c1pk * sk).ResModSigned(pm, qL);
    var pk = new RLWECipher(c0pk, c1pk, pm, t, qL_1, qL, sp * qL);

    var seqRlks = new Dictionary<Rational, RLWECipher>();
    for (int i = 0; i < level; i++)
    {
        var qi = seqMods[i];
        var erlk = RLWE.GenDiscrGauss(n);
        var c1rlk = RLWE.GenUnif(n, sp.Pow(i) * qi);
        var c0rlk = (t * erlk + c1rlk * sk - sp.Pow(i) * sk.Pow(2)).ResModSigned(pm, sp.Pow(i) * qi);
        var rlk = new RLWECipher(c0rlk, c1rlk, pm, t, qi / sp, sp.Pow(i) * qi, sp.Pow(i + 1) * qi);
        seqRlks[qi] = rlk;
    }

    return (pm, sk, t, seqMods, sp, pk, seqRlks);
}

RLWECipher MulRelinBGV(RLWECipher cipher0, RLWECipher cipher1, RLWECipher rlk)
{
    var (pm, t, q0, q1, q2) = cipher0.PM_T_Q;
    var p1 = rlk.Q1 / q1;
    var f = q1 / rlk.Q2;
    var d0 = (cipher0.A * cipher1.A).ResMod(pm, q1);
    var d1 = (cipher0.A * cipher1.B + cipher0.B * cipher1.A).ResMod(pm, q1);
    var d2 = (cipher0.B * cipher1.B).ResMod(pm, q1);

    var ai = (p1 * d0 - d2 * rlk.A).ResModSigned(pm, rlk.Q2);
    var bi = (p1 * d1 - d2 * rlk.B).ResModSigned(pm, rlk.Q2);
    var af = (ai * f).ClosestModulusTo(ai, t).CoefsModSigned(q0);
    var bf = (bi * f).ClosestModulusTo(bi, t).CoefsModSigned(q0);
    return new RLWECipher(af, bf, pm, t, q0 * q1 / q2, q0, q1);
}

void RunLeveledBGV(int N, int level, bool mod2 = false)
{
    var t0 = Primes10000.First(t1 => t1 % N == 1);
    if (mod2)
        t0 = 2;
    
    var (pm, sk, t, seqMods, sp, pk, rlks) = KeyGenBGV(N / 2, t0, level);
    Console.WriteLine($"pm = {pm} T = {t}");
    Console.WriteLine($"sk = {sk}");
    Console.WriteLine($"seqMods:[{seqMods.Glue(" | ")}] [{sp}]");
    Console.WriteLine($"pk => {pk.Params}");
    foreach (var rlk in rlks)
        Console.WriteLine($"rlk[{rlk.Key}] => {rlk.Value.Params}");

    Console.WriteLine();

    var size = 1 << level - 1;
    var n = pm.Degree;

    var seqMsg = size.SeqLazy().Select(_ => RLWE.GenUnif(n, t)).ToArray();
    var seqCipher = seqMsg.Select(m => RLWE.EncryptBGV(m, pk)).ToArray();
    var mul = seqMsg.Aggregate((xi, xj) => (xi * xj).ResModSigned(pm, t));

    var qMul = new Queue<RLWECipher>(seqCipher);
    while (qMul.Count > 1) qMul.Enqueue(MulRelinBGV(qMul.Dequeue(), qMul.Peek(), rlks[qMul.Dequeue().Q1]));
    var cMul = qMul.Dequeue();
    var d_mul = RLWE.DecryptBGV(cMul, sk);
    
    seqMsg.Println($"level:{level} Size:{size}");
    Console.WriteLine(" *  {0}", Enumerable.Repeat('-', seqMsg.Max(l => $"{l}".Length)).Glue());
    Console.WriteLine($" =  {mul}");
    Console.WriteLine($"    {d_mul}");
    if (!d_mul.Equals(mul))
        throw new("fail");
}

{
    var level = 4;
    for (int i = 0; i < 5; i++)
    {
        try
        {
            RunLeveledBGV(N: 32, level, mod2: true);
            Console.WriteLine();
            RunLeveledBGV(N: 32, level);
            Console.WriteLine();
        }
        catch (Exception e)
        {
            Console.WriteLine(e); // crash times to times
        }
        
        Console.WriteLine();
    }
    
    // Symb();
}

// pm = x^16 + 1 T = 2
// sk = x^15 - x^14 + x^13 - x^12 + x^11 + x^9 - x^7 + x^6 + x^4 - x^2 - x - 1
// seqMods:[97 | 9409 | 912673 | 88529281] [97]
// pk => RLWECipher Q2:8587340257    Q1:88529281    Q0:912673    T:2    PM:x^16 + 1
// rlk[97] => RLWECipher Q2:9409    Q1:97    Q0:1    T:2    PM:x^16 + 1
// rlk[9409] => RLWECipher Q2:88529281    Q1:912673    Q0:97    T:2    PM:x^16 + 1
// rlk[912673] => RLWECipher Q2:832972004929    Q1:8587340257    Q0:9409    T:2    PM:x^16 + 1
// rlk[88529281] => RLWECipher Q2:7837433594376961    Q1:80798284478113    Q0:912673    T:2    PM:x^16 + 1
// 
// level:4 Size:8
//     x^15 + x^13 + x^11 + x^9 + x^8 + x^5 + x^3 + x
//     x^14 + x^13 + x^12 + x^10 + x^8 + x^6 + x^5 + x^4 + x^3 + x^2
//     x^15 + x^12 + x^11 + x^10 + x^7 + x^6 + x^5 + x^4 + x^3 + x^2
//     x^15 + x^12 + x^5 + x^4 + x^3 + x^2 + 1
//     x^15 + x^14 + x^13 + x^12 + x^10 + x^9 + x^8 + x^7 + x^6 + x^5 + x^4 + x^3
//     x^15 + x^11 + x^9 + x^8 + x^7 + x^6 + x^5 + x^4 + x^3 + x + 1
//     x^12 + x^10 + x^9 + x^7 + x^6 + x^5 + x^3
//     x^15 + x^13 + x^10 + x^8 + x^7 + x^4 + x^3 + x + 1
//  *  --------------------------------------------------------------------------
//  =  x^14 + x^13 + x^11 + x^5 + x^4 + x^3 + x^2 + 1
//     x^14 + x^13 + x^11 + x^5 + x^4 + x^3 + x^2 + 1
// 
// pm = x^16 + 1 T = 97
// sk = x^15 + x^13 - x^12 + x^11 - x^9 - x^8 + x^7 + x^5 - x^4 + x^3 - x^2 - x + 1
// seqMods:[43457 | 1888510849 | 82069015964993 | 3566473226790700801] [43457]
// pk => RLWECipher Q2:154988227016643484709057    Q1:3566473226790700801    Q0:82069015964993    T:97    PM:x^16 + 1
// rlk[43457] => RLWECipher Q2:1888510849    Q1:43457    Q0:1    T:97    PM:x^16 + 1
// rlk[1888510849] => RLWECipher Q2:3566473226790700801    Q1:82069015964993    Q0:43457    T:97    PM:x^16 + 1
// rlk[82069015964993] => RLWECipher Q2:6735323381462275915001490049    Q1:154988227016643484709057    Q0:1888510849    T:97    PM:x^16 + 1
// rlk[3566473226790700801] => RLWECipher Q2:12719731277414873549711715808702041601    Q1:292696948188206124438219753059393    Q0:82069015964993    T:97    PM:x^16 + 1
// 
// level:4 Size:8
//     -45*x^15 + 34*x^14 + 19*x^13 + 34*x^12 - 30*x^11 - 39*x^10 - 43*x^9 + 7*x^8 - 27*x^7 - 22*x^6 - 23*x^4 - 38*x^3 - 8*x^2 - x - 21
//     -36*x^15 + 44*x^14 - 21*x^13 + 39*x^12 - 43*x^11 - 35*x^10 - 22*x^9 - 23*x^8 - 5*x^7 - 29*x^6 + 2*x^5 - 36*x^4 - 40*x^3 + 14*x^2 - 20*x - 19
//     -28*x^15 - 16*x^14 - 11*x^13 - 44*x^12 - 35*x^11 + 39*x^10 - 9*x^9 - 15*x^8 - 25*x^7 - 32*x^6 - 40*x^5 - 21*x^4 - 23*x^3 + 15*x^2 - 45*x + 9
//     -41*x^15 + 41*x^14 + 11*x^13 - 46*x^12 - 14*x^11 - 33*x^10 - 25*x^9 + 42*x^8 - 41*x^7 + 17*x^6 - 34*x^5 - 5*x^4 - 40*x^3 - 44*x^2 - 27*x + 45
//     -4*x^15 + 23*x^14 - 8*x^13 + 35*x^12 + 19*x^11 - 42*x^10 - 31*x^9 + 19*x^8 - 35*x^7 - 28*x^6 - 33*x^5 - 31*x^4 - 45*x^3 + 6*x^2 + 10*x + 21
//     -29*x^15 - 17*x^14 + 8*x^13 - 35*x^12 + 45*x^11 - 47*x^10 + 16*x^9 - 41*x^8 + 19*x^7 + 41*x^6 + 11*x^5 - 42*x^4 - 15*x^3 + x^2 + 22*x + 38
//     -40*x^15 - 16*x^14 - 27*x^13 - 42*x^12 - 13*x^11 + 20*x^10 - 44*x^9 + 11*x^8 + 10*x^7 - 14*x^6 - 34*x^5 + 27*x^4 - 10*x^3 + 31*x^2 - 30*x - 10
//     -19*x^15 + 33*x^14 + 39*x^13 + 24*x^12 - 16*x^11 - 42*x^10 - 21*x^9 - 37*x^8 - 38*x^7 - 20*x^6 + 34*x^5 + 15*x^4 + 42*x^3 - 26*x^2 + 17*x - 26
//  *  ----------------------------------------------------------------------------------------------------------------------------------------------
//  =  8*x^15 - 23*x^14 + 26*x^13 - 14*x^12 + 30*x^11 - 34*x^10 + 34*x^9 + 47*x^8 - 35*x^7 + 29*x^6 + 22*x^5 - 41*x^4 - 38*x^3 - 44*x^2 + 22*x - 9
//     8*x^15 - 23*x^14 + 26*x^13 - 14*x^12 + 30*x^11 - 34*x^10 + 34*x^9 + 47*x^8 - 35*x^7 + 29*x^6 + 22*x^5 - 41*x^4 - 38*x^3 - 44*x^2 + 22*x - 9
// 
// 
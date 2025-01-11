using System.Numerics;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Lattice;
using Rq = FastGoat.Structures.VecSpace.KPoly<FastGoat.UserGroup.Integers.Rational>;
using static FastGoat.Commons.IntExt;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");
// RngSeed(57419);
Ring.DisplayPolynomial = MonomDisplay.StarCaret;

Rational ScaleInt(Rational x, Rational q, Rational p, int r)
{
    var m = x.Mod(r);
    var x1 = (x * p / q).Floor;
    return r.SeqLazy().Select(i => x1 + i).First(e => e.Mod(r).Equals(m));
}

Rq ScaleRq(Rq a, Rational q, Rational p, int r) => a.Coefs.Select(x => ScaleInt(x, q, p, r)).ToKPoly();

RLWECipher ScaleRLWECipher(RLWECipher cipher, Rational p, int r)
{
    var a = ScaleRq(cipher.A, cipher.Q, p, r);
    var b = ScaleRq(cipher.B, cipher.Q, p, r);
    return new(a, b, cipher.PM, p);
}

void EncodeDecodeRLWE()
{
    var rlwe = new RLWE(16);
    var (n, t0, q0) = (rlwe.N, rlwe.Q, rlwe.Q * Primes10000.First(q0 => q0 > rlwe.Q));
    var (pm, sk, t, q, pk, rlk) = FHE.KeyGenBGV(n, t0, q0);
    FHE.Show(pm, sk, t, q, pk, rlk);

    for (int k = 0; k < 50; ++k)
    {
        var mi = DistributionExt.DiceSample(n, [0, 1]).ToArray();
        var ec = rlwe.Encode(mi);
        var cipher = FHE.EncryptBGV(ec, pm, t, q, pk);
        var dc = FHE.DecryptBGV(cipher, pm, sk, t);
        var mf = rlwe.Decode(dc);
        Console.WriteLine($"mi:[{mi.Glue(", ")}]");
        Console.WriteLine($"  :{ec}");
        cipher.Show("cipher");
        Console.WriteLine($"  :{dc}");
        Console.WriteLine($"mf:[{mf.Glue(", ")}]");
        Console.WriteLine();
        if (!mi.SequenceEqual(mf))
            throw new();
    }
}

RLWECipher NOTBGV(RLWECipher cipher, Rq pm, Rational t, Rational q, RLWECipher rlk, RLWECipher[] encXpow,
    RLWECipher[] exsk)
{
    var t0 = (int)t.Num;
    var th = t0 / 2;
    var thi = t0 - th;
    var z = InvModPbez(thi, t0);
    var f = z * th % t0;
    var cx = FHE.ExtractCoefs(cipher, pm, exsk);
    return f * cx.Zip(encXpow).Select(e => FHE.MulRelinBGV(e.First + thi, e.Second, pm, q, rlk))
        .Aggregate((a0, a1) => a0 + a1);
}

RLWECipher ANDBGV(RLWECipher c1, RLWECipher c2, Rq pm, Rational t, Rational q, RLWECipher rlk, RLWECipher[] encXpow,
    RLWECipher[] exsk)
{
    var t0 = (int)t.Num;
    var th = t0 / 2;
    var err = th * th % t0;
    var z = InvModPbez(err, t0);
    var f = z * th % t0;
    var c1x = FHE.ExtractCoefs(c1, pm, exsk);
    var c2x = FHE.ExtractCoefs(c2, pm, exsk);
    return f * c1x.Zip(c2x).Select(e => FHE.MulRelinBGV(e.First, e.Second, pm, q, rlk))
        .Zip(encXpow).Select(e => FHE.MulRelinBGV(e.First, e.Second, pm, q, rlk))
        .Aggregate((a0, a1) => a0 + a1);
}

RLWECipher NANDBGV(RLWECipher c1, RLWECipher c2, Rq pm, Rational t, Rational q, RLWECipher rlk, RLWECipher[] encXpow,
    RLWECipher[] exsk)
{
    return NOTBGV(ANDBGV(c1, c2, pm, t, q, rlk, encXpow, exsk), pm, t, q, rlk, encXpow, exsk);
}

RLWECipher NORBGV(RLWECipher c1, RLWECipher c2, Rq pm, Rational t, Rational q, RLWECipher rlk, RLWECipher[] encXpow,
    RLWECipher[] exsk)
{
    var nc1 = NOTBGV(c1, pm, t, q, rlk, encXpow, exsk);
    var nc2 = NOTBGV(c2, pm, t, q, rlk, encXpow, exsk);
    return ANDBGV(nc1, nc2, pm, t, q, rlk, encXpow, exsk);
}

RLWECipher ORBGV(RLWECipher c1, RLWECipher c2, Rq pm, Rational t, Rational q, RLWECipher rlk, RLWECipher[] encXpow,
    RLWECipher[] exsk)
{
    return NOTBGV(NORBGV(c1, c2, pm, t, q, rlk, encXpow, exsk), pm, t, q, rlk, encXpow, exsk);
}

RLWECipher XORBGV(RLWECipher c1, RLWECipher c2, Rq pm, Rational t, Rational q, RLWECipher rlk, RLWECipher[] encXpow,
    RLWECipher[] exsk)
{
    var c1_and_nc2 = ANDBGV(c1, NOTBGV(c2, pm, t, q, rlk, encXpow, exsk), pm, t, q, rlk, encXpow, exsk);
    var nc1_and_c2 = ANDBGV(NOTBGV(c1, pm, t, q, rlk, encXpow, exsk), c2, pm, t, q, rlk, encXpow, exsk);
    return ORBGV(c1_and_nc2, nc1_and_c2, pm, t, q, rlk, encXpow, exsk);
}

void ExpandExtract()
{
    var (n, t0, q0) = (8, 17, 289);
    var (pm0, sk0, t, q, pk, rlk) = FHE.KeyGenBGV(n, t0, q0);
    var m = FHE.GenUnif(n, t);
    var c1 = FHE.EncryptBGV(m, pm0, t, q, pk);

    var ind = new Indeterminates<Xi>(MonomOrder.GrLex, new Xi("X"));
    ind.ExtendAppend(n.SeqLazy().Select(i => new Xi($"s{i}")).ToArray());
    ind.ExtendAppend(n.SeqLazy().Select(i => new Xi($"a{i}")).ToArray());
    ind.ExtendAppend(n.SeqLazy().Select(i => new Xi($"b{i}")).ToArray());
    var xis = new Polynomial<Rational, Xi>(ind, Rational.KZero());
    var x = xis.Variables.First();
    var sk = xis.Variables.Skip(1).Take(n).ToVec();
    var ak = xis.Variables.Skip(n + 1).Take(n).ToVec();
    var bk = xis.Variables.Skip(2 * n + 1).Take(n).ToVec();
    var pm = x.Pow(n) + 1;

    var a = ak.Select((ai, i) => ai * x.Pow(i)).Aggregate((ai, aj) => ai + aj);
    var b = bk.Select((bi, i) => bi * x.Pow(i)).Aggregate((bi, bj) => bi + bj);
    var s = sk.Select((si, i) => si * x.Pow(i)).Aggregate((si, sj) => si + sj);

    var prod = (a - b * s).Div(pm).rem;
    var decompProd = Ring.Decompose(prod, "X");

    Console.WriteLine($"s = {sk0}");
    c1.Show($"m = {m}");
    Console.WriteLine();

    Console.WriteLine($"A = {a}");
    Console.WriteLine($"B = {b}");
    Console.WriteLine($"s = {s}");
    Console.WriteLine($"A - s * B = {prod}");
    Console.WriteLine($"A - s * B = {decompProd.Item1.Select(l => $"({l.Value}) * {l.Key}").Glue(" +\n      ")}");
    Console.WriteLine();

    var extr = FHE.Extract(c1, n);

    Console.WriteLine($"{ak.Select((ai, i) => $"{ai} = {c1.A[i],4}").Glue(", ")}");
    Console.WriteLine($"{bk.Select((bi, i) => $"{bi} = {c1.B[i],4}").Glue(", ")}");
    Console.WriteLine($"{sk.Select((si, i) => $"{si} = {sk0[i],4}").Glue(", ")}");
    extr.Println(l => $"<{-l.bi.ToVec()}, s> + {l.ai}", "Extract");
    extr.Select(e => e.ai - e.bi.Zip(sk).Select(f => f.First * f.Second).Aggregate((e1, e2) => e1 + e2))
        .Println("Extract substitute");

    Console.WriteLine($"m = {m}");
    var m2 = extr.Select(e => e.ai - e.bi.Select((bi, i) => bi * sk0[i]).Aggregate((e1, e2) => e1 + e2))
        .Select((e, i) => e.Mod(t) * x.Pow(i)).Aggregate((e1, e2) => e1 + e2);
    Console.WriteLine($"  = {m2}");
}

void TestExtract()
{
    var (n, t0, q0) = (8, 17, 289);
    var (pm, sk, t, q, pk, rlk) = FHE.KeyGenBGV(n, t0, q0);
    FHE.Show(pm, sk, t, q, pk, rlk);
    var (_, esk) = FHE.EXSK(pm, sk, t, q, pk);

    for (int k = 0; k < 50; ++k)
    {
        var m = FHE.GenUnif(n, t);
        var c1 = FHE.EncryptBGV(m, pm, t, q, pk);
        var ciphers = FHE.ExtractCoefs(c1, pm, esk);
        foreach (var (ci, i) in ciphers.Select((ci, i) => (ci, i)))
        {
            var di = FHE.DecryptBGV(ci, pm, sk, t);
            if (!di.Equals(di.One * m[i]))
                throw new($"step[{k}]; m{i} = {m[i]} != {di}");
        }

        var seqDecrypt = ciphers.Select(ci => FHE.DecryptBGV(ci, pm, sk, t)).ToArray();
        Console.WriteLine($"m:{m}");
        Console.WriteLine($" :{seqDecrypt.Reverse().ToVec()}");
        Console.WriteLine();
    }
}

void BGVNot(int n, int nbTrials = 50)
{
    var rlwe = new RLWE(n);
    var (_, t0, q0) = (rlwe.N, rlwe.Q, rlwe.Q * Primes10000.First(q0 => q0 > rlwe.Q));
    var (pm, sk, t, q, pk, rlk) = FHE.KeyGenBGV(n, t0, q0);
    FHE.Show(pm, sk, t, q, pk, rlk);
    var (encXpow, exsk) = FHE.EXSK(pm, sk, t, q, pk);

    for (int k = 0; k < nbTrials; ++k)
    {
        var m1 = DistributionExt.DiceSample(n, [0, 1]).ToArray();
        var ec1 = rlwe.Encode(m1);
        var c1 = FHE.EncryptBGV(ec1, pm, t, q, pk);

        var nmi = m1.Select(e => 1 - e).ToArray();
        var nc1 = NOTBGV(c1, pm, t, q, rlk, encXpow, exsk);

        var dnm = FHE.DecryptBGV(nc1, pm, sk, t);
        var nmf = rlwe.Decode(dnm);

        Console.WriteLine($"ec1 :{ec1}");
        Console.WriteLine($"    :{dnm}");
        Console.WriteLine($"m   :[{m1.Glue()}]");
        Console.WriteLine($"!m  :[{nmi.Glue()}]");
        Console.WriteLine($"    :[{nmf.Glue()}]");
        if (!nmi.SequenceEqual(nmf))
            throw new($"step[{k}]");

        Console.WriteLine();
    }
}

int OpBool(string name, int c1, int c2)
{
    if (name == "and")
        return c1 & c2;
    if (name == "nand")
        return 1 - (c1 & c2);
    if (name == "nor")
        return 1 - (c1 | c2);
    if (name == "or")
        return c1 | c2;
    else
        return c1 ^ c2;
}

RLWECipher OpBGV(string name, RLWECipher c1, RLWECipher c2, Rq pm, Rational t, Rational q, RLWECipher rlk,
    RLWECipher[] encXpow, RLWECipher[] exsk)
{
    if (name == "and")
        return ANDBGV(c1, c2, pm, t, q, rlk, encXpow, exsk);
    if (name == "nand")
        return NANDBGV(c1, c2, pm, t, q, rlk, encXpow, exsk);
    if (name == "nor")
        return NORBGV(c1, c2, pm, t, q, rlk, encXpow, exsk);
    if (name == "or")
        return ORBGV(c1, c2, pm, t, q, rlk, encXpow, exsk);
    else
        return XORBGV(c1, c2, pm, t, q, rlk, encXpow, exsk);
}

void BGVGates(int n, string name, int nbTrials = 50)
{
    var rlwe = new RLWE(n);
    var (_, t0, q0) = (rlwe.N, rlwe.Q, rlwe.Q * Primes10000.First(q0 => q0 > rlwe.Q));
    var (pm, sk, t, q, pk, rlk) = FHE.KeyGenBGV(n, t0, q0);
    FHE.Show(pm, sk, t, q, pk, rlk);
    var (encXpow, exsk) = FHE.EXSK(pm, sk, t, q, pk);

    for (int k = 0; k < nbTrials; ++k)
    {
        var m1 = DistributionExt.DiceSample(n, [0, 1]).ToArray();
        var ec1 = rlwe.Encode(m1);
        var c1 = FHE.EncryptBGV(ec1, pm, t, q, pk);
        var m2 = DistributionExt.DiceSample(n, [0, 1]).ToArray();
        var ec2 = rlwe.Encode(m2);
        var c2 = FHE.EncryptBGV(ec2, pm, t, q, pk);

        var m1m2i = m1.Zip(m2).Select(e => OpBool(name, e.First, e.Second)).ToArray();
        var c1c2 = OpBGV(name, c1, c2, pm, t, q, rlk, encXpow, exsk);

        var d1d2 = FHE.DecryptBGV(c1c2, pm, sk, t);
        var m1m2f = rlwe.Decode(d1d2);

        var fmt = $"{name}(m1,m2)";
        Console.WriteLine($"ec1        :{ec1}");
        Console.WriteLine($"ec2        :{ec2}");
        Console.WriteLine($"           :{d1d2}");
        Console.WriteLine($"m1         :[{m1.Glue()}]");
        Console.WriteLine($"m2         :[{m2.Glue()}]");
        Console.WriteLine($"{fmt,-11}:[{m1m2i.Glue()}]");
        Console.WriteLine($"           :[{m1m2f.Glue()}]");
        if (!m1m2i.SequenceEqual(m1m2f))
            throw new($"step[{k}]");

        Console.WriteLine();
    }
}

{
    // EncodeDecodeRLWE();
    // ExpandExtract();
    
    var n = 32;
    BGVNot(n);
    BGVGates(n, "and");
    BGVGates(n, "nand");
    BGVGates(n, "nor");
    BGVGates(n, "or");
    BGVGates(n, "xor");
}
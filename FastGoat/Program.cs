using System.Numerics;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Floats;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Polynoms;
using static FastGoat.Commons.IntExt;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");

void testPSQL1()
{
    var n = 8;
    var O1 = 25;
    var O2 = O1 + n;
    var pi = BigReal.Pi(O2);
    var beta = BigReal.FromBigIntegerAndExponent(BigInteger.Parse("-1669947371922907049619"), 1, O1);

    var ai = n.Range().Select(k => pi.Pow(k).ToBigReal(O1)).ToKMatrix();
    ai.Coefs[0, n - 1] = beta;

    // var y = 2 / BigCplx.Sqrt(BigCplx.FromBigInteger(3, O2));
    var y = 3 * pi.One / 2;
    // var y = BigCplx.FromBigInteger(2, O2);
    Console.WriteLine(ai);
    Console.WriteLine(y);
    GlobalStopWatch.AddLap();
    var coefs = PSLQM2<Dble, BigReal>.TwoLevelMultipair(ai, y, n - 1).Select(e => e.ToRational).ToArray();
    GlobalStopWatch.Show($"PSLQ");
    var P = -FG.KPoly('π', coefs.SkipLast(1).ToArray()) / coefs.Last();
    Console.WriteLine($"beta = {P}");
    Console.WriteLine($"diff = {P.Substitute(pi).ToBigReal(O1) - beta}");
    Console.WriteLine();
}

void PSLQminPoly(int r, int s, int O)
{
    var n = r * s + 1; // Expected polynomial degree plus one
    var o = BigReal.BrOne(O);
    var alpha = BigReal.NthRoot(3 * o, r) - BigReal.NthRoot(2 * o, s); // a = 3^(1/r) - 2^(1/s)
    var ai = n.Range().Select(k => alpha.Pow(k)).ToKMatrix();

    var O2 = O + n;
    var y = 2 / BigReal.Sqrt(BigReal.FromBigInteger(3, O2));
    // var y = 3 * BigReal.BrOne(O2) / 2;
    // var y = BigReal.FromBigInteger(2, O2);

    GlobalStopWatch.AddLap();
    var coefs = PSLQM2<Dble, BigReal>.TwoLevelMultipair(ai, y, n - 1).Select(e => e.ToRational).ToArray();
    GlobalStopWatch.Show($"PSLQ min poly a = 3^(1/{r}) - 2^(1/{s})");
    var P = FG.KPoly('X', coefs).Monic;
    Console.WriteLine($"P = {P}");
    Console.WriteLine($"P(a) = {P.Substitute(alpha).ToBigReal(O / 2)}");
    Console.WriteLine();
    
    GlobalStopWatch.AddLap();
    var x = P.X;
    Console.WriteLine("From algebraic method");
    Console.WriteLine($"P = {IntFactorisation.PrimitiveElt(x.Pow(r) - 3, x.Pow(s) - 2)[0].F.SubstituteChar('X')}");
    GlobalStopWatch.Show();
    Console.WriteLine();
}

(BigCplx[], BigCplx sum) AlphaBetaPolynomial(BigCplx alpha, BigCplx beta, int d)
{
    if (Logger.Level != LogLevel.Off)
        Console.WriteLine("Start LLL algorithm");

    var z = BigCplx.BcZero(alpha.O);
    var mat = new KMatrix<BigCplx>(z.Zero, 1, d).Zero;
    var ai = alpha.One;
    var aipow = new List<BigCplx>();
    for (int i = 0; i < mat.N - 1; i++)
    {
        aipow.Add(ai);
        mat.Coefs[0, i] = ai;
        ai *= alpha;
    }

    mat.Coefs[0, mat.N - 1] = beta;

    var O = alpha.O;
    var O2 = O + 2 * d;
    // var y = 2 / BigCplx.Sqrt(3 * BigCplx.BcOne(O2));
    // var y = 3 * BigCplx.BcOne(O2) / 2;
    var y = BigCplx.Sqrt(2 * BigCplx.BcOne(O2));
    var rel = PSLQM2<Cplx, BigCplx>.TwoLevelMultipair(mat, y, d - 1);
    var sum = rel.SkipLast(1).Zip(aipow).Aggregate(z.Zero, (acc, v) => acc + v.First * v.Second) / beta;
    return (rel.SkipLast(1).Select(c => c.RoundEven).ToArray(), sum.RoundEven);
}

void ExamplePol(KPoly<Rational> P, int O)
{
    GlobalStopWatch.AddLap();
    Ring.DisplayPolynomial = MonomDisplay.StarCaret;
    Console.WriteLine(P);
    var (X, a) = FG.EPolyXc(P, 'a');
    var l = IntFactorisation.AlgebraicRoots(X.Pow(2) + 1);
    var nRoots = FG.NRoots(P.ToBcPoly(O));
    nRoots.Println("Roots");
    var alpha = nRoots[0];
    Console.WriteLine($"alpha={alpha}");

    var Qi = l.Count != 0;
    Console.WriteLine($"Qi {Qi}");
    var I = Qi ? l.First(e => e.Poly.Substitute(alpha).Equals(alpha.I)) : a.Zero;
    Console.WriteLine($"{I} {I * I} {I.Poly.Substitute(alpha)}");

    EPoly<Rational> Algebraic(BigCplx e) => e.RealPart.ToRational + e.ImaginaryPart.ToRational * I;
    var set = new HashSet<EPoly<Rational>>() { a };
    var If = alpha.I;
    var eqCplx = EqualityComparer<BigCplx>.Create(
        (a0, a1) => (a0 - a1).IsZero4d(), a0 => a0.O
    );
    var remRoots = nRoots.Skip(1).ToHashSet(eqCplx);
    var dicRoots = nRoots.ToDictionary(e => e.ToBigCplx(O / 2), e => e, eqCplx);
    var alphaOpp = (-alpha).ToBigCplx(O / 2);
    if (dicRoots.ContainsKey(alphaOpp))
    {
        remRoots.Remove(dicRoots[alphaOpp]);
        set.Add(-a);
    }

    while (true)
    {
        Console.WriteLine($"remRoots:{remRoots.Count}");
        if (remRoots.Count == 0)
            break;

        var beta = remRoots.First();
        var (col, sum) = AlphaBetaPolynomial(alpha, beta, P.Degree + 1);
        if (!Qi && sum.RealPart.IsZero())
        {
            sum /= If;
            col = col.Select(e => e / If).ToArray();
        }

        var sumAlg = Algebraic(sum);
        var relAlg = col.Select((c, i) => a.Pow(i) * Algebraic(c)).Aggregate((c0, c1) => c0 + c1);
        var betaAlg = relAlg / sumAlg;
        var Pbeta = P.Substitute(betaAlg);
        Console.WriteLine($"[{col.Glue(",")}]");
        Console.WriteLine($"sum = {sum.RoundEven}");
        Console.WriteLine($"beta={betaAlg}");
        Console.WriteLine($"P(beta)={Pbeta}");
        Console.WriteLine();

        if (Pbeta.IsZero())
        {
            set.UnionWith(Group.KAut(set.Append(betaAlg).ToArray()).Select(e => e.E));
            remRoots.ExceptWith(set.Select(e => e.Poly.Substitute(alpha).ToBigCplx(O / 2))
                .Select(e => dicRoots[e]));
        }
    }

    var P2 = set.Select(e => X - e).Aggregate((x0, x1) => x0 * x1);
    if (!P2.Equals(P.Substitute(X)))
        throw new();

    DisplayGroup.HeadElementsNames(Group.KAut(set.ToArray()));
    GlobalStopWatch.Show($"P = {P}");
    Console.WriteLine();
}

void StartingTest()
{
    for (int k = 0; k < 10; ++k)
        testPSQL1();

    // MinPoly of a = 3^(1/r) - 2^(1/s)
    PSLQminPoly(r: 2, s: 2, O: 20);
    PSLQminPoly(r: 3, s: 3, O: 40);
    PSLQminPoly(r: 2, s: 5, O: 50);
    PSLQminPoly(r: 3, s: 4, O: 70);
    PSLQminPoly(r: 2, s: 7, O: 90);
    PSLQminPoly(r: 3, s: 5, O: 90);
    PSLQminPoly(r: 4, s: 4, O: 90);
    PSLQminPoly(r: 4, s: 5, O: 120);
    PSLQminPoly(r: 5, s: 5, O: 180);
    PSLQminPoly(r: 5, s: 6, O: 250);
    PSLQminPoly(r: 6, s: 6, O: 320);
}

void ExamplesGaloisPolynomial()
{
    var x = FG.QPoly();
    ExamplePol(x.Pow(4) - 2 * x.Pow(2) + 9, 20); // C2 x C2
    ExamplePol(x.Pow(4) + x.Pow(3) + x.Pow(2) + x + 1, 20); // C4
    
    ExamplePol(x.Pow(6) + 2 * x.Pow(5) - 2 * x.Pow(3) + 7 * x.Pow(2) + 8 * x + 13, 20); // C6
    ExamplePol(x.Pow(6) + 108, 30); // S3
    
    ExamplePol(x.Pow(4) - 8 * x.Pow(3) + 20 * x.Pow(2) - 16 * x + 2, 30); // C4
    ExamplePol(x.Pow(8) - 8 * x.Pow(6) + 20 * x.Pow(4) - 16 * x.Pow(2) + 2, 30); // C8
    ExamplePol(x.Pow(8) - x.Pow(4) + 1, 30); // C2 x C2 x C2
    ExamplePol(x.Pow(8) + 1, 30); // C2 x C4
    ExamplePol(x.Pow(8) + 4 * x.Pow(6) + 2 * x.Pow(4) + 28 * x.Pow(2) + 1, 30); // D8
    ExamplePol(x.Pow(8) - 12 * x.Pow(6) + 36 * x.Pow(4) - 36 * x.Pow(2) + 9, 30); // Q8;
    
    ExamplePol(x.Pow(10) - 2 * x.Pow(9) - 2 * x.Pow(8) + 6 * x.Pow(7) + 14 * x.Pow(6) - 26 * x.Pow(5) + 15 * x.Pow(4) -
        4 * x.Pow(3) + 58 * x.Pow(2) - 40 * x + 89, 60); // C10
    
    ExamplePol(x.Pow(10) + 10 * x.Pow(8) + 125 * x.Pow(6) + 500 * x.Pow(4) + 2500 * x.Pow(2) + 4000, 100); // D10
    
    ExamplePol(x.Pow(12) + 6 * x.Pow(8) + 26 * x.Pow(6) - 63 * x.Pow(4) + 162 * x.Pow(2) + 81, 100); // A4
    
    ExamplePol(x.Pow(18) + 171 * x.Pow(12) + 5130 * x.Pow(6) + 27, 120); // C3 x: C6
    
    ExamplePol(x.Pow(20) + 2500 * x.Pow(10) + 50000, 130); // C5 x: C4
}

void LQtest()
{
    var I = Cplx.I;
    for (int k = 0; k < 100; ++k)
        foreach (var n in 5.Range(2))
        {
            var A = (n * n).Range().Select(_ => Rng.Next(-9, 10) + I * Rng.Next(-9, 10)).ToKMatrix(n);
            var detA = A.Det;
            if (detA.IsZero())
                continue;
            
            var L = PSLQ.LQ(A);
            var Q = L.Inv() * A;
            Console.WriteLine("A");
            Console.WriteLine(A);
            Console.WriteLine("L");
            Console.WriteLine(L);
            Console.WriteLine("Q");
            Console.WriteLine(Q);
            var QQc = Q * Q.Conj();
            var QQc0 = QQc.Select(e => e.RoundEven).ToKMatrix(n);
            var check1 = QQc0.Equals(QQc0.One);
            var detL = L.Det;
            var c = detA / detL;
            var check2 = (c * c.Conj - 1).Absolute2.IsZero();
            var check3 = (c + (-1).Pow(n)).Absolute2.IsZero();
            Console.WriteLine($"Q x Qconj = I{n} {QQc0.Equals(QQc0.One)}");
            Console.WriteLine(QQc);
            Console.WriteLine($"Det A = {detA}");
            Console.WriteLine($"Det L = {detL} = {c} * Det A");
            Console.WriteLine();
            if (!check1 || !check2)
                throw new($"check1:{check1} check2:{check2} check3:{check3}");
        }
}

{
    GlobalStopWatch.Restart();
    Ring.DisplayPolynomial = MonomDisplay.StarCaret;
    
    LQtest();
    StartingTest();
    ExamplesGaloisPolynomial();
    
    Console.Beep();
}
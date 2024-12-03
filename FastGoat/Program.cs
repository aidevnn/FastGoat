using System.CodeDom;
using System.Collections;
using System.Diagnostics;
using System.Numerics;
using System.Reflection.Emit;
using FastGoat.Commons;
using FastGoat.Examples;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.Naming;
using FastGoat.Structures.Subgroups;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Characters;
using FastGoat.UserGroup.DatabaseSmallGroups;
using FastGoat.UserGroup.EllCurve;
using FastGoat.UserGroup.Floats;
using FastGoat.UserGroup.GModuleN;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Matrix;
using FastGoat.UserGroup.GModuleN;
using FastGoat.UserGroup.Perms;
using FastGoat.UserGroup.Polynoms;
using FastGoat.UserGroup.Words;
using static FastGoat.Commons.IntExt;
using static FastGoat.Commons.EnumerableExt;
using FastGoat.UserGroup.Padic;
using FastGoat.UserGroup.Words.Tools;
using GroupRegX = System.Text.RegularExpressions;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");

void IterOneLevelMultipair(KMatrix<BigCplx> H, KMatrix<BigCplx> A, KMatrix<BigCplx> B,
    KMatrix<BigCplx> T, KMatrix<BigCplx> y, (int i, BigCplx yi)[] gamma_pow, bool imq = false)
{
    // 1. Sort the entries of the (n − 1)-long vector {γ^i |Hii |} in decreasing order,
    // producing the sort indices.
    // 2. Beginning at the sort index m1 corresponding to the largest γ^i |Hii |, select
    // pairs of indices (mi , mi + 1), where mi is the sort index. If at any step
    // either mi or mi + 1 has already been selected, pass to the next index in
    // the list. Continue until either βn pairs have been selected, or the list
    // is exhausted. Let p denote the number of pairs actually selected in this
    // manner.
    // 3. For i := 1 to p, exchange the entries of y indexed mi and mi + 1, and the
    // corresponding rows of A, B and H; endfor.
    var m = gamma_pow.Select(e => (e.i, (e.yi * H[e.i, e.i].Absolute)))
        .OrderByDescending(e => e.Item2).ToArray();
    var n = H.M;
    var list = new List<(int, int)>(n);
    var selected = new HashSet<int>(n);
    for (int i = 0; i < n - 1; i++)
    {
        var (b, c) = (m[i].i, m[i].i + 1);
        if (selected.Contains(b) || selected.Contains(c))
            continue;

        list.Add((b, c));
        selected.UnionWith([b, c]);
        if (imq || list.Count + 1 > 0.4 * n)
            break;
    }

    foreach (var (mi, _) in list)
    {
        // Console.WriteLine($"Swap {(mi, mi + 1)}");
        (y.Coefs[mi, 0], y.Coefs[mi + 1, 0]) = (y[mi + 1, 0], y[mi, 0]);
        for (int k = 0; k < n; k++)
        {
            (A.Coefs[mi, k], A.Coefs[mi + 1, k]) = (A.Coefs[mi + 1, k], A.Coefs[mi, k]);
            (B.Coefs[mi, k], B.Coefs[mi + 1, k]) = (B.Coefs[mi + 1, k], B.Coefs[mi, k]);
            if (k < n - 1)
                (H.Coefs[mi, k], H.Coefs[mi + 1, k]) = (H.Coefs[mi + 1, k], H.Coefs[mi, k]);
        }
    }

    // 4. Remove corners on H diagonal: For j := 1 to p: if mj ≤ n − 2 then set
    // t0 := Sqrt(H2mj_mj + H2mj_mj1), t1 := Hmj_mj /t0 and t2 := Hmj_mj1 /t0 ; for 
    // i := mj to n: set t3 := Hi_mj ; t4 := Hi_mj1 ; Hi_mj := t*1t3 + t*2t4 ; and
    // Hi_mj1 := −t2t3 + t1t4 ; endfor; endif; endfor.
    for (int j = 0; j < list.Count; j++)
    {
        var mj = list[j].Item1;
        if (mj < n - 2)
        {
            var (b, c) = (H[mj, mj], H[mj, mj + 1]);
            var t0 = BigCplx.Sqrt(b.Absolute2 + c.Absolute2);
            var t1 = b / t0;
            var t2 = c / t0;
            var t1c = t1.Conj;
            var t2c = t2.Conj;
            for (int i = mj; i < n; i++)
            {
                var t3 = H[i, mj];
                var t4 = H[i, mj + 1];
                H.Coefs[i, mj] = t1c * t3 + t2c * t4;
                H.Coefs[i, mj + 1] = -t2 * t3 + t1 * t4;
            }
        }
    }

    // 5. Reduce H: For i := 2 to n: for j := 1 to n − i + 1: set l := i + j − 1; for
    // k := j +1 to l −1: set Hlj := Hlj −Tlk Hkj ; endfor; set Tlj := nint(Hlj /Hjj )
    // and Hlj := Hlj − Tlj Hjj ; endfor; endfor. [Note that the n x (n − 1) integer
    // array T is set before it is used.]
    for (int i = 1; i < n; i++)
    {
        for (int j = 0; j < n - i; j++)
        {
            var l = i + j;
            for (int k = j + 1; k < l; k++)
                H.Coefs[l, j] -= T[l, k] * H[k, j];

            T.Coefs[l, j] = (H[l, j] / H[j, j]).RoundEven;
            H.Coefs[l, j] -= T[l, j] * H[j, j];
        }
    }

    // 6. Update y: For j := 1 to n − 1: for i := j + 1 to n: set yj := yj + Tij yi ;
    // endfor; endfor.
    // 7. Update A and B: For k := 1 to n: for j := 1 to n − 1: for i := j + 1 to n:
    // set Aik := Aik − Tij Ajk and Bjk := Bjk + Tij Bik ; endfor; endfor; endfor.
    for (int j = 0; j < n - 1; j++)
    {
        for (int i = j + 1; i < n; i++)
            y.Coefs[j, 0] += T[i, j] * y[i, 0];
    }

    for (int k = 0; k < n; k++)
    {
        for (int j = 0; j < n - 1; j++)
        {
            for (int i = j + 1; i < n; i++)
            {
                A.Coefs[i, k] -= T[i, j] * A[j, k];
                B.Coefs[j, k] += T[i, j] * B[i, k];
            }
        }
    }
}

BigCplx[] PSLQM1(KMatrix<BigCplx> x, BigCplx gamma, int checkPos)
{
    // Initialize:
    // 1. For j := 1 to n: for i := 1 to n: if i = j then set Aij := 1 and Bij := 1
    // else set Aij := 0 and Bij := 0; endfor; endfor.
    var n = x.N;
    var z = BigCplx.BcZero(x.KOne.O);
    var A = new KMatrix<BigCplx>(z, n, n).One;
    var B = new KMatrix<BigCplx>(z, n, n).One;
    var T = new KMatrix<BigCplx>(z, n, n - 1);

    // 2. For k := 1 to n: set sk :=Sqrt(Sum[j=k to n] xj x*j) ; endfor; set t = 1/s1 ; for k := 1 to
    // n: set yk := txk ; sk := tsk ; endfor.
    var s = n.Range().Select(j => (n - j).Range(j).Aggregate(x.KZero, (acc, k) => acc + x[0, k].Absolute2))
        .Select(e => BigCplx.Sqrt(e)).ToArray();
    var t = s[0].Inv();
    var y = x.Select(xi => xi * t).ToKMatrix(n);
    s = s.Select(si => si * t).ToArray();

    // 3. Initial H: For j := 1 to n − 1: for i := 1 to j − 1: set Hij := 0; endfor;
    // set Hjj := sj+1 /sj ; for i := j + 1 to n: set Hij := −y*i yj /(sj sj+1 ); endfor;
    // endfor.
    var H = new KMatrix<BigCplx>(z, n, n - 1);
    for (int j = 0; j < n - 1; j++)
    {
        for (int i = 0; i < j - 1; i++)
            H.Coefs[i, j] = z;

        H.Coefs[j, j] = s[j + 1] / s[j];
        for (int i = j + 1; i < n; i++)
            H.Coefs[i, j] = (-y[i, 0].Conj * y[j, 0]) / (s[j] * s[j + 1]);
    }

    var step = 0;
    var O1 = z.O;
    var gamma_pow = (n - 1).Range().Select(i => (i, yi: gamma.Pow(i + 1).ToBigCplx(O1))).ToArray();

    // Iteration: Repeat the following steps until precision has been exhausted or a
    // relation has been detected.
    while (true)
    {
        ++step;
        IterOneLevelMultipair(H, A, B, T, y, gamma_pow);

        // 8. Norm bound: Compute M := 1/ maxj |Hjj |. Then there can exist no
        // relation vector whose Euclidean norm is less than M .
        // 9. Termination test: If the largest entry of A exceeds the level of numeric
        // precision used, then precision is exhausted. If the smallest entry of the y
        // vector is less than the detection epsilon, a relation has been detected and
        // is given in the corresponding row of B.
        if (A.Any(e => e.MaxV >= O1 + 2))
            throw new("Precision is exhausted");

        if (y.Select((yk, k) => (yk, k)).Any(e => e.yk.ToBigCplx(O1 - n / 2).IsZero() && !B[e.k, checkPos].IsZero()))
        {
            if (Logger.Level != LogLevel.Off)
                Console.WriteLine($"Possible Solution step:{step}");

            var (_, idx) = y.Select((yk, k) => (yk, k))
                .First(e => e.yk.ToBigCplx(O1 - n / 2).IsZero() && !B[e.k, checkPos].IsZero());
            
            // Console.WriteLine(B);
            // Console.WriteLine(y);
            // Console.WriteLine($"RowIdx:{idx}");
            
            Console.WriteLine("End PSLQ algorithm");
            Console.WriteLine($"Step {step}");
            return B.GetRow(idx).Select(c => c.RoundEven).ToArray();
        }
    }
    
    throw new("#1 Problem");
}

void testPSQL1()
{
    var n = 8;
    var O1 = 25;
    var O2 = O1 + n;
    var pi = BigCplx.Pi(O2);
    var beta = BigCplx.FromBigReal(
        BigReal.FromBigIntegerAndExponent(BigInteger.Parse("-1669947371922907049619"), 1, O1));

    var ai = n.Range().Select(k => pi.Pow(k).ToBigCplx(O1)).ToKMatrix();
    ai.Coefs[0, n - 1] = beta;

    // var y = 2 / BigCplx.Sqrt(BigCplx.FromBigInteger(3, O2));
    var y = 3 * pi.One / 2;
    // var y = BigCplx.FromBigInteger(2, O2);
    Console.WriteLine(ai);
    Console.WriteLine(y);
    GlobalStopWatch.AddLap();
    var coefs = PSLQM1(ai, y, n - 1).Select(e => e.RealPart.ToRational).ToArray();
    GlobalStopWatch.Show($"PSLQ");
    var P = -FG.KPoly('π', coefs.SkipLast(1).ToArray()) / coefs.Last();
    Console.WriteLine($"beta = {P}");
    Console.WriteLine($"diff = {P.Substitute(pi).ToBigCplx(O1) - beta}");
    Console.WriteLine();
}

void PSLQminPoly(int r, int s, int O)
{
    var n = r * s + 1; // Expected polynomial degree plus one
    var o = BigCplx.BcOne(O);
    var alpha = BigCplx.NthRoot(3 * o, r) - BigCplx.NthRoot(2 * o, s); // a = 3^(1/r) - 2^(1/s)
    var ai = n.Range().Select(k => alpha.Pow(k)).ToKMatrix();

    var O2 = O + n;
    // var y = 2 / BigCplx.Sqrt(BigCplx.FromBigInteger(3, O2));
    var y = 3 * BigCplx.BcOne(O2) / 2;
    // var y = BigCplx.FromBigInteger(2, O2);

    GlobalStopWatch.AddLap();
    var coefs = PSLQM1(ai, y, n - 1).Select(e => e.RealPart.ToRational).ToArray();
    GlobalStopWatch.Show($"PSLQ min poly a = 3^(1/{r}) - 2^(1/{s})");
    var P = FG.KPoly('X', coefs).Monic;
    Console.WriteLine($"P = {P}");
    Console.WriteLine($"P(a) = {P.Substitute(alpha).ToBigCplx(O / 2)}");
    Console.WriteLine();
    // GlobalStopWatch.AddLap();
    // AlgebraicIntegerRelationLLL.AlphaBetaPolynomial(alpha, alpha.Pow(n - 1), n, O);
    // GlobalStopWatch.Show("LLL");
    // Console.ReadLine();
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
    var O2 = O + d;
    // var y = 2 / BigCplx.Sqrt(3 * BigCplx.BcOne(O2));
    // var y = 3 * BigCplx.BcOne(O2) / 2;
    var y = BigCplx.Sqrt(2 * BigCplx.BcOne(O2));
    var rel = PSLQM1(mat, y, d - 1);
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

        if (Pbeta.IsZero())
        {
            Console.WriteLine($"[{col.Glue(",")}]");
            Console.WriteLine($"sum = {sum.RoundEven}");
            Console.WriteLine($"beta={betaAlg}");
            Console.WriteLine($"P(beta)={Pbeta}");
            Console.WriteLine();
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
}

void ExamplesGaloisPolynomial()
{
    var x = FG.QPoly();
    ExamplePol(x.Pow(4) - 2 * x.Pow(2) + 9, 20); // C2 x C2
    ExamplePol(x.Pow(4) + x.Pow(3) + x.Pow(2) + x + 1, 20); // C4

    ExamplePol(x.Pow(6) + 2 * x.Pow(5) - 2 * x.Pow(3) + 7 * x.Pow(2) + 8 * x + 13, 20); // C6
    ExamplePol(x.Pow(6) + 108, 40); // S3

    ExamplePol(x.Pow(4) - 8 * x.Pow(3) + 20 * x.Pow(2) - 16 * x + 2, 30); // C4
    ExamplePol(x.Pow(8) - 8 * x.Pow(6) + 20 * x.Pow(4) - 16 * x.Pow(2) + 2, 30); // C8
    ExamplePol(x.Pow(8) - x.Pow(4) + 1, 30); // C2 x C2 x C2
    ExamplePol(x.Pow(8) + 1, 30); // C2 x C4
    ExamplePol(x.Pow(8) + 4 * x.Pow(6) + 2 * x.Pow(4) + 28 * x.Pow(2) + 1, 30); // D8
    ExamplePol(x.Pow(8) - 12 * x.Pow(6) + 36 * x.Pow(4) - 36 * x.Pow(2) + 9, 30); // Q8;

    ExamplePol(x.Pow(10) - 2 * x.Pow(9) - 2 * x.Pow(8) + 6 * x.Pow(7) + 14 * x.Pow(6) - 26 * x.Pow(5) + 15 * x.Pow(4) -
        4 * x.Pow(3) + 58 * x.Pow(2) - 40 * x + 89, 60); // C10
    ExamplePol(x.Pow(10) + 10 * x.Pow(8) + 125 * x.Pow(6) + 500 * x.Pow(4) + 2500 * x.Pow(2) + 4000, 60); // D10

    ExamplePol(x.Pow(12) + 6 * x.Pow(8) + 26 * x.Pow(6) - 63 * x.Pow(4) + 162 * x.Pow(2) + 81, 100); // A4

    // ExamplePol(x.Pow(18) + 171 * x.Pow(12) + 5130 * x.Pow(6) + 27, 120); // C3 x: C6
    //
    // ExamplePol(x.Pow(20) + 2500 * x.Pow(10) + 50000, 120); // C5 x: C4
}

{
    GlobalStopWatch.Restart();
    Ring.DisplayPolynomial = MonomDisplay.StarCaret;

    StartingTest();
    ExamplesGaloisPolynomial();
    
    Console.Beep();
}
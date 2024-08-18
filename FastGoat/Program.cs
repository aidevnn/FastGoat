using System.CodeDom;
using System.Collections;
using System.Diagnostics;
using System.Numerics;
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

/*
The two-level multipair PSLQ algorithm
David H. Bailey
May 2, 2024
 */
Rational[] PSLQonelvl(KMatrix<BigReal> x, BigReal gamma)
{
    // Initialize:
    // 1. For j := 1 to n: for i := 1 to n: if i = j then set Aij := 1 and Bij := 1
    // else set Aij := 0 and Bij := 0; endfor; endfor.
    var n = x.N;
    var z = BigReal.BrZero(x.KOne.O);
    var A = new KMatrix<BigReal>(z, n, n).One;
    var B = new KMatrix<BigReal>(z, n, n).One;
    var T = new KMatrix<BigReal>(z, n, n);

    // 2. For k := 1 to n: set sk :=Sqrt(Sum[j=k to n] xj^2) ; endfor; set t = 1/s1 ; for k := 1 to
    // n: set yk := txk ; sk := tsk ; endfor.
    var s = n.Range().Select(j => (n - j).Range(j).Aggregate(x.KZero, (acc, k) => acc + x[0, k].Pow(2)))
        .Select(e => BigReal.Sqrt(e)).ToArray();
    var t = s[0].Inv();
    var y = x.Select(xi => xi * t).ToArray();
    s = s.Select(si => si * t).ToArray();

    // 3. Initial H: For j := 1 to n − 1: for i := 1 to j − 1: set Hij := 0; endfor;
    // set Hjj := sj+1 /sj ; for i := j + 1 to n: set Hij := −yi yj /(sj sj+1 ); endfor;
    // endfor.
    var H = new KMatrix<BigReal>(z, n, n - 1);
    for (int j = 0; j < n - 1; j++)
    {
        for (int i = 0; i < j - 1; i++)
            H.Coefs[i, j] = z;

        H.Coefs[j, j] = s[j + 1] / s[j];
        for (int i = j + 1; i < n; i++)
            H.Coefs[i, j] = (-y[i] * y[j]) / (s[j] * s[j + 1]);
    }

    var step = 0;
    var O1 = z.O;
    var O2 = gamma.O;
    var gamma_pow = (n - 1).Range().Select(i => (i, yi: gamma.Pow(i + 1))).ToArray();

    // Iteration: Repeat the following steps until precision has been exhausted or a
    // relation has been detected.
    while (true)
    {
        ++step;
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
        var m = gamma_pow.Select(e => (e.i, (e.yi * H[e.i, e.i].Absolute.ToBigReal(O2)).ToBigReal(O1)))
            .OrderByDescending(e => e.Item2).ToArray();
        var list = new List<(int, int)>();
        var selected = new HashSet<int>();
        for (int i = 0; i < n - 1; i++)
        {
            var (b, c) = (m[i].i, m[i].i + 1);
            if (selected.Contains(b) || selected.Contains(c))
                continue;

            list.Add((b, c));
            selected.UnionWith([b, c]);
            if (list.Count > 0.4 * n)
                break;
        }

        foreach (var (mi, _) in list)
        {
            (y[mi], y[mi + 1]) = (y[mi + 1], y[mi]);
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
        // i := mj to n: set t3 := Hi_mj ; t4 := Hi_mj1 ; Hi_mj := t1t3 + t2t4 ; and
        // Hi_mj1 := −t2t3 + t1t4 ; endfor; endif; endfor.
        for (int j = 0; j < list.Count; j++)
        {
            var mj = list[j].Item1;
            if (mj < n - 2)
            {
                var t0 = BigReal.Sqrt(H[mj, mj].Pow(2) + H[mj, mj + 1].Pow(2));
                var t1 = H[mj, mj] / t0;
                var t2 = H[mj, mj + 1] / t0;
                for (int i = mj; i < n; i++)
                {
                    var t3 = H[i, mj];
                    var t4 = H[i, mj + 1];
                    H.Coefs[i, mj] = t1 * t3 + t2 * t4;
                    H.Coefs[i, mj + 1] = -t2 * t3 + t1 * t4;
                }
            }
        }

        // 5. Reduce H: For i := 2 to n: for j := 1 to n − i + 1: set l := i + j − 1; for
        // k := j +1 to l −1: set Hlj := Hlj −Tlk Hkj ; endfor; set Tlj := nint(Hlj /Hjj )
        // and Hlj := Hlj − Tlj Hjj ; endfor; endfor. [Note that the n × (n − 1) integer
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

        for (int j = 1; j <= n - 1; j++)
        {
            for (int i = j + 1; i <= n; i++)
                y[j - 1] += T[i - 1, j - 1] * y[i - 1];
        }

        for (int k = 1; k <= n; k++)
        {
            for (int j = 1; j <= n - 1; j++)
            {
                for (int i = j + 1; i <= n; i++)
                {
                    A.Coefs[i - 1, k - 1] -= T[i - 1, j - 1] * A[j - 1, k - 1];
                    B.Coefs[j - 1, k - 1] += T[i - 1, j - 1] * B[i - 1, k - 1];
                }
            }
        }

        // 8. Norm bound: Compute M := 1/ maxj |Hjj |. Then there can exist no
        // relation vector whose Euclidean norm is less than M .
        // 9. Termination test: If the largest entry of A exceeds the level of numeric
        // precision used, then precision is exhausted. If the smallest entry of the y
        // vector is less than the detection epsilon, a relation has been detected and
        // is given in the corresponding row of B.

        var M = (n - 1).Range().Max(j => H[j, j]).Inv();
        if (A.Any(e => e.V >= O1 + 2))
            throw new("Precision is exhausted");

        if (y.Any(e => e.ToBigReal(O1 - 4).IsZero()))
        {
            var ym = y.Select((yi, i) => (e: yi, i)).OrderBy(c => c.e.Absolute).First();
            if (Logger.Level != LogLevel.Off)
            {
                Console.WriteLine($"Possible Solution step:{step}");
                Console.WriteLine(B.GetRow(ym.i).Select(c => c.RoundEven.ToRational).ToKMatrix());
                Console.WriteLine();
            }

            return B.GetRow(ym.i).Select(c => c.RoundEven.ToRational).ToArray();
        }
    }

    throw new();
}

void testPSQL1()
{
    var n = 8;
    var O1 = 20;
    var O2 = O1 + n;
    var pi = BigReal.Pi(O2);
    var beta = BigReal.FromBigIntegerAndExponent(BigInteger.Parse("-1669947371922907049619"), 1, O2);

    var ai = n.Range().Select(k => pi.Pow(k)).ToKMatrix();
    ai.Coefs[0, n - 1] = beta;

    var gamma = 2 / BigReal.Sqrt(BigReal.FromBigInteger(3, O2));
    // var gamma = BigReal.FromBigInteger(3, O2) / 2;

    GlobalStopWatch.AddLap();
    var coefs = PSLQonelvl(ai, gamma);
    Console.WriteLine(coefs.ToKMatrix());
    GlobalStopWatch.Show($"multipair PSLQ");
    var P = FG.KPoly('π', coefs.SkipLast(1).ToArray());
    Console.WriteLine($"beta = {-P / coefs.Last()}");
}

KPoly<Rational> LLLminPoly(int r, int s, int O)
{
    var n = r * s + 1; // Expected polynomial degree plus one 
    var alpha = BigReal.NthRoot(3, r, O) - BigReal.NthRoot(2, s, O); // a = 3^(1/r) - 2^(1/s)

    if (Logger.Level != LogLevel.Off)
        GlobalStopWatch.AddLap();
    
    var coefs = AlgebraicIntegerRelationLLL.AlphaBetaPolynomial(alpha, alpha.Pow(n - 1), n, O);
    var x = FG.QPoly('X');
    var P = x.Pow(n - 1) - coefs.Select((c, k) => c * x.Pow(k)).Aggregate((a, b) => a + b);
    if (Logger.Level != LogLevel.Off)
    {
        GlobalStopWatch.Show($"LLL min poly a = 3^(1/{r}) - 2^(1/{s})");
        Console.WriteLine($"P = {P} P(a) = {P.Substitute(alpha).ToBigReal(3 * O / 4)}");
        Console.WriteLine();
    }
    return P;
}

KPoly<Rational> PSLQonelvlminPoly(int r, int s, int O)
{
    var n = r * s + 1; // Expected polynomial degree plus one
    var alpha = BigReal.NthRoot(3, r, O) - BigReal.NthRoot(2, s, O); // a = 3^(1/r) - 2^(1/s)
    var ai = n.Range().Select(k => alpha.Pow(k)).ToKMatrix();

    var gamma = 2 / BigReal.Sqrt(BigReal.FromBigInteger(3, O));
    // var y = BigReal.FromBigInteger(3, O2) / 2;
    // var y = BigReal.FromBigInteger(2, O2);

    GlobalStopWatch.AddLap();
    var coefs = PSLQonelvl(ai, gamma);
    var P = FG.KPoly('X', coefs).Monic;
    GlobalStopWatch.Show($"PSLQ multipair min poly a = 3^(1/{r}) - 2^(1/{s})");
    if (Logger.Level != LogLevel.Off)
    {
        Console.WriteLine($"P = {P} and P(a) = {P.Substitute(alpha).ToBigReal(3 * O / 4)}");
        Console.WriteLine();
    }

    return P;
}

{
    GlobalStopWatch.Restart();
    Ring.DisplayPolynomial = MonomDisplay.StarCaret;

    testPSQL1();
    testPSQL1();
    
    // MinPoly of a = 3^(1/r) - 2^(1/s)
    List<(int r, int s, int O)> rsO = new()
    {
        (2, 2, 20),
        (3, 3, 40),
        (2, 5, 50),
        (3, 4, 70),
        (2, 7, 90),
        (3, 5, 90),
        (4, 4, 90),
        (4, 5, 120),
        (5, 5, 180),
    };

    Logger.Level = LogLevel.Level1;
    foreach (var (r, s, O) in rsO)
    {
        Console.WriteLine($"Min poly of a = 3^(1/{r}) - 2^(1/{s}) with LLL and multipair PSLQ algorithms on {O} digits");
        Console.WriteLine();

        var P1 = LLLminPoly(r, s, O);
        var P2 = PSLQonelvlminPoly(r, s, O);
        
        if (!P1.Equals(P2))
            throw new();
    }
    
    // Min poly of a = 3^(1/5) - 2^(1/5) with LLL and multipair PSLQ algorithms on 180 digits
    // 
    // # LLL  min poly a = 3^(1/5) - 2^(1/5) Time:1m19s
    // P = X^25 - 5*X^20 + 3760*X^15 + 11240*X^10 + 116255*X^5 - 1 P(a) = 0
    // 
    // # PSLQ multipair min poly a = 3^(1/5) - 2^(1/5) Time:1m46s
    // P = X^25 - 5*X^20 + 3760*X^15 + 11240*X^10 + 116255*X^5 - 1 and P(a) = 0
    // 
}
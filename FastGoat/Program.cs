using System.Diagnostics.Tracing;
using System.Globalization;
using System.Numerics;
using System.Security.Principal;
using FastGoat.Commons;
using FastGoat.Examples;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Characters;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Matrix;
using FastGoat.UserGroup.Perms;
using FastGoat.UserGroup.Polynoms;
using FastGoat.UserGroup.Words;
using FastGoat.UserGroup.Words.ToddCoxeter;
using static FastGoat.Commons.IntExt;
using static FastGoat.Commons.EnumerableExt;
using FastGoat.UserGroup.Padic;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");


KMatrix<BigReal> Hermite(KMatrix<BigReal> xi)
{
    var n = xi.N;
    var si = n.Range().Select(j => (n - j).Range(j).Aggregate(xi.KZero, (acc, k) => acc + xi[0, k].Pow(2)))
        .Select(e => BigReal.Sqrt(e)).ToArray();

    var H = new KMatrix<BigReal>(xi.KZero, n, n - 1);
    for (int i = 0; i < H.M; i++)
    {
        for (int j = 0; j < H.N; j++)
        {
            if (i == j)
                H.Coefs[i, j] = si[i + 1] / si[i];
            else if (j < i)
                H.Coefs[i, j] = -xi[0, i] * xi[0, j] / (si[j] * si[j + 1]);
        }
    }

    return H;
}

void ReducedHermite(KMatrix<BigReal> H0, KMatrix<BigReal> xi0, KMatrix<BigReal> A0, KMatrix<BigReal> B0)
{
    var n = H0.M;
    for (int i = 1; i < n; i++)
    {
        for (int j = i - 1; j >= 0; j--)
        {
            var q = (H0[i, j] / H0[j, j]).Round0;
            xi0.Coefs[0, j] += q * xi0[0, i];
            for (int k = 0; k < n; k++)
            {
                A0.Coefs[i, k] -= q * A0[j, k];
                B0.Coefs[k, j] += q * B0[k, i];

                if (k <= j)
                    H0.Coefs[i, k] -= q * H0[j, k];
            }
        }
    }
}

/*
ANALYSIS OF PSLQ, AN INTEGER
RELATION FINDING ALGORITHM
Helaman R. P. Ferguson
David H. Bailey
Steve Arno
03 July 1997
 */
Rational[] PSLQ(KMatrix<BigReal> ai, BigReal y)
{
    // Initialize
    var n = ai.N;
    var O1 = ai.KOne.O;
    var O2 = y.O;

    var A = new KMatrix<BigReal>(ai.KZero, n, n).One;
    var B = new KMatrix<BigReal>(ai.KZero, n, n).One;

    var nai = BigReal.NormN(2, ai.ToArray());
    var xi = ai / nai;
    var H = Hermite(xi);
    ReducedHermite(H, xi, A, B);

    var step = 0;
    var seq = (n - 1).Range();
    var ypow = seq.Select(i => (i, yi: y.Pow(i + 1))).ToArray();
    
    while (true)
    {
        ++step;
        
        // Step 1: Exchange
        var H0 = H;
        
        var (r, _) = ypow.Select(e => (e.i, (e.yi * H0[e.i, e.i].Absolute.ToBigReal(O2)).ToBigReal(O1))).MaxBy(e => e.Item2);
        var (b, c) = r < n - 2 ? (H[r + 1, r], H[r + 1, r + 1]) : (H.KZero, H.KZero);
        
        (xi.Coefs[0, r], xi.Coefs[0, r + 1]) = (xi[0, r + 1], xi[0, r]);
        for (int k = 0; k < n; k++)
        {
            if (k < n - 1)
                (H.Coefs[r, k], H.Coefs[r + 1, k]) = (H[r + 1, k], H[r, k]);
            
            (A.Coefs[r, k], A.Coefs[r + 1, k]) = (A[r + 1, k], A[r, k]);
            (B.Coefs[k, r], B.Coefs[k, r + 1]) = (B[k, r + 1], B[k, r]);
        }

        // Step 2: Corner
        if (r < n - 2)
        {
            var d = BigReal.Sqrt(b * b + c * c);
            var (e00, e01) = (b / d, -c / d);
            var (e10, e11) = (c / d, b / d);
            for (int k = r; k < n; k++)
            {
                var (h0, h1) = (H[k, r], H[k, r + 1]);
                (H.Coefs[k, r], H.Coefs[k, r + 1]) = (h0 * e00 + h1 * e10, h0 * e01 + h1 * e11);
            }
        }
        
        // Step 3: Reduction
        ReducedHermite(H, xi, A, B);
        
        // Step 4: Termination
        if (xi.Any(e => e.ToBigReal(O1 - 4).IsZero()))
            break;

        var H1 = H;
        if (seq.Select(k => H1[k, k]).Any(hkk => hkk.V > O1 || hkk.V < -O1))
            throw new("#0 Relation not found");
    }

    Console.WriteLine("xi");
    Console.WriteLine(xi.T);
    var ri = xi.ToList().FindIndex(e => e.ToBigReal(O1 - 4).IsZero());
    Console.WriteLine("H");
    Console.WriteLine(H);
    Console.WriteLine("B");
    Console.WriteLine(B.Select(b => b.ToRational.RoundEven).ToKMatrix(n));
    Console.WriteLine($"END step {step}");
    Console.WriteLine();
    var sol = B.GetCol(ri).Select(b => -b.ToRational.RoundEven).ToKMatrix();
    sol *= sol[0, n - 1].Sign;
    Console.WriteLine("Possible Solution");
    Console.WriteLine(sol);
    Console.WriteLine();
    return sol.ToArray();
}

void testPSQL1()
{
    var n = 8; 
    var O1 = 25;
    var O2 = O1 + n;
    var pi = BigReal.Pi(O2);
    var beta = BigReal.FromBigIntegerAndExponent(BigInteger.Parse("-1669947371922907049619"), 1, O2);
    
    var ai = n.Range().Select(k => pi.Pow(k).ToBigReal(O1)).ToKMatrix();
    ai.Coefs[0, n - 1] = beta;
    
    // var y = 2 / BigReal.Sqrt(BigReal.FromBigInteger(3, O2));
    var y = BigReal.FromBigInteger(3, O2) / 2;
    // var y = BigReal.FromBigInteger(2, O2);
    
    GlobalStopWatch.AddLap();
    var coefs = PSLQ(ai, y);
    GlobalStopWatch.Show($"PSLQ");
    Console.WriteLine($"beta = {FG.KPoly('π', coefs).Monic}");
}

void PSLQminPoly(int r, int s, int O)
{
    var n = r * s + 1; // Expected polynomial degree plus one
    var alpha = BigReal.NthRoot(3, r, O) - BigReal.NthRoot(2, s, O); // a = 3^(1/r) - 2^(1/s)
    var ai = n.Range().Select(k => alpha.Pow(k)).ToKMatrix();
    
    var O2 = O + n;
    
    // var y = 2 / BigReal.Sqrt(BigReal.FromBigInteger(3, O2));
    var y = BigReal.FromBigInteger(3, O2) / 2;
    // var y = BigReal.FromBigInteger(2, O2);
    
    GlobalStopWatch.AddLap();
    var coefs = PSLQ(ai, y);
    GlobalStopWatch.Show($"PSLQ min poly a = 3^(1/{r}) - 2^(1/{s})");
    Console.WriteLine($"P = {FG.KPoly('X', coefs)}");
    // GlobalStopWatch.AddLap();
    // AlgebraicIntegerRelationLLL.AlphaBetaPolynomial(alpha, alpha.Pow(n - 1), n, O);
    // GlobalStopWatch.Show("LLL");
    // Console.ReadLine();
}

{
    GlobalStopWatch.Restart();
    Ring.DisplayPolynomial = MonomDisplay.StarCaret;

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

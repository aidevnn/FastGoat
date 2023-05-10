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

(KMatrix<BigReal> H, KMatrix<BigReal> D, KMatrix<BigReal> Di) ReducedHermite(KMatrix<BigReal> H)
{
    var n = H.M;
    var D = new KMatrix<BigReal>(H.KZero, n, n).One;
    var H0 = new KMatrix<BigReal>(H.Coefs);
    for (int i = 1; i < n; i++)
    {
        for (int j = i - 1; j >= 0; j--)
        {
            var q = (H0[i, j] / H0[j, j]).Round0;
            // D.Coefs[i, j] -= q;
            // H0.Coefs[i, j] -= q * H0[j, j];
            for (int k = 0; k < n; k++)
            {
                D.Coefs[i, k] -= q * D[j, k];

                if (k <= j)
                    H0.Coefs[i, k] -= q * H0[j, k];
            }
        }
    }

    return (D * H, D, D.Inv());
}

KMatrix<BigReal> Corner(KMatrix<BigReal> H, int r)
{
    var n = H.M;
    var Q = new KMatrix<BigReal>(H.KZero, n - 1, n - 1).Zero;
    if (r == n - 2)
        return Q.One;

    var (b, c) = (H[r + 1, r], H[r + 1, r + 1]);
    var d = BigReal.Sqrt(b * b + c * c);

    for (int i = 0; i < Q.M; i++)
    {
        if (i != r && i != r + 1)
            Q.Coefs[i, i] = b.One;
    }

    Q.Coefs[r, r] = b / d;
    Q.Coefs[r, r + 1] = -c / d;
    Q.Coefs[r + 1, r] = c / d;
    Q.Coefs[r + 1, r + 1] = b / d;

    return Q;
}

/*
ANALYSIS OF PSLQ, AN INTEGER
RELATION FINDING ALGORITHM
Helaman R. P. Ferguson
David H. Bailey
Steve Arno
03 July 1997
 */
void PSLQ(KMatrix<BigReal> ai)
{
    // Initialize
    var n = ai.N;
    var O1 = ai.KOne.O;
    var O2 = O1 + 2 * n;
    var y = 2 / BigReal.Sqrt(BigReal.FromBigInteger(3, O2));
    // var y = BigReal.FromBigInteger(2, O2);

    var A = new KMatrix<BigReal>(ai.KZero, n, n).One;
    var B = new KMatrix<BigReal>(ai.KZero, n, n).One;

    var nai = BigReal.NormN(2, ai.ToArray());
    var xi = ai / nai;
    var (H, D, Di) = ReducedHermite(Hermite(xi));
    xi = xi * Di;
    A = D * A;
    B = B * Di;

    var step = 0;
    while (true)
    {
        ++step;
        Console.WriteLine(new { step });

        // Step 1: Exchange
        var H0 = H;
        var (r, _) = (n - 1).Range().Select(j => (j, y.Pow(j + 1).ToBigReal(O1) * H0[j, j].Absolute)).MaxBy(e => e.Item2);
        var R = new KMatrix<BigReal>(H.KZero, n, n).One;
        R.Coefs[r, r] = R.Coefs[r + 1, r + 1] = H.KZero;
        R.Coefs[r, r + 1] = R.Coefs[r + 1, r] = H.KOne;

        xi = xi * R;
        H = R * H;
        A = R * A;
        B = B * R;

        // Step 2: Corner
        var Q = Corner(H0, r);
        H = H * Q;

        // Step 3: Reduction
        (H, D, Di) = ReducedHermite(H);
        xi = xi * Di;
        A = D * A;
        B = B * Di;

        var H1 = H;
        // Step 4: Termination
        if (xi.Any(e => e.IsZero()) || (n - 1).Range().Any(k => H1[k, k].IsZero()))
            break;

        // Checking Properties
        if ((n - 1).Range().Grid2D().Any(e => e.t1 < e.t2 && !H1[e.t1, e.t2].IsZero()))
            throw new("#1 lower diagonal");

        if ((n - 1).Range().Grid2D().Any(e => e.t1 > e.t2 && H1[e.t1, e.t2].Absolute.CompareTo((H1[e.t2, e.t2] / 2).Absolute) > 0))
            throw new("#2 |Hij| < 1/2 * Hjj");

        // Console.WriteLine("H");
        // Console.WriteLine(H);
        // Console.WriteLine("B");
        // Console.WriteLine(B);
        // Console.WriteLine();
    }

    Console.WriteLine("H");
    Console.WriteLine(H);
    Console.WriteLine("B");
    Console.WriteLine(B);
    Console.WriteLine($"END step {step}");
    Console.WriteLine();
}

void TestPSLQ(int r, int O)
{
    // BigReal.Display = BigReal.DigitsForm.SciForm;
    var n = r * r + 1; // Expected polynomial degree plus one
    var alpha = BigReal.NthRoot(3, r, O) - BigReal.NthRoot(2, r, O); // a = 3^(1/r) - 2^(1/r)
    // AlgebraicIntegerRelationLLL.AlphaBetaPolynomial(alpha, alpha.Pow(n - 1), n, O1);
    var ai = n.Range().Select(k => alpha.Pow(k)).ToKMatrix();
    PSLQ(ai);
}

{
    TestPSLQ(2, 20); // step 11
    TestPSLQ(3, 40); // step 188
    TestPSLQ(4, 80); // step 1391
}

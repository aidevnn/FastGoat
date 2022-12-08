using System.ComponentModel.DataAnnotations;
using System.Diagnostics;
using System.IO.IsolatedStorage;
using System.Linq.Expressions;
using FastGoat.Commons;
using FastGoat.Examples;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Matrix;
using FastGoat.UserGroup.Perms;
using FastGoat.UserGroup.Polynoms;
using FastGoat.UserGroup.Words;
using FastGoat.UserGroup.Words.ToddCoxeter;
using static FastGoat.Commons.IntExt;
using static FastGoat.Commons.EnumerableExt;
using System.Numerics;
using System.Threading.Channels;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");
Random rnd = new Random();

KMatrix<Rational> RandMatrixRational(int n)
{
    var m = new KMatrix<Rational>(Rational.KZero(), n, n);
    while (true)
    {
        var m0 = Ring.Matrix(n, Rational.KZero(), (n * n).Range().Select(i => rnd.Next(n + 1)).ToArray());
        m = new(m0);
        if (!m.Det.IsZero())
            return m;
    }
}

Rational RoundDefault(Rational e)
{
    var (num, denom) = e;
    var (q, r) = BigInteger.DivRem(num, denom);
    var rs = r.Sign;
    var r0 = r * rs * 2;
    if (r0 < denom)
        return new(q, 1);

    return new(q + rs, 1);
}

Rational SquareNorm2(KMatrix<Rational> v)
{
    if (v.M == 1)
        return (v * v.T)[0, 0];
    else if (v.N == 1)
        return (v.T * v)[0, 0];

    throw new ArgumentException();
}

void SwapRows<T>(int i, int j, T[,] A)
{
    var cols = A.GetLength(1);
    for (int k = 0; k < cols; k++)
        (A[i, k], A[j, k]) = (A[j, k], A[i, k]);
}

Rational Round(Rational e)
{
    var (num, denom) = e;
    var (q, r) = BigInteger.DivRem(num, denom);
    var rs = r.Sign;
    var r0 = r * rs * 2;
    if (r0 < denom || (r0 == denom && BigInteger.IsEvenInteger(q)))
        return new(q, 1);

    return new(q + rs, 1);
}

KMatrix<Rational> LLL(KMatrix<Rational> v)
{
    var n = v.N;
    var w = v.Cols;
    var (Ws, M) = Ring.GramSchmidt(v);
    var v0 = new KMatrix<Rational>(Ws.Coefs);
    var ws = Ws.Cols;
    // Console.WriteLine((M * KMatrix<Rational>.MergeSameRows(ws).T).Equals(KMatrix<Rational>.MergeSameRows(w).T));
    var N = M.Coefs;
    int i = 1;
    while (i < n)
    {
        for (int j = i - 1; j >= 0; j--)
        {
            var ruij = Round(N[i, j]);
            w[i] -= ruij * w[j];
            for (int k = 0; k <= j; k++)
            {
                N[i, k] -= ruij * N[j, k];
            }
        }

        if (i >= 1)
        {
            var wsip2 = SquareNorm2(ws[i - 1]);
            var wsi2 = SquareNorm2(ws[i]);
            if (wsip2.CompareTo(2 * wsi2) > 0)
            {
                var a = N[i, i - 1];
                var b = a * wsip2 / (wsi2 + a.Pow(2) * wsip2);
                (ws[i - 1], ws[i]) = (ws[i] + a * ws[i - 1], ws[i - 1] - b * (ws[i] + a * ws[i - 1]));
                (w[i - 1], w[i]) = (w[i], w[i - 1]);
                SwapRows(i - 1, i, N);
                for (int k = i - 1; k < n; k++)
                {
                    (N[k, i - 1], N[k, i]) = (b * N[k, i - 1] + (1 - a * b) * N[k, i], N[k, i - 1] - a * N[k, i]);
                }

                i--;
            }
            else
            {
                i++;
            }
        }
        else
        {
            i++;
        }
    }

    // var B = KMatrix<Rational>.MergeSameRows(w);
    // var Bs = KMatrix<Rational>.MergeSameRows(ws);
    // Console.WriteLine("M * Ws.T = W.T {0}", (M * Bs.T).Equals(B.T));
    // var seq = n.Range().Grid2D(n.Range()).Where(e => e.t1 <= e.t2).ToArray();
    // Console.WriteLine(
    //     $"Norms Prop Before : {seq.All(e => SquareNorm2(v.GetCol(e.t1)).CompareTo(2.Pow(e.t2 - 1) * SquareNorm2(v0.GetCol(e.t2))) < 1)}");
    // Console.WriteLine(
    //     $"Norms Prop After  : {seq.All(e => SquareNorm2(w[e.t1]).CompareTo(2.Pow(e.t2 - 1) * SquareNorm2(ws[e.t2])) < 1)}");
    return KMatrix<Rational>.MergeSameRows(w);
}

{
    var m = Ring.Matrix(2, Rational.KZero(), 6, 8, 4, 4);
    var A = new KMatrix<Rational>(m);
    Console.WriteLine("Matrix A and LLL(A)");
    Console.WriteLine(A);
    Console.WriteLine(LLL(A));
}

{
    var m = Ring.Matrix(3, Rational.KZero(), 1, -1, 3, 1, 0, 5, 1, 2, 6);
    var A = new KMatrix<Rational>(m);
    Console.WriteLine("Matrix A and LLL(A)");
    Console.WriteLine(A);
    Console.WriteLine(LLL(A));
}

for (int k = 0; k < 20; ++k)
{
    var n = rnd.Next(10) + 2;
    Ring.MatrixDisplayForm = Ring.MatrixDisplay.SquareBracket;
    var A = RandMatrixRational(n);
    var detA = Rational.Abs(A.Det);
    var B = LLL(A);
    var detB = Rational.Abs(B.Det);
    Console.WriteLine("Matrix A");
    Console.WriteLine(A);
    Console.WriteLine("LLL(A)");
    Console.WriteLine(B);
    Console.WriteLine();
    Console.WriteLine("|Det(A)| = {0}; |Det(LLL(A))| = {1} {2}", detA, detB, detA.Equals(detB));
    Console.WriteLine("LLL(A) = LLL(LLL(A)) {0}", LLL(B).Equals(B));
}
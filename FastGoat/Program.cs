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

(KMatrix<K> O, KMatrix<K> U) GramSchmidt<K>(KMatrix<K> A, bool details = false)
    where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>, IFloatElt<K>
{
    var (m, n) = A.Dim;
    var vs = new KMatrix<K>[n];
    var u = Ring.Diagonal(A.KOne, n);
    var allCols = A.Cols;
    for (int i = 0; i < n; i++)
    {
        vs[i] = allCols[i];
        for (int j = 0; j < i; j++)
        {
            var sum0 = A.KZero;
            var sum1 = A.KZero;
            for (int k = 0; k < m; k++)
            {
                sum0 += vs[i].Coefs[k, 0] * vs[j].Coefs[k, 0].Conj;
                sum1 += vs[j].Coefs[k, 0] * vs[j].Coefs[k, 0].Conj;
            }

            // var uij = u[i, j] = (ScalarProduct(vs[i], vs[j]) / ScalarProduct(vs[j], vs[j]));
            var uij = u[i, j] = sum0 / sum1;

            // vs[i] -= uij * vs[j];
            for (int k = 0; k < m; k++)
                vs[i].Coefs[k, 0] -= uij * vs[j].Coefs[k, 0];
        }
    }

    return (KMatrix<K>.MergeSameRows(vs), new(u));
}

KMatrix<K> LLL<K>(KMatrix<K> A) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>, IFloatElt<K>
{
    var n = A.N;
    var w = A.Cols;
    var (Ws, M) = GramSchmidt(A);
    var ws = Ws.Cols;
    var N = new KMatrix<K>(M.Coefs);
    int i = 1;
    while (i < n)
    {
        var wi = w[i];
        for (int j = i - 1; j >= 0; j--)
        {
            var ruij = N[i, j].RoundEven;
            // w[i] -= ruij * w[j];
            var wj = w[j];
            for (int k = 0; k < n; k++)
                wi.Coefs[k, 0] -= ruij * wj.Coefs[k, 0];

            for (int k = 0; k <= j; k++)
                N.Coefs[i, k] -= ruij * N[j, k];
        }

        if (i >= 1)
        {
            var wsip2 = A.KZero;
            var wsi2 = A.KZero;
            var wsip = ws[i - 1];
            var wsi = ws[i];
            var wip = w[i - 1];
            for (int k = 0; k < n; k++)
            {
                wsip2 += wsip.Coefs[k, 0].Absolute2;
                wsi2 += wsi.Coefs[k, 0].Absolute2;
            }

            if (K.Abs(wsip2).CompareTo(K.Abs(2 * wsi2)) > 0)
            {
                var a = N[i, i - 1];
                var b = a * wsip2 / (wsi2 + a.Absolute2 * wsip2);
                for (int k = 0; k < n; k++)
                {
                    var tmp = wsi.Coefs[k, 0] + a * wsip.Coefs[k, 0];
                    (wsip.Coefs[k, 0], wsi.Coefs[k, 0]) = (tmp, wsip.Coefs[k, 0] - b * tmp);
                    (wip.Coefs[k, 0], wi.Coefs[k, 0]) = (wi.Coefs[k, 0], wip.Coefs[k, 0]);
                }

                for (int k = 0; k < n; k++)
                    (N.Coefs[i - 1, k], N.Coefs[i, k]) = (N[i, k], N[i - 1, k]);

                var coef = 1 - a * b;
                for (int k = i - 1; k < n; k++)
                    (N.Coefs[k, i - 1], N.Coefs[k, i]) = (b * N[k, i - 1] + coef * N[k, i], N[k, i - 1] - a * N[k, i]);

                i--;
            }
            else
                i++;
        }
        else
            i++;
    }

    return KMatrix<K>.MergeSameRows(w);
}

BigCplx[] AlphaBetaPolynomial(BigCplx alpha, BigCplx beta, int d, int O)
{
    if (Logger.Level != LogLevel.Off)
        Console.WriteLine("Start LLL algorithm");
    
    var N = BigCplx.FromBigReal(BigReal.Pow10(O, O));
    var mat = new KMatrix<BigCplx>(N.Zero, d, d).Zero;
    var ai = alpha.One;
    var aipow = new List<BigCplx>();
    for (int i = 0; i < mat.M - 1; i++)
    {
        aipow.Add(ai);
        mat.Coefs[i, i] = N.One;
        mat.Coefs[i, mat.N - 1] = (ai * N).RoundEven;
        ai *= alpha;
    }
    
    mat.Coefs[mat.N - 1, mat.N - 1] = (beta * N).RoundEven;
        
    var lll = IntFactorisation.LLL(mat.T);

    Console.WriteLine(mat);
    Console.WriteLine();
    Console.WriteLine(lll);

    var col = lll.Cols
        .OrderBy(l => l.SkipLast(1).Zip(aipow).Aggregate(-beta, (acc, v) => acc + v.First * v.Second).Magnitude2)
        .First();
    
    Console.WriteLine("End LLL algorithm");
    Console.WriteLine("Possible Solution");
    Console.WriteLine(col.T);
    Console.WriteLine();

    return col.SkipLast(1).Select(c => -c.RoundEven).ToArray();
}

void TestGramSchimdt()
{
    RngSeed(2514);
    {
        Rational[] A0 = ["1", "2", "3", "4"];
        var A = A0.ToKMatrix(2);
        Console.WriteLine("A");
        Console.WriteLine(A);
        Console.WriteLine();
        Console.WriteLine(A.Cols.Select(e => $"Matrix([{e.Glue(",")}])").Glue(", ")); // export sympy
        Console.WriteLine();

        Console.WriteLine("GS");
        var (O, U) = GramSchmidt(A);
        Console.WriteLine("O");
        Console.WriteLine(O);
        Console.WriteLine("U");
        Console.WriteLine(U);
        Console.WriteLine();
    }

    {
        Rational[] A0 = ["3", "2", "1", "2"];
        var A = A0.ToKMatrix(2);
        Console.WriteLine("A");
        Console.WriteLine(A);
        Console.WriteLine();
        Console.WriteLine(A.Cols.Select(e => $"Matrix([{e.Glue(",")}])").Glue(", ")); // export sympy
        Console.WriteLine();

        Console.WriteLine("GS");
        var (O, U) = GramSchmidt(A);
        Console.WriteLine("O");
        Console.WriteLine(O);
        Console.WriteLine("U");
        Console.WriteLine(U);
        Console.WriteLine();
    }

    {
        var i = Cplx.I;
        var A = 4.Range().Select(_ => Rng.Next(-1, 2) + Rng.Next(-1, 2) * i).ToKMatrix(2);
        Console.WriteLine("A");
        Console.WriteLine(A);
        Console.WriteLine();
        Console.WriteLine(A.Cols.Select(e => $"Matrix([{e.Glue(",")}])").Glue(", ")); // export sympy
        Console.WriteLine();

        Console.WriteLine("GS");
        var (O, U) = GramSchmidt(A);
        Console.WriteLine("O");
        Console.WriteLine(O);
        Console.WriteLine("U");
        Console.WriteLine(U);
        Console.WriteLine();
    }


    {
        var i = Cplx.I;
        var A = 9.Range().Select(_ => Rng.Next(-1, 2) + Rng.Next(-1, 2) * i).ToKMatrix(3);
        Console.WriteLine("A");
        Console.WriteLine(A);
        Console.WriteLine();
        Console.WriteLine(A.Cols.Select(e => $"Matrix([{e.Glue(",")}])").Glue(", ")); // export sympy
        Console.WriteLine();

        Console.WriteLine("GS");
        var (O, U) = GramSchmidt(A);
        Console.WriteLine("O");
        Console.WriteLine(O);
        Console.WriteLine("U");
        Console.WriteLine(U);
        Console.WriteLine();
    }
}

void TestLLL()
{
    {
        var o = Dble.DbleOne();
        var A = new[] { 1, -1, 3, 1, 0, 5, 1, 2, 6 }.Select(e => e * o).ToKMatrix(3);
        Console.WriteLine("A");
        Console.WriteLine(A);
        Console.WriteLine("LLL");
        Console.WriteLine(IntFactorisation.LLL(A));
        Console.WriteLine("LLL2");
        Console.WriteLine(LLL(A));
        Console.WriteLine();
    }

    {
        var i = Cplx.I;
        var o = Cplx.COne;
        Cplx[] A0 =
        [
            -2 + 2 * i, 7 + 3 * i, 7 + 3 * i, -5 + 4 * i,
            3 + 3 * i, -2 + 4 * i, 6 + 2 * i, -1 + 4 * i,
            2 + 2 * i, -8 * o, -9 + i, -7 + 5 * i,
            8 + 2 * i, -9 * o, 6 + 3 * i, -4 + 4 * i
        ];

        var A = A0.ToKMatrix(4);
        Console.WriteLine("A");
        Console.WriteLine(A);
        Console.WriteLine("LLL2");
        Console.WriteLine(LLL(A));
        Console.WriteLine();
    }
}

void ExamplePol()
{
    Ring.DisplayPolynomial = MonomDisplay.StarCaret;
    var x = FG.QPoly();
    var P = x.Pow(6) + 108;
    Console.WriteLine(P);

    var O1 = 50; // rounding digits
    var O2 = 60; // maximum precision digits

    var nRoots = FG.NRoots(P.ToBcPoly(O2));
    nRoots.Println();
    var alpha = nRoots[0];
    Logger.Level = LogLevel.Level1;
    foreach (var beta in nRoots.Skip(1))
    {
        var rel = AlphaBetaPolynomial(alpha, beta, P.Degree + 1, O1);
        var sum = rel.Select((c, i) => c * alpha.Pow(i)).Aggregate((a0, a1) => a0 + a1);
        Console.WriteLine($"alpha={alpha}");
        var fact = (sum / beta).ToBigCplx(O1 - 10);
        var test = fact.IsReal() && (fact.RealPart - fact.RealPart.RoundEven).IsZero4d();
        Console.WriteLine($"beta={beta} sum={sum} fact={fact}");
        var r = fact.RealPart.RoundEven.ToRational;
        var P1 = rel.Select((c, i) => c.RealPart.ToRational * x.Pow(i)).Aggregate((a0, a1) => a0 + a1) / r;
        if (test)
            Console.WriteLine($"beta = P1(alpha) with P1={P1}");
        else
            Console.WriteLine("#### Problem");
        Console.WriteLine();
    }
}

{
    TestGramSchimdt();
    TestLLL();
    ExamplePol();
}
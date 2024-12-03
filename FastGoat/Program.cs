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

K InnerProduct<K>(IEnumerable<K> A, IEnumerable<K> B) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>, IFloatElt<K>
{
    return A.Zip(B).Select(e => e.First * e.Second.Conj).Aggregate((a0, a1) => a0 + a1);
}

KMatrix<K> LLL<K>(KMatrix<K> A) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>, IFloatElt<K>
{
    var n = A.N;
    var w = A.Cols;
    var (Ws, M) = GramSchmidt(A);
    var ws = Ws.Cols;
    var N = new KMatrix<K>(M.Coefs);
    int i = 1;
    var d = A.KOne / 2;
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

            var a = N[i, i - 1];
            if (K.Abs(wsip2) > 2 * K.Abs(wsi2))
            {
                var b = (a * wsip2 / (wsi2 + a.Absolute2 * wsip2)).Conj; // Solve everything
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

(BigCplx[] col, BigCplx sum)[] AlphaBetaPolynomial(BigCplx alpha, BigCplx beta, int d)
{
    if (Logger.Level != LogLevel.Off)
        Console.WriteLine("Start LLL algorithm");
    
    var z = BigCplx.BcZero(alpha.O);
    var N = BigCplx.FromBigReal(BigReal.Pow10(alpha.O - 10, alpha.O));
    var mat = new KMatrix<BigCplx>(z.Zero, d, d).Zero;
    var ai = alpha.One;
    var aipow = new List<BigCplx>();
    for (int i = 0; i < mat.M - 1; i++)
    {
        aipow.Add(ai);
        mat.Coefs[i, i] = z.One;
        mat.Coefs[i, mat.N - 1] = (ai * N).RoundEven;
        ai *= alpha;
    }

    mat.Coefs[mat.M - 1, mat.N - 1] = (beta * N).RoundEven;

    var lll = LLL(mat.T);

    Console.WriteLine(mat);
    Console.WriteLine();
    Console.WriteLine(lll);

    var cols = lll.Cols
        .Select(
            l => (l, sum: l.SkipLast(1).Zip(aipow).Aggregate(z.Zero, (acc, v) => acc + v.First * v.Second) / beta)
        )
        .OrderBy(e => e.l.SkipLast(1).Aggregate(z.Zero, (acc, f) => acc + f.Absolute2))
        .Where(e => BigCplx.Abs(e.sum) < double.Exp10(alpha.O / 10.0) && !e.sum.ToBigCplx(alpha.O / 2).IsZero() && e.sum.IsGaussian())
        .ToArray();
    Console.WriteLine("End LLL algorithm");
    return cols.Select(e => (col: e.l.SkipLast(1).Select(c => c.RoundEven).ToArray(), e.sum.RoundEven)).ToArray();
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

void ExamplePol(KPoly<Rational> P, int O)
{
    Ring.DisplayPolynomial = MonomDisplay.StarCaret;
    Console.WriteLine(P);
    var (X, a) = FG.EPolyXc(P, 'a');
    var l = IntFactorisation.AlgebraicRoots(X.Pow(2) + 1);
    var nRoots = FG.NRoots(P.ToBcPoly(O));
    nRoots.Println();
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
    
    while(true)
    {
        var sz = 0;
        Console.WriteLine($"remRoots:{remRoots.Count}");
        if (remRoots.Count == 0)
            break;
        
        var beta = remRoots.First();
        var rels = AlphaBetaPolynomial(alpha, beta, P.Degree + 1);
        Console.WriteLine($"Possibles Solutions:{rels.Length}");
        foreach (var cs in rels)
        {
            var (col, sum) = cs;
            if (sum.RealPart.IsZero())
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
                Console.WriteLine(sum);
                Console.WriteLine($"beta={betaAlg}");
                Console.WriteLine($"P(beta)={Pbeta}");
                Console.WriteLine();
                set.UnionWith(Group.KAut(set.Append(betaAlg).ToArray()).Select(e => e.E));
                remRoots.ExceptWith(set.Select(e => e.Poly.Substitute(alpha).ToBigCplx(O / 2))
                    .Select(e => dicRoots[e]));
            }
        }
    }
    
    var P2 = set.Select(e => X - e).Aggregate((x0, x1) => x0 * x1);
    if (!P2.Equals(P.Substitute(X)))
        throw new();
    
    DisplayGroup.HeadElementsNames(Group.KAut(set.ToArray()));
    Console.WriteLine();
}

{
    // TestGramSchimdt();
    // TestLLL();

    var x = FG.QPoly();
    ExamplePol(x.Pow(4) - 2 * x.Pow(2) + 9, 20); // C2 x C2
    ExamplePol(x.Pow(4) + x.Pow(3) + x.Pow(2) + x + 1, 20); // C4
    ExamplePol(x.Pow(6) + 108, 20); // S3
    ExamplePol(x.Pow(8) + 4 * x.Pow(6) + 2 * x.Pow(4) + 28 * x.Pow(2) + 1, 30); // D8
    
    ExamplePol(x.Pow(10) + 10 * x.Pow(8) + 125 * x.Pow(6) + 500 * x.Pow(4) + 2500 * x.Pow(2) + 4000, 60); // D10
    ExamplePol(x.Pow(12) + 6 * x.Pow(8) + 26 * x.Pow(6) - 63 * x.Pow(4) + 162 * x.Pow(2) + 81, 100); // A4
    ExamplePol(x.Pow(18) + 171 * x.Pow(12) + 5130 * x.Pow(6) + 27, 120); // C3 x: C6
    ExamplePol(x.Pow(20) + 2500 * x.Pow(10) + 50000, 120); // C5 x: C4
}
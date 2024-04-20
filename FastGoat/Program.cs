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
using FastGoat.UserGroup.GModuleN;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Matrix;
using FastGoat.UserGroup.Perms;
using FastGoat.UserGroup.Polynoms;
using FastGoat.UserGroup.Words;
using static FastGoat.Commons.IntExt;
using static FastGoat.Commons.EnumerableExt;
using FastGoat.UserGroup.Padic;
using FastGoat.UserGroup.Words.Tools;
using GroupRegX = System.Text.RegularExpressions.Group;
using OrdMats = System.Collections.Generic.Dictionary<int, FastGoat.UserGroup.Matrix.Mat[]>;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");

void SymbolicDeterminantN(int n)
{
    var n2 = n * n;
    var xs = Ring.Polynomial(Rational.KZero(), MonomOrder.Lex, (n2, "a"));
    var M = xs.ToKMatrix(n);
    var det = Ring.Determinant(M.Coefs, M.KZero);
    // var cof = Ring.CoMatrix(M.Coefs, M.KZero);

    var rg = n.Range();
    Console.WriteLine($"int Det{n}x{n}(int[] mat)");
    Console.WriteLine("{");
    for (int i = 0; i < n; i++)
    {
        var x0s = xs.Skip(i * n).Take(n).ToArray();
        var i0 = i;
        Console.WriteLine($"    var ({x0s.Glue(", ")}) = ({rg.Select(j => i0 * n + j).Glue(", ", "mat[{0}]")});");
    }

    Console.WriteLine($"    var det = {det};");
    Console.WriteLine("    return ModP(det);");
    Console.WriteLine("}");
    Console.WriteLine();
}

void SymbolicInverseN(int n)
{
    var n2 = n * n;
    var xs = Ring.Polynomial(Rational.KZero(), MonomOrder.Lex, (n2, "a"));
    var M = xs.ToKMatrix(n);
    // var det = Ring.Determinant(M.Coefs, M.KZero);
    var com = Ring.CoMatrix(M.T.Coefs, M.KZero);

    var rg = n.Range();
    Console.WriteLine($"int Inv{n}x{n}(int[] mat, int[] inv)");
    Console.WriteLine("{");
    Console.WriteLine($"    var det = Det{n}x{n}(mat);");
    Console.WriteLine("    var idet = UnInvertible[det];");
    for (int i = 0; i < n; i++)
    {
        var x0s = xs.Skip(i * n).Take(n).ToArray();
        var i0 = i;
        Console.WriteLine($"    var ({x0s.Glue(", ")}) = ({rg.Select(j => i0 * n + j).Glue(", ", "mat[{0}]")});");
    }

    Console.WriteLine();
    var ys = xs.Select(x => $"{x}".Replace('a', 'b')).ToArray();
    for (int k = 0; k < n2; ++k)
    {
        var (i, j) = (k / n, k % n);
        Console.WriteLine($"    var {ys[k]} = inv[{k}] = ModP(({com[i, j]}) * idet);");
    }

    var s = $"{ys.Last()}";
    for (int i = n2 - 2; i >= 0; i--)
        s = i == 0 ? $"{ys[i]} + P * {s}" : $"({ys[i]} + P * {s})";
    
    Console.WriteLine($"    return {s};");
    Console.WriteLine("}");
    Console.WriteLine();
}

void GenerateCode(int max = 5)
{
    Ring.DisplayPolynomial = MonomDisplay.StarCaret;
    for (int i = 2; i <= max; ++i)
        SymbolicDeterminantN(i);
    
    for (int i = 2; i <= max; ++i)
        SymbolicInverseN(i);
}

Ring.MatrixDisplayForm = Ring.MatrixDisplay.SquareBracketNoFmt;

{
    var sl = FG.SL2p(3);
    DisplayGroup.HeadGenerators(sl);
}

void Sym6InGL62()
{
    var gl = new GL(6, 2);
    var a = gl[
        0, 1, 0, 0, 0, 0,
        1, 0, 0, 0, 0, 0,
        0, 0, 1, 0, 0, 0,
        0, 0, 0, 1, 0, 0,
        0, 0, 0, 0, 1, 0,
        0, 0, 0, 0, 0, 1
    ];
    var b = gl[
        0, 1, 0, 0, 0, 0,
        0, 0, 1, 0, 0, 0,
        0, 0, 0, 1, 0, 0,
        0, 0, 0, 0, 1, 0,
        0, 0, 0, 0, 0, 1,
        1, 0, 0, 0, 0, 0
    ];

    var s6 = Group.Generate("G", gl, a, b);
    DisplayGroup.HeadOrdersNames(s6);
}

{
    // GenerateCode(6);
    GlobalStopWatch.Bench(5, "Sym6InGL62", () => Sym6InGL62());
    GlobalStopWatch.Bench(5, "Sym6", () => DisplayGroup.HeadOrdersNames(FG.Symmetric(6)));
}
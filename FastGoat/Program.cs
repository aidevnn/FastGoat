using System.Globalization;
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
using System.Reflection;
using System.Runtime.InteropServices;
using System.Runtime.Intrinsics.Arm;
using System.Security.Cryptography.X509Certificates;
using System.Xml;
using FastGoat.UserGroup.Padic;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");

(EPoly<K> W, int l) Primitive3<K>(EPoly<K> U, EPoly<K> V) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
{
    var n = U.F.Degree;
    var vecV = V.Poly.ToVMatrix(n);
    for (int l = 1; l < 50; l++)
    {
        var W = U + l * V;
        var M = KMatrix<K>.MergeSameRows(n.Range().Select(i => W.Pow(i).Poly.ToVMatrix(n)).ToArray());
        var vM = KMatrix<K>.MergeSameRows(vecV, M);
        var dimKerM = M.NullSpace().nullity;
        var dimKerVM = vM.NullSpace().nullity;
        if (dimKerM == dimKerVM)
        {
            return (W, l);
        }
    }

    throw new();
}

// {
//     Ring.MatrixDisplayForm = Ring.MatrixDisplay.SquareBracket;
//     var x = FG.QPoly();
//     var (z, i) = FG.NumberFieldQ((x.Pow(2) - 3, "√3"), (x.Pow(2) + 1, "i"));
//     var one = z.One;
//     var gl = FG.GLnK($"Q({z})", 2, z);
//     var j = (1 + z * i) / 2;
//     var j2 = (1 - z * i) / 2;
//     
//     var A = gl[j, 0, 0, j2];
//     var B = gl[0, 1, 1, 0];
//     var C = gl[one / 2, z / 2, -z / 2, one / 2];
//     var D = gl[1, 0, 0, -1];
//     var GAB = Group.Generate("D12-AB", gl, A, B);
//     var GCD = Group.Generate("D12-CD", gl, C, D);
//     DisplayGroup.HeadElements(GAB);
//     DisplayGroup.HeadElements(GCD);
//
//     var d12 = FG.DihedralWg(6);
//     var a = d12["a"];
//     var b = d12["b"];
//     DisplayGroup.HeadElements(d12);
//
//     var pMapAB = Group.PartialMap((a, A), (b, B));
//     var homAB = Group.Hom(d12, Group.HomomorphismMap(d12, GAB, pMapAB));
//     Console.WriteLine(homAB.HomMap.Glue("\n"));
//     DisplayGroup.AreIsomorphics(d12, GAB);
//     Console.WriteLine();
//
//     var pMapCD = Group.PartialMap((a, C), (b, D));
//     var homCD = Group.Hom(d12, Group.HomomorphismMap(d12, GCD, pMapCD));
//     Console.WriteLine(homCD.HomMap.Glue("\n"));
//
//     DisplayGroup.AreIsomorphics(GAB, GCD);
// }

// {
//     Ring.MatrixDisplayForm = Ring.MatrixDisplay.SquareBracket;
//     var x = FG.QPoly();
//     var (z, i, X) = FG.NumberFieldQ((x.Pow(2) - 3, "√3"), (x.Pow(2) + 1, "i"), "X");
//     var one = z.One;
//     var gl = FG.GLnK($"Q({z})", 2, z);
//     var j = (1 + z * i) / 2;
//     var j2 = (1 - z * i) / 2;
//     
//     var A = gl[j, 0, 0, j2];
//     var B = gl[0, 1, 1, 0];
//     var C = gl[one / 2, z / 2, -z / 2, one / 2];
//     var D = gl[1, 0, 0, -1];
//
//     var matX = Ring.Diagonal(X, 2).ToKMatrix();
//     Console.WriteLine(A - matX);
//     Console.WriteLine((A - matX).Det);
//     Console.WriteLine(B - matX);
//     Console.WriteLine((B - matX).Det);
//     Console.WriteLine(C - matX);
//     Console.WriteLine((C - matX).Det);
//     Console.WriteLine(D - matX);
//     Console.WriteLine((D - matX).Det);
//
// }

{
    DisplayGroup.HeadElements(FG.AbelianWg(2, 2));
    DisplayGroup.HeadElements(FG.AbelianWg(2, 3));
    DisplayGroup.HeadElements(FG.AbelianWg(2, 3, 5));
    DisplayGroup.HeadElements(FG.AbelianWg(2, 3, 2, 3));

    // var a4 = FG.Alternate(4);
    // Group.DisplayOrbx(a4, Group.ByConjugate(a4));
    //
    // var a5 = FG.Alternate(5);
    // Group.DisplayOrbx(a5, Group.ByConjugate(a5));
    //
    // var a6 = FG.Alternate(6);
    // Group.DisplayOrbx(a6, Group.ByConjugate(a6));

    // var (X, T) = Ring.EPolynomial("X", "T", Rational.KZero(), MonomOrder.Lex);
    // var matB = Ring.Matrix(2, X, -1, -1, 1, 0).ToKMatrix();
    // var matX = Ring.Diagonal(X, 2).ToKMatrix();
    // Console.WriteLine(matB);
    // Console.WriteLine(matB - matX);
    // Console.WriteLine((matB - matX).Det);
    // Console.WriteLine(Ring.Determinant((matB - matX).Coefs, X));
}
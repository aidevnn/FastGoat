using System.Collections;
using System.Diagnostics;
using FastGoat.Commons;
using FastGoat.Examples;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Perms;
using FastGoat.UserGroup.Polynoms;
using FastGoat.UserGroup.Words;
using static FastGoat.Commons.IntExt;
using static FastGoat.Commons.EnumerableExt;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");
{
    var coefs = Ring.Polynomial(ZnInt.KZero(), "abcdefghijklmnop".ToArray());
    var z0 = Ring.PolynomialZero(ZnInt.KZero());
    var mat = Ring.Matrix(4, coefs);
    Ring.DisplayMatrix(mat);
    Console.WriteLine("Det = {0}", Ring.Determinant(mat, z0));
    Console.WriteLine();
}

{
    var coefs = Ring.Polynomial(ZnInt.KZero(), "abcdefghijklmnopqrstuvwxy".ToArray());
    var z0 = Ring.PolynomialZero(ZnInt.KZero());
    var mat = Ring.Matrix(5, coefs);
    Ring.DisplayMatrix(mat);
    Console.WriteLine("Det = {0}", Ring.Determinant(mat, z0));
    Console.WriteLine();
}
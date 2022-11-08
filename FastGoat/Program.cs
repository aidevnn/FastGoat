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
    var a = Ring.QPolynomial('a');
    var p = a.Pow(2) + 1;
    Console.WriteLine(p);

    var (x, a1, k) = Ring.ExtPolynomial(p, 'x');
    var x_ai = Enumerable.Range(0, k + 1).Select(e => x - a1.Pow(e)).ToArray();
    Console.WriteLine(x_ai.Glue("; "));

    Console.WriteLine(x_ai.Aggregate((xi, xj) => xi * xj));
}
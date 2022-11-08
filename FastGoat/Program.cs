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
    var p = a.Pow(2) + 5; // a*a-2, a*a+5
    Console.WriteLine(p);

    var (x, a1) = Ring.ExtPolynomial(p, 'x');
    Console.WriteLine(x - a1);
    Console.WriteLine(x + a1);
    Console.WriteLine((x - a1) * (x + a1));
}
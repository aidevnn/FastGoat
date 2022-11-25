using System.Diagnostics;
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

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");

{
    var x = FG.QPoly();
    Console.WriteLine(Ring.Discriminant(x.Pow(2) + 3 * x - 5));
    Console.WriteLine(Ring.Discriminant(x.Pow(3) + 4 * x + 12));
    Console.WriteLine(Ring.Discriminant(x.Pow(4) + 3 * x.Pow(2) + 1));
    Console.WriteLine(Ring.Discriminant(x.Pow(12) + x.Pow(11) - x.Pow(9) - 2 * x.Pow(8) + x.Pow(5) + x.Pow(4)));
    Console.WriteLine(Ring.Discriminant(x.Pow(12) - 3 * x.Pow(7) + 4));
}
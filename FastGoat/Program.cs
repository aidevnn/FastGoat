using System.ComponentModel.DataAnnotations;
using System.Diagnostics;
using System.Numerics;
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

{
    var x = Ring.Polynomial(Rational.KZero(), "x")[0];
    Console.WriteLine(Ring.LcmPolynomial((x + 1) * (x + 5).Pow(3) * (x - 6), (x + 2) * (x + 5).Pow(2)));
    Console.WriteLine((x + 1) * (x + 5).Pow(3) * (x + 2) * (x - 6));
    Console.WriteLine();
    Console.WriteLine(Ring.GcdPolynomial((x + 1) * (x + 5).Pow(3) * (x - 6), (x + 2) * (x + 5).Pow(2)));
    Console.WriteLine((x + 5).Pow(2));
    Console.WriteLine();
}

{
    var x = FG.QPoly();
    Console.WriteLine(Ring.Lcm((x + 1) * (x + 5).Pow(3) * (x - 6), (x + 2) * (x + 5).Pow(2)).Monic);
    Console.WriteLine((x + 1) * (x + 5).Pow(3) * (x + 2) * (x - 6));
    Console.WriteLine();
    Console.WriteLine(Ring.Gcd((x + 1) * (x + 5).Pow(3) * (x - 6), (x + 2) * (x + 5).Pow(2)).Monic);
    Console.WriteLine((x + 5).Pow(2));
    Console.WriteLine();
}

{
    var (x, y) = Ring.Polynomial(Rational.KZero(), "x", "y").Deconstruct();
    Console.WriteLine(Ring.LcmPolynomial((x + 1) * (y + 5).Pow(3) * (x - 6), (x + 2) * (y + 5).Pow(2)));
    Console.WriteLine(Ring.LcmPolynomial((x + 2) * (y + 5).Pow(2), (x + 1) * (y + 5).Pow(3) * (x - 6)));
    Console.WriteLine((x + 1) * (y + 5).Pow(3) * (x + 2) * (x - 6));
    Console.WriteLine();
    Console.WriteLine(Ring.GcdPolynomial((x + 1) * (y + 5).Pow(3) * (x - 6), (x + 2) * (y + 5).Pow(2)));
    Console.WriteLine(Ring.GcdPolynomial((x + 2) * (y + 5).Pow(2), (x + 1) * (y + 5).Pow(3) * (x - 6)));
    Console.WriteLine((y + 5).Pow(2));
    Console.WriteLine();
}
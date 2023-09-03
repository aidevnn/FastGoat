using System.ComponentModel.DataAnnotations;
using System.Diagnostics;
using System.Numerics;
using System.Security.Cryptography;
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
Ring.DisplayPolynomial = MonomDisplay.StarCaret;

{
    var x = FG.QPoly();
    var P = x.Pow(5) - x.Pow(4) - 4 * x.Pow(3) + 3 * x.Pow(2) + 3 * x - 1;
    Console.WriteLine(P);
    var P1 = P.Substitute(5 * x / 12);
    var P2 = P1.PrimitiveZPoly();
    Console.WriteLine(P1);
    Console.WriteLine(P2);
    Console.WriteLine(IntFactorisation.ConstCoef(P1));
    Console.WriteLine(IntFactorisation.ConstCoef(P2, monic: true));
}

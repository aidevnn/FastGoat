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
//
// {
//     Ring.DisplayPolynomial = MonomDisplay.Caret;
//     var x = FG.QPoly('X');
//     var (X, _) = FG.EPolyXc(x.Pow(2) - 2, 'a');
//     var (minPoly, a0, b0) = IntFactorisation.PrimitiveElt(X.Pow(2) - 3);
//     var roots = IntFactorisation.AlgebraicRoots(minPoly);
//     Console.WriteLine("Q(√2, √3) = Q(α)");
//     var gal = GaloisTheory.GaloisGroup(roots);
//     DisplayGroup.AreIsomorphics(gal, FG.Abelian(2, 2));
// }

// {
//     Ring.DisplayPolynomial = MonomDisplay.Caret;
//     var x = FG.QPoly('X');
//     var (X, i) = FG.EPolyXc(x.Pow(2) + 1, 'i');
//     var (minPoly, _, _) = IntFactorisation.PrimitiveElt(X.Pow(4) - 2);
//     var roots = IntFactorisation.AlgebraicRoots(minPoly);
//     Console.WriteLine("With α^4-2 = 0, Q(α, i) = Q(β)");
//     var gal = GaloisTheory.GaloisGroup(roots, 'β');
//     DisplayGroup.AreIsomorphics(gal, FG.Dihedral(4));
// }
//
// {
//     Ring.DisplayPolynomial = MonomDisplay.Caret;
//     var x = FG.QPoly('X');
//     var (X, i) = FG.EPolyXc(x.Pow(2) + 1, 'i');
//     var (minPoly, _, _) = IntFactorisation.PrimitiveElt(X.Pow(4) - 2);
//     var roots = IntFactorisation.AlgebraicRoots(minPoly);
//     Console.WriteLine("With α^4-2 = 0, Q(α, i) = Q(β)");
//     var gal = GaloisTheory.GaloisGroup(roots, 'β');
//     DisplayGroup.AreIsomorphics(gal, FG.Dihedral(4));
// }

// {
//     Ring.DisplayPolynomial = MonomDisplay.Caret;
//     
//     var x = FG.QPoly('X');
//     // var P = x.Pow(5) + x.Pow(4) - 4 * x.Pow(3) - 3 * x.Pow(2) + 3 * x + 1;
//     // var P = x.Pow(8) + 1;
//     // var P = x.Pow(9) + 6 * x.Pow(8) - 6 * x.Pow(7) - 53 * x.Pow(6) + 45 * x.Pow(5) + 135 * x.Pow(4) - 197 * x.Pow(3)
//     //     + 66 * x.Pow(2) + 3 * x - 1;
//     // var P = x.Pow(7) + x.Pow(6) - 18 * x.Pow(5) - 35 * x.Pow(4) + 38 * x.Pow(3) + 104 * x.Pow(2) + 7 * x - 49;
//     // var P = x.Pow(8) - x.Pow(4) + 1;
//     // var P = x.Pow(8) - x.Pow(7) - 7 * x.Pow(6) + 6 * x.Pow(5) + 15 * x.Pow(4) - 10 * x.Pow(3) - 10 * x.Pow(2) + 4 * x + 1;
//     // var P = x.Pow(8) + 28*x.Pow(4) + 2500;
//     // var P = x.Pow(10) - 2 * x.Pow(9) - 20 * x.Pow(8) + 2 * x.Pow(7) + 69 * x.Pow(6) - x.Pow(5) - 69 * x.Pow(4)
//     //     + 2 * x.Pow(3) + 20 * x.Pow(2) - 2 * x - 1;
//     var P = x.Pow(8) + 24 * x.Pow(4) + 16;
//     // var P = x.Pow(8) + 4 * x.Pow(6) + 2 * x.Pow(4) + 28 * x.Pow(2) + 1;
//     // var P = x.Pow(8) - 12 * x.Pow(6) + 36 * x.Pow(4) - 36 * x.Pow(2) + 9;
//     // var P = x.Pow(9) + x.Pow(8) - 8 * x.Pow(7) - 7 * x.Pow(6) + 21 * x.Pow(5) + 15 * x.Pow(4)
//     //     - 20 * x.Pow(3) - 10 * x.Pow(2) + 5 * x + 1;
//
//     Console.WriteLine($"Polynomial P={P}");
//     Console.WriteLine();
//     var roots = IntFactorisation.AlgebraicRoots(P, details: true);
//     var gal = GaloisTheory.GaloisGroup(roots);
//     DisplayGroup.AreIsomorphics(gal, FG.Dihedral(4));
//     // Console.WriteLine();
// }


{
    Ring.DisplayPolynomial = MonomDisplay.Caret;
    var x = FG.QPoly('X');
    var P = x.Pow(5) + x.Pow(4) - 4 * x.Pow(3) - 3 * x.Pow(2) + 3 * x + 1;
    var roots = IntFactorisation.AlgebraicRoots(P, details: true);
    var gal = GaloisTheory.GaloisGroup(roots);
    DisplayGroup.AreIsomorphics(gal, FG.Abelian(5));
}

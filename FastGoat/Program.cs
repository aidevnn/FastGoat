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
    var x = FG.QPoly();
    // var (X, y) = FG.EPolyXc(x.Pow(5) + x.Pow(4) - 4 * x.Pow(3) - 3 * x.Pow(2) + 3 * x + 1, 'y');
    // var (X, y) = FG.EPolyXc(x.Pow(8) + 1, 'y');
    // var (X, y) = FG.EPolyXc(x.Pow(8) - x.Pow(7) - 7 * x.Pow(6) + 6 * x.Pow(5) + 15 * x.Pow(4)
    //     - 10 * x.Pow(3) - 10 * x.Pow(2) + 4 * x + 1, 'y');
    // var (X, y) = FG.EPolyXc(x.Pow(9) + 6 * x.Pow(8) - 6 * x.Pow(7) - 53 * x.Pow(6) + 45 * x.Pow(5) + 135 * x.Pow(4) - 197 * x.Pow(3)
    //     + 66 * x.Pow(2) + 3 * x - 1, 'y');
    // var (X, y) = FG.EPolyXc(x.Pow(7) + x.Pow(6) - 18 * x.Pow(5) - 35 * x.Pow(4) + 38 * x.Pow(3) + 104 * x.Pow(2) + 7 * x - 49, 'y');
    // var (X, y) = FG.EPolyXc(x.Pow(8) - x.Pow(4) + 1, 'y');
    // var (X, y) = FG.EPolyXc(x.Pow(8) - x.Pow(7) - 7 * x.Pow(6) + 6 * x.Pow(5) + 15 * x.Pow(4)
    //     - 10 * x.Pow(3) - 10 * x.Pow(2) + 4 * x + 1, 'y');
    // var (X, y) = FG.EPolyXc(x.Pow(8) + 28*x.Pow(4) + 2500, 'y');
    var (X, y) = FG.EPolyXc(x.Pow(10) - 2 * x.Pow(9) - 20 * x.Pow(8) + 2 * x.Pow(7) + 69 * x.Pow(6) - x.Pow(5) - 69 * x.Pow(4)
        + 2 * x.Pow(3) + 20 * x.Pow(2) - 2 * x - 1, 'y');
    // var (X, y) = FG.EPolyXc(x.Pow(8) + 24 * x.Pow(4) + 16, 'y');
    // var (X, y) = FG.EPolyXc(x.Pow(8) - 12 * x.Pow(6) + 36 * x.Pow(4) - 36 * x.Pow(2) + 9, 'y');
    // var (X, y) = FG.EPolyXc(x.Pow(9) + x.Pow(8) - 8 * x.Pow(7) - 7 * x.Pow(6) + 21 * x.Pow(5) + 15 * x.Pow(4) - 20 * x.Pow(3)
    //     - 10 * x.Pow(2) + 5 * x + 1, 'y');

    // IntFactorisation.AlgebraicFactors(X.Pow(8) - 2, true);
    var roots = IntFactorisation.AlgebraicRoots(y.F.Substitute(X), true);
    var gal = GaloisTheory.GaloisGroup(roots);
    // DisplayGroup.AreIsomorphics(gal, FG.Symmetric(3));
    // Console.WriteLine();
}
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

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");

{
    Console.WriteLine("P-Adic Rationals");
    var q1 = new Rational(313, 3000);
    var q2 = new Rational(-313, 3000);
    var q3 = new Rational(53, 3750);
    var q4 = new Rational(-61, 2250);
    foreach (var q in new[] { q1, q2, q3, q4 })
        Console.WriteLine("{0} = {1}", q, Padic.Convert(5, 7, q));
    
    foreach (var q in new[] { q1, q2, q3, q4 })
        Console.WriteLine("{0} = {1}", q, Padic.Convert(3, 5, q));
    
    // 313/3000 = [2244.343434(5⁷)]
    // -313/3000 = [3200.10101(5⁷)]
    // 53/3750 = [32404.040404(5⁷)]
    // -61/2250 = [3423.411033(5⁷)]
    // 313/3000 = [12.1111(3⁵)]
    // -313/3000 = [20.1111(3⁵)]
    // 53/3750 = [10.1202(3⁵)]
    // -61/2250 = [222(3⁵)]

    // gap> fam:=PurePadicNumberFamily(5,7);;
    // gap> PadicNumber(fam, 313/3000);PadicNumber(fam, -313/3000);PadicNumber(fam, 53/3750);PadicNumber(fam, -61/2250);
    // 2244.343(5)
    // 3200.101(5)
    // 32404.04(5)
    // 3423.411(5)
    //
    // gap> fam:=PurePadicNumberFamily(3,5);;
    // gap> PadicNumber(fam, 313/3000);PadicNumber(fam, -313/3000);PadicNumber(fam, 53/3750);PadicNumber(fam, -61/2250);
    // 12.111(3)
    // 20.111(3)
    // 10.12(3)
    // 222(3)

}
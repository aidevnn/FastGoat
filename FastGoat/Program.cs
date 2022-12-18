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
    Console.WriteLine("Qp P-Adic Numbers");
    Console.WriteLine(new QpAdic(5, 3));
    Console.WriteLine();

    var a1 = new QpAdic(5, 3, 101);
    var a2 = new QpAdic(5, 3, 75);
    var a3 = new QpAdic(5, 3, 15);
    Console.WriteLine(a1);
    Console.WriteLine(a1 * 5.Pow(2));
    Console.WriteLine(a1 / 5.Pow(2));
    Console.WriteLine(a1 * 5.Pow(5));
    Console.WriteLine(a1 / 5.Pow(5));
    Console.WriteLine(a2);
    Console.WriteLine(a2 * 5);
    Console.WriteLine(a2 / 5);
    Console.WriteLine(a3);
    Console.WriteLine(a3 * 5);
    Console.WriteLine(a3 / 5);
    Console.WriteLine();

    var b1 = new QpAdic(5, 3, 171);
    var b2 = new QpAdic(5, 3, -171);
    var b3 = new QpAdic(5, 6, 171);
    var b4 = new QpAdic(5, 6, -171);
    Console.WriteLine(b1);
    Console.WriteLine(b2);
    Console.WriteLine(b2.Equals(b1.Opp()));
    Console.WriteLine(b3);
    Console.WriteLine(b4);
    Console.WriteLine(b3.Equals(b4.Opp()));
    Console.WriteLine(b1 * 5.Pow(4));
    Console.WriteLine(b1 / 5.Pow(4));
    Console.WriteLine();

    var c1 = new QpAdic(3, 3, 101);
    var c2 = new QpAdic(3, 6, 101);
    var c3 = new QpAdic(3, 4, -101);
    var c4 = c3.Inv();
    Console.WriteLine(c1);
    Console.WriteLine(c1 * 3.Pow(2));
    Console.WriteLine(c1 / 3.Pow(2));
    Console.WriteLine(c2);
    Console.WriteLine(c2 * 3.Pow(2));
    Console.WriteLine(c2 / 3.Pow(2));
    Console.WriteLine(c3);
    Console.WriteLine(c3 * 3.Pow(2));
    Console.WriteLine(c3 / 3.Pow(2));
    Console.WriteLine((c3 * c4).Equals(c3.One));
    Console.WriteLine((c3 * c3 * c3).Equals(c3.Pow(3)));
    Console.WriteLine((c4 * c4).Equals(c3.Pow(-2)));
    Console.WriteLine((c4 * c4).Equals(c3.Pow(2).Inv()));

    var d1 = new QpAdic(3, 7, new Rational(101));
    var d2 = new QpAdic(3, 7, new Rational(909));
    var d3 = new QpAdic(3, 7, new Rational(101, 9));
    Console.WriteLine(d1);
    Console.WriteLine(d1 * 3.Pow(2));
    Console.WriteLine(d1 / 3.Pow(2));
    Console.WriteLine(d2);
    Console.WriteLine(d2 * 3.Pow(3));
    Console.WriteLine(d2 / 3.Pow(3));
    Console.WriteLine(d3);
    Console.WriteLine(d3 * 3);
    Console.WriteLine(d3 / 3);
    Console.WriteLine(d2.Equals(9 * d1));
    Console.WriteLine(d3.Equals(d1 / 9));
    Console.WriteLine();

    var e1 = new QpAdic(3, 7, new Rational(909) + new Rational(101, 9));
    var e2 = new QpAdic(3, 11, new Rational(909) + new Rational(101, 9));
    var e3 = new QpAdic(5, 7, new Rational(313, 3000));
    Console.WriteLine(e1);
    Console.WriteLine(e1 * 3.Pow(5));
    Console.WriteLine(e1 / 3.Pow(5));
    Console.WriteLine();
    Console.WriteLine(e2 * 3.Pow(4));
    Console.WriteLine(e2 / 3.Pow(4));
    Console.WriteLine(e3);
    Console.WriteLine(e3 * 5.Pow(2));
    Console.WriteLine(e3 / 5.Pow(2));
    Console.WriteLine(e1.Equals(d1 / 9 + d2));
    Console.WriteLine(e1.Equals(d2 + d3));
    Console.WriteLine();

    var f1 = new QpAdic(5, 7, new Rational(313, 75)) * 5.Pow(4);
    var f2 = new QpAdic(5, 7, 415) / 5.Pow(2);
    Console.WriteLine(f1);
    Console.WriteLine(f2);
    Console.WriteLine(f1 + f2);
    Console.WriteLine((f1 + f2) * 125);
    Console.WriteLine(f1 * 125);
    Console.WriteLine((f1 + f2) * 3125);
    Console.WriteLine(f1 * 3125);
    Console.WriteLine(f1 * 3125 + f2);
    Console.WriteLine((f1 * 125 + f2 * 125).Equals((f1 + f2) * 125));
    Console.WriteLine();
    
    Console.WriteLine(f1);
    Console.WriteLine(f1.Norm);
    Console.WriteLine(f1.Normalized);
    Console.WriteLine(f1.Normalized / f1.Norm);
    Console.WriteLine(f1 * f1.Norm);
    Console.WriteLine(f1.Normalized.Equals(f1 * f1.Norm));
    Console.WriteLine();
}


{
    var q1 = new Rational(313, 3000);
    var q2 = new Rational(-313, 3000);
    var q3 = new Rational(53, 3750);
    var q4 = new Rational(-61, 2250);
    foreach (var q in new[] { q1, q2, q3, q4 })
        Console.WriteLine("{0} = {1}", q, new QpAdic(5, 7, q));
    
    foreach (var q in new[] { q1, q2, q3, q4 })
        Console.WriteLine("{0} = {1}", q, new QpAdic(3, 5, q));
    
    // 313/3000 = [2244.343(5)~7]
    // -313/3000 = [3200.101(5)~7]
    // 53/3750 = [32404.04(5)~7]
    // -61/2250 = [3423.411(5)~7]
    // 313/3000 = [12.111(3)~5]
    // -313/3000 = [20.111(3)~5]
    // 53/3750 = [10.12(3)~5]
    // -61/2250 = [222(3)~5]
}
//
// {
//     Console.WriteLine("P-Adic Rationals");
//     var q1 = new Rational(313, 3000);
//     var q2 = new Rational(-313, 3000);
//     var q3 = new Rational(53, 3750);
//     var q4 = new Rational(-61, 2250);
//     foreach (var q in new[] { q1, q2, q3, q4 })
//         Console.WriteLine("{0} = {1}", q, Padic.Convert(5, 7, q));
//     
//     foreach (var q in new[] { q1, q2, q3, q4 })
//         Console.WriteLine("{0} = {1}", q, Padic.Convert(3, 5, q));
//     
//     // 313/3000 = [2244.343434(5⁷)]
//     // -313/3000 = [3200.10101(5⁷)]
//     // 53/3750 = [32404.040404(5⁷)]
//     // -61/2250 = [3423.411033(5⁷)]
//     // 313/3000 = [12.1111(3⁵)]
//     // -313/3000 = [20.1111(3⁵)]
//     // 53/3750 = [10.1202(3⁵)]
//     // -61/2250 = [222(3⁵)]
//     
//     // gap> fam:=PurePadicNumberFamily(5,7);;
//     // gap> PadicNumber(fam, 313/3000);PadicNumber(fam, -313/3000);PadicNumber(fam, 53/3750);PadicNumber(fam, -61/2250);
//     // 2244.343(5)
//     // 3200.101(5)
//     // 32404.04(5)
//     // 3423.411(5)
//     //
//     // gap> fam:=PurePadicNumberFamily(3,5);;
//     // gap> PadicNumber(fam, 313/3000);PadicNumber(fam, -313/3000);PadicNumber(fam, 53/3750);PadicNumber(fam, -61/2250);
//     // 12.111(3)
//     // 20.111(3)
//     // 10.12(3)
//     // 222(3)
//
// }
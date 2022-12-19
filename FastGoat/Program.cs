using System.Globalization;
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
using System.Xml;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");

{
     Console.WriteLine("Qp Zealous P-Adic Numbers");
     Console.WriteLine(new PadicZealous(5, 3));
     Console.WriteLine();

     Console.WriteLine("Serie a");
     var a1 = new PadicZealous(5, 17, 86);
     Console.WriteLine(a1);
     Console.WriteLine(a1.Inv());
     Console.WriteLine(a1.Inv() * a1);
     Console.WriteLine(a1 * 5.Pow(2));
     Console.WriteLine(a1 / 5.Pow(2));
     Console.WriteLine(a1 * 5.Pow(5)); // one incorrect digit, 86~O(5^2) and 5^5~O(5^5) and the limit is O(5^7)
     Console.WriteLine(a1 / 5.Pow(5)); // one incorrect digit
     var a2 = new PadicZealous(5, 3, 75);
     Console.WriteLine(a2);
     Console.WriteLine(a2 * 5);
     Console.WriteLine(a2 / 5);
     var a3 = new PadicZealous(5, 3, 15);
     Console.WriteLine(a3);
     Console.WriteLine(a3 * 5);
     Console.WriteLine(a3 / 5);
     Console.WriteLine();

     Console.WriteLine("Serie b");
     var b1 = new PadicZealous(5, 3, 171);
     var b2 = new PadicZealous(5, 3, -171);
     var b3 = new PadicZealous(5, 6, 171);
     var b4 = new PadicZealous(5, 6, -171);
     Console.WriteLine(b1);
     Console.WriteLine(b2);
     Console.WriteLine(b2.Equals(b1.Opp()));
     Console.WriteLine(b3);
     Console.WriteLine(b4);
     Console.WriteLine(b3.Equals(b4.Opp()));
     Console.WriteLine(b1 * 5.Pow(4));
     Console.WriteLine(b1 / 5.Pow(4));
     Console.WriteLine();
     
     Console.WriteLine("Serie c");
     var c1 = new PadicZealous(3, 3, 101);
     var c2 = new PadicZealous(3, 6, 101);
     var c3 = new PadicZealous(3, 4, -101);
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
     Console.WriteLine();
     
     Console.WriteLine("Serie d");
     var d1 = new PadicZealous(3, 7, new Rational(101));
     var d2 = new PadicZealous(3, 7, new Rational(909));
     var d3 = new PadicZealous(3, 7, new Rational(101, 9));
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
     
     Console.WriteLine("Serie e");
     var e1 = new PadicZealous(3, 7, new Rational(909) + new Rational(101, 9));
     var e2 = new PadicZealous(3, 11, new Rational(909) + new Rational(101, 9));
     var e3 = new PadicZealous(5, 17, new Rational(313, 3000));
     Console.WriteLine(e1);
     Console.WriteLine(d1 / 9 + d2);
     Console.WriteLine(e1 * 3);
     Console.WriteLine(e1 * 9);
     Console.WriteLine(e1 / 3);
     Console.WriteLine(e1 / 9);
     Console.WriteLine(e2 * 3.Pow(4));
     Console.WriteLine(e2 / 3.Pow(4));
     Console.WriteLine(e3);
     Console.WriteLine(e3 * 5.Pow(2));
     Console.WriteLine(e3 / 5.Pow(2));
     Console.WriteLine(e1.Equals(d1 / 9 + d2));
     Console.WriteLine(e1.Equals(d2 + d3));
     Console.WriteLine();

    Console.WriteLine("Serie f");
    var f1 = new PadicZealous(5, 17, new Rational(313, 75)) * 5.Pow(4);
    var f2 = new PadicZealous(5, 17, 415) / 5.Pow(2);
    Console.WriteLine(f1);
    Console.WriteLine(f2);
    Console.WriteLine(f1 + f2);
    Console.WriteLine((f1 + f2) * 125);
    Console.WriteLine(f1 * 125);
    Console.WriteLine((f1 + f2) * 3125);
    Console.WriteLine(f1 * 3125);
    Console.WriteLine(f1 * 3125 + f2);

    var f3 = f1 * 125 + f2 * 125;
    var f4 = (f1 + f2) * 125;
    Console.WriteLine(f3);
    Console.WriteLine(f4);
    Console.WriteLine(f3.PNVSM.Equals(f4.PNVSM));
    Console.WriteLine((f3 - f4));
    Console.WriteLine((f3 - f4).Norm.IsZero());
    Console.WriteLine();

    Console.WriteLine(f1);
    Console.WriteLine(f1.Norm);
    Console.WriteLine(f1.Normalized);
    Console.WriteLine(f1.Normalized / f1.Norm);
    Console.WriteLine(f1 * f1.Norm);
    Console.WriteLine(f1.Normalized.Equals(f1 * f1.Norm));
    Console.WriteLine();
}

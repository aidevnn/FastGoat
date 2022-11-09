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
    var Q1 = new Rational();
    var x = Ring.Polynomial(new Rational(0, 1));
    Console.WriteLine(x);
    Console.WriteLine(x + 4);
    Console.WriteLine(-5 + 3 * x + 4);
    Console.WriteLine((x + 1) * (x - 3) / 5);
    Console.WriteLine((x + 1).Pow(2).Div(x - 1));

    var p1 = (x + 1) * (x + 3) * (x - 7);
    var p2 = (x * x + 5) * (x + 1);
    Console.WriteLine(p1);
    Console.WriteLine(p2);
    var (x0, y0) = Ring.Bezout(p1, p2);
    Console.WriteLine(new { x0, y0 });
    Console.WriteLine(p1 * x0 + p2 * y0);
    Console.WriteLine(Ring.Gcd(p1, p2));
}

{
    var x = Ring.Polynomial(new ZnInt(17, 0));
    Console.WriteLine(x);
    Console.WriteLine(x + 4);
    Console.WriteLine(-5 + 3 * x + 4);
    Console.WriteLine((x + 1) * (x - 3) / 5);
    Console.WriteLine((x + 1).Pow(2).Div(x - 1));

    var p1 = (x + 1) * (x + 3) * (x - 7);
    var p2 = (x * x + 5) * (x + 1);
    Console.WriteLine(p1);
    Console.WriteLine(p2);
    var (x0, y0) = Ring.Bezout(p1, p2);
    Console.WriteLine(new { x0, y0 });
    Console.WriteLine(p1 * x0 + p2 * y0);
    Console.WriteLine(Ring.Gcd(p1, p2));
}

{
    var Q1 = new Rational();
    var (x, y) = Ring.Polynomial('X', 'Y', Q1.Zero);
    Console.WriteLine(x);
    Console.WriteLine(y);
    Console.WriteLine(x + y);
    Console.WriteLine((x + y).Pow(2));
    Console.WriteLine((x + y).Pow(2).Div(x + 1));

    var a0 = (x + y).Pow(2);
    var (q, r) = a0.Div(x + 1);
    var a1 = q * (x + 1) + r;
    Console.WriteLine(a0);
    Console.WriteLine(new { q, r });
    Console.WriteLine(a1);
    Console.WriteLine(a1.Equals(a0));
    Console.WriteLine((x + y).Pow(3));
    Console.WriteLine((x + y).Pow(3).Div(x + 2));
    Console.WriteLine((x + y).Pow(3).Div(x * y)); // Throw exception 
    
}

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
    var (x, y) = Ring.Polynomial('X', 'Y', ZnInt.KZero());
    var a = Ring.Polynomial('a', ZnInt.KZero());
    Console.WriteLine(new { x, y, a });
    Console.WriteLine(new[] { x, y, a }.SelectMany(e => e.Indeterminates).Glue(" "));
    Console.WriteLine((x + a).Pow(2));
    Console.WriteLine((y + 1 + a).Pow(2));

    var p = (y + 1 + a).Pow(2) * (x - 4);
    Console.WriteLine(p);
    Console.WriteLine(p.Substitue(y, y.Zero));
    Console.WriteLine((x + a + 3).Pow(2).Substitue(x, x.Zero));
    Console.WriteLine((x + a + 3).Pow(2));
    Console.WriteLine((x + a + 3).Pow(2).Substitue(x, a - 1));
    Console.WriteLine((x + a + 3).Pow(2).Substitue(x, a));
}
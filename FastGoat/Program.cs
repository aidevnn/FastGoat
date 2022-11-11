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
    var x = Ring.Polynomial(ZnInt.KZero());
    var f = x.Pow(3) - x * x + 2 * x - 3;
    var y = Ring.Polynomial('Y', ZnInt.KZero());
    var g = y * y + y + 1;
    
    Console.WriteLine(f);
    Console.WriteLine(f.D(x));
    
    Console.WriteLine(g);
    Console.WriteLine(g.D(y));
}

{
    var x = Ring.Polynomial(ZnInt.KZero());
    var X = x.Indeterminates.First();
    var f = 2 * x * x + 3 * x + 4;
    var g = 5 * x + 6;
    var S = Ring.SylvesterMatrix(f, X, g, X);
    Console.WriteLine("Det = {0}", Ring.Determinant(S, f.Zero));
    Ring.DisplayMatrix(S);
}

{
    var x = Ring.Polynomial(ZnInt.KZero());
    var X = x.Indeterminates.First();
    var f = 2 * x * x + 3 * x + 4;
    var S = Ring.SylvesterMatrix(f, X, f, X);
    Console.WriteLine("Det = {0}", Ring.Determinant(S, f.Zero));
    Ring.DisplayMatrix(S);
}

{
    var x = Ring.Polynomial(ZnInt.KZero());
    var X = x.Indeterminates.First();
    var g = 5 * x + 6;
    var S = Ring.SylvesterMatrix(g, X, g, X);
    Console.WriteLine("Det = {0}", Ring.Determinant(S, g.Zero));
    Ring.DisplayMatrix(S);
}

{
    var x = Ring.Polynomial(ZnInt.KZero());
    var X = x.Indeterminates.First();
    var f = x.Pow(3) + 2 * x * x + 3 * x + 4;
    var S = Ring.SylvesterMatrix(f, X, f, X);
    Console.WriteLine("Det = {0}", Ring.Determinant(S, f.Zero));
    Ring.DisplayMatrix(S);
}

{
    var x = Ring.Polynomial(ZnInt.KZero());
    var X = x.Indeterminates.First();
    var f = x.Pow(3) + 2 * x * x + 3 * x + 4;
    var S = Ring.SylvesterMatrix(f, X, f + 5 * x.Pow(4), X);
    Console.WriteLine("Det = {0}", Ring.Determinant(S, f.Zero));
    Ring.DisplayMatrix(S);
}

{
    var x = Ring.Polynomial(ZnInt.KZero());
    var X = x.Indeterminates.First();
    var f = 2 * x.Pow(2) + 3 * x + 1;
    var g = 7 * x.Pow(2) + x + 3;
    var S = Ring.SylvesterMatrix(f, X, g, X);
    Console.WriteLine("Det = {0}", Ring.Determinant(S, f.Zero));
    Ring.DisplayMatrix(S);
}

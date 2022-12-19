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
using FastGoat.UserGroup.Padic;

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
     var a1 = new PadicZealous(5, 7, 86);
     Console.WriteLine(a1);
     Console.WriteLine(a1.Inv());
     Console.WriteLine(a1.Inv() * a1);
     Console.WriteLine(a1 * 5.Pow(2));
     Console.WriteLine(a1 / 5.Pow(2));
     Console.WriteLine(a1 * 5.Pow(5)); // one incorrect digit, 86~O(5^2) and 5^5~O(5^5) and the limit is O(5^7)
     Console.WriteLine(a1 / 5.Pow(5)); // one incorrect digit

     Console.WriteLine(86 * a1.Ppow(5));
     Console.WriteLine(86 / a1.Ppow(5));
}

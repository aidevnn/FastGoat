using System.Numerics;
using System.Text;
using Examples;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.EllCurve;
using FastGoat.UserGroup.Floats;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Polynoms;
using static FastGoat.Commons.IntExt;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");
RecomputeAllPrimesUpTo(32000000);

{
    Console.WriteLine(NextPrimes(2, 20).Glue(", "));
    Console.WriteLine(NextPrimes(1000, 10).Glue(", "));
    Console.WriteLine(NextPrime(BigInteger.Pow(10, 6)));
    Console.WriteLine(NextPrime(BigInteger.Pow(10, 10)));
    Console.WriteLine(NextPrime(BigInteger.Pow(10, 14)));
    Console.WriteLine(NextPrime(BigInteger.Pow(10, 15)));
}
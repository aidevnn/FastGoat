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

    Console.WriteLine("Comparison of two implementations of Padic");
    Console.WriteLine("Integer digits fixed and floating digits non fixed");
    {
        var a = new Padic(5, 7, 5);
        var a0 = a;
        for (int i = 0; i < 10; i++)
        {
            Console.WriteLine($"{a0}");
            a0 = a.Add(Padic.Convert(5, 7, new Rational(51, BigInteger.Pow(5, i))));
        }
    }
    Console.WriteLine();
    
    Console.WriteLine("Integer digits fixed and floating digits non fixed");
    {
        var a = new Padic(5, 7, 25);
        var a0 = a;
        for (int i = 0; i < 10; i++)
        {
            Console.WriteLine($"{a0}");
            a0 = a.Add(a0.CompleteDivision(a));
        }
    }
    Console.WriteLine();


    Console.WriteLine("Integer and floating digits fixed");
    {
        var a = new QpAdic(5, 7, 5);
        var a0 = a;
        for (int i = 0; i < 10; i++)
        {
            Console.WriteLine($"{a0}");
            a0 = a + (new QpAdic(5, 7, new Rational(51, BigInteger.Pow(5, i))));
        }
    }
    Console.WriteLine();
    
    Console.WriteLine("Integer and floating digits fixed");
    {
        var a = new QpAdic(5, 7, 25);
        var a0 = a;
        for (int i = 0; i < 10; i++)
        {
            Console.WriteLine($"{a0}");
            a0 = a + a0 / a;
        }
    }
    Console.WriteLine();
}
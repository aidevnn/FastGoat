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
    Console.WriteLine("Inverse floating padic");
    {
        var a = Padic.Convert(2, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1);
        var a0 = a.Add(a.Ppow(2));
        Console.WriteLine($"c = {a}");
        Console.WriteLine($"a0 = {a0}");
        for (int i = 0; i < 7; i++)
        {
            Console.WriteLine($"Step {i} : {a0}");
            a0 = a0.Mul(2).Sub(a.Mul(a0.Pow(2)));
        }
        Console.WriteLine($"Check a^-1 = {a.Inv()} => {a0.Inv().Equals(a)}");
    }
    Console.WriteLine();

    Console.WriteLine("Inverse zealous");
    {
        var a = QpAdic.Convert(2, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1);
        var a0 = a + a.Ppow(2);
        Console.WriteLine($"c = {a}");
        Console.WriteLine($"a0 = {a0}");
        for (int i = 0; i < 7; i++)
        {
            Console.WriteLine($"Step {i} : {a0}");
            a0 = 2 * a0 - a * a0.Pow(2);
        }
        Console.WriteLine($"Check a^-1 = {a.Inv()} => {a0.Inv().Equals(a)}");
    }
    Console.WriteLine();
}

{
    Console.WriteLine("Square Root floating padic");
    {
        var a = Padic.Convert(2, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1);
        var a0 = a.Add(a.Ppow(2));
        Console.WriteLine($"c = {a}");
        Console.WriteLine($"a0 = {a0}");
        for (int i = 0; i < 7; i++)
        {
            Console.WriteLine($"Step {i} : {a0}");
            a0 = a0.Add(a.CompleteDivision(a0)).CompleteDivision(a.One.Mul(2));
        }
        Console.WriteLine($"Check ai^2 = {a0.Pow(2)} => {a0.Pow(2).Equals(a)}");
    }
    Console.WriteLine();

    Console.WriteLine("Square Root zealous");
    {
        var a = QpAdic.Convert(2, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1);
        var a0 = a + a.Ppow(2);
        Console.WriteLine($"c = {a}");
        Console.WriteLine($"a0 = {a0}");
        for (int i = 0; i < 7; i++)
        {
            Console.WriteLine($"Step {i} : {a0}");
            a0 = (a0 + a / a0) / 2;
        }
        Console.WriteLine($"Check ai^2 = {a0.Pow(2)} => {a0.Pow(2).Equals(a)}");
    }
    Console.WriteLine();
}

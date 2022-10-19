using System.Collections;
using FastGoat;
using FastGoat.Examples;
using FastGoat.Gp;
using FastGoat.UserGroup;
using static FastGoat.IntExt;
using static FastGoat.EnumerableExt;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");

{
    foreach (var e in YieldAllCombinations(4))
        Console.WriteLine(e.Select(b => b ? 1 : 0).Glue());

    Console.WriteLine();
    
    foreach (var e in YieldAllPermutations(4))
        Console.WriteLine(e.Glue());
}
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
    var l0 = new[] { 1, 2, 3 };
    var l1 = new[] { 4, 5 };

    foreach (var e in l0.MultiLoopWith(l1, l0))
    {
        Console.WriteLine(e.Glue(" "));
    }
}
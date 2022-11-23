using System.Diagnostics;
using FastGoat.Commons;
using FastGoat.Examples;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Matrix;
using FastGoat.UserGroup.Perms;
using FastGoat.UserGroup.Polynoms;
using FastGoat.UserGroup.Words;
using FastGoat.UserGroup.Words.ToddCoxeter;
using static FastGoat.Commons.IntExt;
using static FastGoat.Commons.EnumerableExt;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");

{
    for (int i = 2; i < 51; i++)
    {
        var dicn = Group.DiCyclicSdp(i);
        DisplayGroup.HeadOrders(dicn);
    }
}

{
    GlobalStopWatch.Restart();
    for (int i = 2; i < 51; i++)
    {
        Console.WriteLine($"Go{i}");
        var g1 = Group.DiCyclicSdp(i); // dynamic type
        var g2 = Group.DiCyclic(i);
        if (g2.IsIsomorphicTo(g1))
            continue;
        
        DisplayGroup.AreIsomorphics(g2, g1); 
        DisplayGroup.HeadOrders(g1);
        DisplayGroup.HeadOrders(g2);
    }
    GlobalStopWatch.Show("DiCyclic Order <= 200"); // Time:133652 ms
}
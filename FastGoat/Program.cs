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
    for (int i = 2; i < 51; i++)
    {
        Console.WriteLine($"Go{i}");
        var g1 = Group.DiCyclicSdp(i);
        var g2 = Group.DiCyclic(i);
        if (g1.IsIsomorphicTo(g2))
            continue;
        
        DisplayGroup.AreIsomorphics(g2, g1); 
        DisplayGroup.HeadOrders(g1);
        DisplayGroup.HeadOrders(g2);
        // Dic30 && Dic42 errors
    }
}
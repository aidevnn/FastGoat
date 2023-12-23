using System.Collections;
using System.ComponentModel;
using System.Diagnostics;
using System.Numerics;
using FastGoat.Commons;
using FastGoat.Examples;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.Naming;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Characters;
using FastGoat.UserGroup.GModuleN;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Matrix;
using FastGoat.UserGroup.Perms;
using FastGoat.UserGroup.Polynoms;
using FastGoat.UserGroup.Words;
using FastGoat.UserGroup.Words.ToddCoxeter;
using static FastGoat.Commons.IntExt;
using static FastGoat.Commons.EnumerableExt;
using FastGoat.UserGroup.Padic;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");

void Ord16()
{
    GlobalStopWatch.Restart();
    FG.AllExtensions(
            (FG.Abelian(8), FG.Abelian(2)),
            (FG.Abelian(2, 4), FG.Abelian(2)),
            (FG.Abelian(2, 2, 2), FG.Abelian(2)))
        .Take(14)
        .NamingExts()
        .DisplayExts();
    
    GlobalStopWatch.Show();
    Console.Beep();
}

void Ord20_40()
{
    GlobalStopWatch.Restart();
    var all20 = FG.AllExtensions(
            (FG.Abelian(5), FG.Abelian(4)),
            (FG.Abelian(2, 5), FG.Abelian(2)))
        .Take(5)
        .NamingExts()
        .ToArray();

    var all40 = FG.AllExtensions(
            (FG.Abelian(2, 5), FG.Abelian(4)),
            (FG.Abelian(2, 5), FG.Abelian(2, 2)))
        .Take(14)
        .NamingExts()
        .ToArray();
    
    all20.DisplayExts();
    all40.DisplayExts();

    GlobalStopWatch.Show();
    Console.Beep();
}

void Ord32()
{
    GlobalStopWatch.Restart();
    var nbOpsMax = 10000;
    FG.AllExtensions(
            (nbOpsMax, FG.Abelian(16), FG.Abelian(2)),
            (nbOpsMax, FG.Abelian(2, 8), FG.Abelian(2)),
            (nbOpsMax, FG.Abelian(4, 4), FG.Abelian(2)),
            (nbOpsMax, FG.Abelian(2, 2, 4), FG.Abelian(2)),
            (nbOpsMax, FG.Abelian(2, 4), FG.Abelian(4)),
            (16, FG.Abelian(2, 4), FG.Abelian(2, 2)),
            (1, FG.Abelian(2, 2, 2), FG.Abelian(2, 2)))
        .Take(51)
        .NamingExts()
        .DisplayExts();

    GlobalStopWatch.Show();
    Console.Beep();
}

{
    // Ord16();
    // Ord32();
    Ord20_40();
}
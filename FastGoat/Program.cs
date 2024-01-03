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

{
    GlobalStopWatch.Restart();
    var ord24 = FG.AllExtensions(
            (FG.Abelian(12).ToCGTable(), FG.Abelian(2).ToCGTable()),
            (FG.Abelian(2, 6).ToCGTable(), FG.Abelian(2).ToCGTable()),
            (FG.Alternate(4).ToCGTable(), FG.Abelian(2).ToCGTable()),
            (FG.Quaternion(8).ToCGTable(), FG.Abelian(3).ToCGTable())
        )
        .Select(e => e.allSubs.ToTable())
        .FilterIsomorphic()
        .Naming()
        .ToArray();

    var listByC2 = ord24.SelectMany(e => FG.AllExtensions((e.subsg.Parent, FG.Abelian(2)))).Select(e => e.allSubs.ToTable());
    var listByC3 = FG.AllAbelianGroupsOfOrder(16).SelectMany(e => Group.AllSDPFilterLazy(e, FG.Abelian(3)));
    listByC2.AppendIsomorphic(listByC3)
        .Take(52)
        .Naming()
        .DisplayNames();

    GlobalStopWatch.Show();
    Console.Beep();
    // Total Groups:52
    // #  Time:1m16s
    // 
}

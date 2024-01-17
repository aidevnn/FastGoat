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

void AllGroupsOrderLess32()
{
    GlobalStopWatch.Restart();
    var allAb = 32.Range(1).SelectMany(k => FG.AllAbelianGroupsOfOrder(k)).Select(e => e.ToCGW().AllSubgroups());
    var allMtCycSdp = 32.Range(1).SelectMany(k => FG.MetaCyclicSdp(k)).Select(e => e.ToCGW().AllSubgroups());
    var ext8 = FG.AllExtensions((FG.Abelian(2), FG.Abelian(2, 2))).Select(e => e.allSubs.ToGroupWrapper());
    var ext12 = FG.AllExtensions((FG.Abelian(2, 2), FG.Abelian(3))).Select(e => e.allSubs.ToGroupWrapper());
    var ext16 = FG.AllExtensions((FG.Abelian(8), FG.Abelian(2)), (FG.Abelian(4, 2), FG.Abelian(2)))
        .Select(e => e.allSubs.ToGroupWrapper());
    var ext18 = FG.AllExtensions((FG.Abelian(3, 3), FG.Abelian(2))).Select(e => e.allSubs.ToGroupWrapper());
    var ext24 = FG.AllExtensions(
            (FG.Abelian(12).ToCGW(), FG.Abelian(2)),
            (FG.Abelian(2, 6).ToCGW(), FG.Abelian(2)),
            (FG.Alternate(4).ToCGW(), FG.Abelian(2)),
            (FG.Quaternion(8).ToCGW(), FG.Abelian(3)))
        .Select(e => e.allSubs.ToGroupWrapper());
    var ext27 = FG.AllExtensions((FG.Abelian(3, 3), FG.Abelian(3))).Select(e => e.allSubs.ToGroupWrapper());
    var ext32 = FG.AllExtensions(
            (FG.Abelian(16), FG.Abelian(2)),
            (FG.Abelian(8, 2), FG.Abelian(2)),
            (FG.Abelian(4, 4), FG.Abelian(2)),
            (FG.Abelian(4, 2), FG.Abelian(4)),
            (FG.Abelian(4, 2), FG.Abelian(2, 2)))
        .Select(e => e.allSubs.ToGroupWrapper());

    allAb.AppendIsomorphic(allMtCycSdp, ext8, ext12, ext16, ext18, ext27, ext24, ext32)
        .Naming()
        .DisplayNames();

    GlobalStopWatch.Show();
    Console.Beep();
    
    // Total Groups:144
    // #  Time:25.414s
    // 
}

{
    AllGroupsOrderLess32();
}
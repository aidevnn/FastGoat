using System.Collections;
using System.ComponentModel;
using System.ComponentModel.DataAnnotations;
using System.Diagnostics;
using System.Numerics;
using FastGoat.Commons;
using FastGoat.Examples;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.Naming;
using FastGoat.Structures.Subgroups;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Characters;
using FastGoat.UserGroup.GModuleN;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Matrix;
using FastGoat.UserGroup.Perms;
using FastGoat.UserGroup.Polynoms;
using FastGoat.UserGroup.Words;
using static FastGoat.Commons.IntExt;
using static FastGoat.Commons.EnumerableExt;
using FastGoat.UserGroup.Padic;
using FastGoat.UserGroup.Words.Tools;
using GroupRegX = System.Text.RegularExpressions.Group;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");

void SL25ExtensionsByC2()
{
    var sl25 = FG.SL2p(5);

    Logger.Level = LogLevel.Level1;
    FG.AllExtensions((sl25, FG.Abelian(2)))
        .Select(e => e.allSubs.ToGroupWrapper())
        .FilterIsomorphic()
        .Naming()
        .DisplayNames();
}

void Ord48()
{
    Logger.Level = LogLevel.Level1;
    GlobalStopWatch.Restart();
    var ord24 = FG.AllExtensions(
            (FG.Abelian(12).ToCGW(), FG.Abelian(2)),
            (FG.Abelian(2, 6).ToCGW(), FG.Abelian(2)),
            (FG.Alternate(4).ToCGW(), FG.Abelian(2)),
            (FG.Quaternion(8).ToCGW(), FG.Abelian(3))
        )
        .Select(e => e.allSubs.ToGroupWrapper())
        .FilterIsomorphic()
        .Take(GroupExt.A000001[24])
        .Naming()
        .ToArray();

    var listByC2 = ord24.SelectMany(e => FG.AllExtensions((e.subsg.Parent, FG.Abelian(2)))).Select(e => e.allSubs.ToGroupWrapper());
    var listByC3 = FG.AllAbelianGroupsOfOrder(16)
        .SelectMany(e => FG.AllSDPFilterLazy(e, FG.Abelian(3)))
        .Select(e => e.AllSubgroups().ToGroupWrapper());
    listByC2.AppendIsomorphic(listByC3)
        .Take(GroupExt.A000001[48])
        .Naming()
        .DisplayNames();

    GlobalStopWatch.Show();
    Console.Beep();
    // Total Groups:52
}

{
    // SL25ExtensionsByC2();
    Ord48();
}
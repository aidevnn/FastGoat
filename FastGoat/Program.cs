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

void Ord36()
{
    Logger.Level = LogLevel.Level1;
    GlobalStopWatch.Restart();
    
    FG.AllExtensions(
            (FG.Abelian(9), FG.Abelian(4)),
            (FG.Abelian(9), FG.Abelian(2, 2)),
            (FG.Abelian(3, 3), FG.Abelian(4)),
            (FG.Abelian(3, 3), FG.Abelian(2, 2)),
            (FG.Abelian(2, 6), FG.Abelian(3))
        )
        .Select(e => e.allSubs.ToGroupWrapper())
        .FilterIsomorphic()
        .Take(GroupExt.A000001[36])
        .Naming()
        .DisplayNames();

    GlobalStopWatch.Show("End");
    Console.Beep();
}

void Ord40()
{
    Logger.Level = LogLevel.Level1;
    GlobalStopWatch.Restart();
    
    FG.AllExtensions(
            (FG.Abelian(10), FG.Abelian(4)),
            (FG.Abelian(10), FG.Abelian(2, 2))
        )
        .Select(e => e.allSubs.ToGroupWrapper())
        .FilterIsomorphic()
        .Take(GroupExt.A000001[40])
        .Naming()
        .DisplayNames();

    GlobalStopWatch.Show("End");
    Console.Beep();
}

void Ord54()
{
    Logger.Level = LogLevel.Level1;
    GlobalStopWatch.Restart();
    
    var ord = 54;
    var allAb = FG.AllAbelianGroupsOfOrder(ord).Select(g => g.AllSubgroups().ToGroupWrapper());
    var allMetaCyc = FG.MetaCyclicSdp(ord).Select(g => g.AllSubgroups().ToGroupWrapper());
    var allExt = FG.AllExtensions(
            (FG.Abelian(3, 6), FG.AbelianPerm(3)),
            (FG.Abelian(3, 3), FG.Symmetric(3))
        )
        .Select(e => e.allSubs.ToGroupWrapper());
    
    allAb.AppendIsomorphic(allMetaCyc, allExt)
        .Take(GroupExt.A000001[ord])
        .Naming()
        .DisplayNames();

    GlobalStopWatch.Show("End");
    Console.Beep();
}

void Ord56()
{
    Logger.Level = LogLevel.Level1;
    GlobalStopWatch.Restart();

    var ord = 56;
    var allAb = FG.AllAbelianGroupsOfOrder(ord).Select(g => g.AllSubgroups().ToGroupWrapper());
    var allMetaCyc = FG.MetaCyclicSdp(ord).Select(g => g.AllSubgroups().ToGroupWrapper());
    var sdp = FG.AllSDPFilter(FG.Abelian(2, 2, 2), FG.Abelian(7)).Select(g => g.AllSubgroups().ToGroupWrapper());
    var allExt = FG.AllExtensions(
            (FG.Abelian(28), FG.Abelian(2)),
            (FG.Abelian(2, 14), FG.Abelian(2))
        )
        .Select(e => e.allSubs.ToGroupWrapper());
    
    allAb.AppendIsomorphic(allMetaCyc, sdp, allExt)
        .Take(GroupExt.A000001[ord])
        .Naming()
        .DisplayNames();

    GlobalStopWatch.Show("End");
    Console.Beep();
}

void Ord60()
{
    Logger.Level = LogLevel.Level1;
    GlobalStopWatch.Restart();

    var ord = 60;
    var allAb = FG.AllAbelianGroupsOfOrder(ord).Select(g => g.AllSubgroups().ToGroupWrapper());
    var allMetaCyc = FG.MetaCyclicSdp(ord).Select(g => g.AllSubgroups().ToGroupWrapper());
    var a5 = FG.Alternate(5).AllSubgroups().ToGroupWrapper();
    var allExt = FG.AllExtensions(
            (FG.Abelian(30), FG.Abelian(2)),
            (FG.Abelian(15), FG.Abelian(2, 2)),
            (FG.Abelian(2, 10), FG.Abelian(3))
        )
        .Select(e => e.allSubs.ToGroupWrapper());
    
    allAb.AppendIsomorphic(allMetaCyc, [a5], allExt)
        .Take(GroupExt.A000001[ord])
        .Naming()
        .DisplayNames();

    GlobalStopWatch.Show("End");
    Console.Beep();
}

IEnumerable<(AllSubgroups<WElt> subsg, ANameElt[] names)> GetAllProds(int ord)
{
    var allAb = FG.AllAbelianGroupsOfOrder(ord).Select(g => g.AllSubgroups().ToGroupWrapper());
    var allMetaCyc = FG.MetaCyclicSdp(ord).Select(g => g.AllSubgroups().ToGroupWrapper());
    return allAb.AppendIsomorphic(allMetaCyc)
        .Take(GroupExt.A000001[ord])
        .Naming();
}

{
    // SL25ExtensionsByC2();
    // Ord48();
    // Ord80();

    Logger.Level = LogLevel.Level1;
    GlobalStopWatch.Restart();

    31.Range(33).SelectMany(o => GetAllProds(o))
        .DisplayNames()
        .CheckMissings();

    GlobalStopWatch.Show("End");
    Console.Beep();
}
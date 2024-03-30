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

IEnumerable<AllSubgroups<WElt>> Ord48()
{
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

    var ord = 48;
    return listByC2.AppendIsomorphic(listByC3, GetAllProds(ord)).Take(GroupExt.A000001[ord]);
}

IEnumerable<AllSubgroups<WElt>> Ord36()
{
    var ord = 36;
    return FG.AllExtensions(
            (FG.Abelian(9), FG.Abelian(4)),
            (FG.Abelian(9), FG.Abelian(2, 2)),
            (FG.Abelian(3, 3), FG.Abelian(4)),
            (FG.Abelian(3, 3), FG.Abelian(2, 2)),
            (FG.Abelian(2, 6), FG.Abelian(3))
        )
        .Select(e => e.allSubs.ToGroupWrapper())
        .AppendIsomorphic(GetAllProds(ord)).Take(GroupExt.A000001[ord]);
}

IEnumerable<AllSubgroups<WElt>> Ord40()
{
    var ord = 48;
    return FG.AllExtensions(
            (FG.Abelian(10), FG.Abelian(4)),
            (FG.Abelian(10), FG.Abelian(2, 2))
        )
        .Select(e => e.allSubs.ToGroupWrapper())
        .AppendIsomorphic(GetAllProds(ord)).Take(GroupExt.A000001[ord]);
}

IEnumerable<AllSubgroups<WElt>> Ord50()
{
    var ord = 50;
    return FG.AllSDPFilter(FG.Abelian(5, 5), FG.Abelian(2))
        .Select(g => g.AllSubgroups().ToGroupWrapper())
        .AppendIsomorphic(GetAllProds(ord)).Take(GroupExt.A000001[ord]);
}

IEnumerable<AllSubgroups<WElt>> Ord54()
{
    var ord = 54;
    return FG.AllExtensions(
            (FG.Abelian(3, 6), FG.AbelianPerm(3)),
            (FG.Abelian(3, 3), FG.Symmetric(3))
        )
        .Select(e => e.allSubs.ToGroupWrapper())
        .AppendIsomorphic(GetAllProds(ord)).Take(GroupExt.A000001[ord]);
}

IEnumerable<AllSubgroups<WElt>> Ord56()
{
    var ord = 56;
    var sdp = FG.AllSDPFilter(FG.Abelian(2, 2, 2), FG.Abelian(7)).Select(g => g.AllSubgroups().ToGroupWrapper());
    return FG.AllExtensions(
            (FG.Abelian(28), FG.Abelian(2)),
            (FG.Abelian(2, 14), FG.Abelian(2))
        )
        .Select(e => e.allSubs.ToGroupWrapper())
        .AppendIsomorphic(sdp, GetAllProds(ord)).Take(GroupExt.A000001[ord]);
}

IEnumerable<AllSubgroups<WElt>> Ord60()
{
    var ord = 60;
    var a5 = FG.Alternate(5).AllSubgroups().ToGroupWrapper();
    var allExt = FG.AllExtensions(
            (FG.Abelian(30), FG.Abelian(2)),
            (FG.Abelian(15), FG.Abelian(2, 2)),
            (FG.Abelian(2, 10), FG.Abelian(3))
        )
        .Select(e => e.allSubs.ToGroupWrapper());

    return allExt.AppendIsomorphic([a5], GetAllProds(ord)).Take(GroupExt.A000001[ord]);
}

IEnumerable<AllSubgroups<WElt>> GetAllProds(int ord)
{
    var allAb = FG.AllAbelianGroupsOfOrder(ord).Select(g => g.AllSubgroups().ToGroupWrapper());
    var allMetaCyc = FG.MetaCyclicSdp(ord).Select(g => g.AllSubgroups().ToGroupWrapper());
    return allAb.Concat(allMetaCyc);
}

{
    // SL25ExtensionsByC2();
    // Ord48();
    // Ord80();

    Logger.Level = LogLevel.Level1;
    GlobalStopWatch.Restart();

    31.Range(33).SelectMany(o => GetAllProds(o))
        .FilterIsomorphic()
        .Concat([..Ord36(), ..Ord40(), ..Ord48(), ..Ord50(), ..Ord54(), ..Ord56(), ..Ord60()])
        .FilterIsomorphic()
        .Naming()
        .DisplayNames()
        .CheckMissings();

    GlobalStopWatch.Show("End");
    Console.Beep();
}
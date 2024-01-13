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

void ord64()
{
    GlobalStopWatch.Restart();

    var g3232 = FG.AllExtensions((FG.Abelian(4, 4), FG.Abelian(2))).Select(e => e.allSubs).Naming()
        .First(e => e.subsg.Infos.ToTuples() == (34, 28, 22)).subsg.Parent.ToCGW();

    var sdps = FG.AllSDPFilterLazy(FG.Abelian(4), FG.Abelian(8))
        .AppendIsomorphic(FG.AllSDPFilterLazy(FG.Abelian(8), FG.Abelian(4)))
        .Naming()
        .Select(e => (e.subsg.Parent, FG.Abelian(2)))
        .ToArray();
    var sdp1 = Group.AllSemiDirectProd(FG.Abelian(4, 2), FG.Abelian(4))
        .First(e => e.AllSubgroups().Infos.ToTuples() == (50, 38, 26))
        .ToCGW();
    var sdp2 = Product.Generate(FG.Abelian(2), FG.ModularMaxSdp(4)).ToCGW();
    var listByC2 = sdps.Concat([
        (sdp1, FG.Abelian(2)),
        (sdp2, FG.Abelian(2)),
        (g3232, FG.Abelian(2)),
        (FG.Abelian(8, 4).ToCGW(), FG.Abelian(2)),
        (FG.Abelian(4, 4, 2).ToCGW(), FG.Abelian(2)),
        (FG.Abelian(32).ToCGW(), FG.Abelian(2))
    ]);
    var listByC4_C2C2 = new[]
    {
        (FG.Abelian(8, 2).ToCGW(), FG.Abelian(4)),
        (FG.Abelian(4, 4).ToCGW(), FG.Abelian(4)),
        (FG.Abelian(4, 2, 2).ToCGW(), FG.Abelian(4)),
        (FG.Abelian(8, 2).ToCGW(), FG.Abelian(2, 2))
    };

    var ord8 = new[]
    {
        FG.Abelian(8).ToCGW(),
        FG.Abelian(4, 2).ToCGW(),
        FG.ElementaryAbelian(8).ToCGW(),
        FG.DihedralSdp(4).ToCGW(),
        FG.Quaternion(8).ToCGW()
    };

    var g1612 = Product.Generate(FG.Quaternion(8), FG.Abelian(2));
    var listSdp64a = FG.MetaCyclicSdp(64).Select(e => e.ToGroupWrapper().AllSubgroups())
        .Concat(ord8.Grid2D(ord8).Select(e => Product.Generate(e.t1, e.t2).ToGroupWrapper().AllSubgroups()))
        .Concat(FG.AllSDPFilter(g1612, FG.Abelian(2, 2)))
        .Concat(FG.AllSDPFilter(FG.Abelian(4, 2), FG.DihedralSdp(4)))
        .Concat(FG.AllSDPFilter(FG.Abelian(4, 2, 2), FG.Abelian(2, 2)));

    var listSdp64b = FG.AllSDPFilter(FG.Abelian(4, 4), FG.Abelian(2))
        .Concat(FG.AllSDPFilter(FG.ModularMaxSdp(4), FG.Abelian(2)))
        .Append(Product.Generate(FG.Abelian(4), FG.Quaternion(8)).AllSubgroups().ToGroupWrapper())
        .FilterIsomorphic().Naming().ToArray()
        .SelectMany(e => FG.AllSDPFilter(e.subsg.Parent, FG.Abelian(2), trivial: true));

    var listAb64 = FG.AllAbelianGroupsOfOrder(64).Select(e => e.ToCGW().AllSubgroups());
    var listExts = FG.AllExtensions([..listByC2, ..listByC4_C2C2]).Select(e => e.allSubs.ToGroupWrapper());

    var allOrd64 = listAb64
        .Concat(listSdp64a)
        .Concat(listSdp64b)
        .Concat(listExts)
        .FilterIsomorphic()
        .Take(GroupExt.A000001[64])
        .Naming()
        .ToArray();

    allOrd64.DisplayNames();

    var listIds = FG.AllIds(64).ToList();
    foreach (var sub in allOrd64)
    {
        var id = listIds.Find(e => e.Infos == sub.subsg.Infos);
        listIds.Remove(id);
    }

    var pos = FG.AllIds(64).Where(e => listIds.Any(f => e.Infos == f.Infos)).ToList();
    pos.Println($"Remaining {listIds.Count} groups, possibles {pos.Count}");

    GlobalStopWatch.Show("End Naming");
    Console.Beep();
}

// Remaining 0 groups, possibles 0
// 
// # End Naming Time:16m57s

void ord48()
{
    GlobalStopWatch.Restart();
    var ord24 = FG.AllExtensions(
            (FG.Abelian(12).ToCGW(), FG.Abelian(2)),
            (FG.Abelian(2, 6).ToCGW(), FG.Abelian(2)),
            (FG.Alternate(4).ToCGW(), FG.Abelian(2)),
            (FG.Quaternion(8).ToCGW(), FG.Abelian(3))
        )
        .Select(e => e.allSubs)
        .FilterIsomorphic()
        .Naming()
        .ToArray();

    var listByC2 = ord24.SelectMany(e => FG.AllExtensions((e.subsg.Parent, FG.Abelian(2)))).Select(e => e.allSubs.ToGroupWrapper());
    var listByC3 = FG.AllAbelianGroupsOfOrder(16).Where(e => Group.AbelianGroupType(e).Length < 4)
        .SelectMany(e => FG.AllSDPFilter(e, FG.Abelian(3)))
        .Concat(FG.AllSDPFilterLazy(FG.Abelian(2, 2), FG.Alternate(4)));
    listByC2.AppendIsomorphic(listByC3)
        .Take(GroupExt.A000001[48])
        .Naming()
        .DisplayNames();

    GlobalStopWatch.Show();
    Console.Beep();
    // Total Groups:52
    // #  Time:13.357s
    // 
}

void ord32()
{
    GlobalStopWatch.Restart();
    FG.AllExtensions(
            (FG.Abelian(16), FG.Abelian(2)),
            (FG.Abelian(8, 2), FG.Abelian(2)),
            (FG.Abelian(4, 4), FG.Abelian(2)),
            (FG.Abelian(4, 2), FG.Abelian(4)),
            (FG.Abelian(4, 2), FG.Abelian(2, 2)),
            (FG.Abelian(2, 2, 2), FG.Abelian(2, 2))
        )
        .Select(e => e.allSubs)
        .FilterIsomorphic()
        .Take(GroupExt.A000001[32])
        .Naming()
        .DisplayNames();

    GlobalStopWatch.Show();
    Console.Beep();
    // Total Groups:51
    // #  Time:19.091s
    // 
}

{
    // ord32();
    // ord32();
    // ord32();
    // ord48();
    // ord48();
    // ord48();
    ord64();
}
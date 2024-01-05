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

    var c2 = FG.Abelian(2);
    var sm3232 = FG.AllExtensions((FG.Abelian(4, 4), FG.Abelian(2))).Select(e => e.allSubs).Naming()
        .First(e => e.names[0].IsExtension).subsg.Parent.ToCGTable();

    var sdps = FG.AllSDPFilterLazy(FG.Abelian(4), FG.Abelian(8))
        .AppendIsomorphic(FG.AllSDPFilterLazy(FG.Abelian(8), FG.Abelian(4)))
        .Naming()
        .Select(e => (e.subsg.Parent, c2))
        .ToArray();
    var sdp1 = Group.AllSemiDirectProd(FG.Abelian(4, 2), FG.Abelian(4))
        .First(e => e.AllSubgroups().Infos.ToTuples() == (50, 38, 26))
        .ToCGTable();
    var sdp2 = Product.Generate(c2, FG.ModularMaxSdp(4)).ToCGTable();
    var listByC2 = sdps.Concat([
        (sdp1, c2),
        (sdp2, c2),
        (sm3232, c2),
        (FG.Abelian(8, 4).ToCGTable(), c2),
        (FG.Abelian(32).ToCGTable(), c2)
    ]);
    var listByC4_C2C2 = new[]
    {
        (FG.Abelian(8, 2).ToCGTable(), FG.Abelian(4)),
        (FG.Abelian(4, 4).ToCGTable(), FG.Abelian(4)),
        (FG.Abelian(8, 2).ToCGTable(), FG.Abelian(2, 2))
    };

    var ord8 = new[]
    {
        FG.Abelian(8).ToCGTable(),
        FG.Abelian(4, 2).ToCGTable(),
        FG.ElementaryAbelian(8).ToCGTable(),
        FG.DihedralSdp(4).ToCGTable(),
        FG.Quaternion(8).ToCGTable()
    };

    var sdp32 = FG.AllSDPFilter(FG.Abelian(4, 2, 2), FG.Abelian(2)).Select(e => e.Parent).ToArray();
    var listSdp64a = sdp32.SelectMany(e => FG.AllSDPFilter(e, FG.Abelian(2), trivial: true))
        .Concat(FG.AllSDPFilter(FG.Abelian(4, 2, 2), FG.Abelian(4)))
        .Concat(ord8.Grid2D(ord8).SelectMany(e => FG.AllSDPFilter(e.t1, e.t2, trivial: true)));

    var listSdp64b = FG.AllExtensions(
            (FG.Abelian(2, 4), FG.Abelian(4)),
            (FG.Abelian(4, 4), FG.Abelian(2)))
        .Select(e => e.allSubs)
        .FilterIsomorphic()
        .Naming()
        .ToArray()
        .SelectMany(e => FG.AllSDPFilter(e.subsg.Parent, c2, trivial: true));

    var listAb64 = FG.AllAbelianGroupsOfOrder(64).Select(e => e.ToCGTable().AllSubgroups());
    var listExts = FG.AllExtensions([..listByC2, ..listByC4_C2C2]).Select(e => e.allSubs.ToTable());

    var allOrd64 = listAb64
        .Concat(listSdp64a)
        .Concat(listExts)
        .Concat(listSdp64b)
        .FilterIsomorphic()
        .Take(267);
        // .ToArray(); // Memory issue 

    var listIds = FG.AllIds(64).ToList();
    var nb = 0;
    foreach (var sub in allOrd64)
    {
        var id = listIds.Find(e => e.Infos == sub.Infos);
        listIds.Remove(id);
        FG.DisplayBox(sub, ++nb);
    }

    var pos = FG.AllIds(64).Where(e => listIds.Any(f => e.Infos == f.Infos)).ToList();
    pos.Println($"Remaining {listIds.Count} groups, possibles {pos.Count}");

    Console.WriteLine();
    Console.WriteLine($"Total Groups:{nb}");
    GlobalStopWatch.Show();
    Console.Beep();
}

// Remaining 0 groups, possibles 0
// 
// 
// Total Groups:267
// #  Time:20m53s
// 

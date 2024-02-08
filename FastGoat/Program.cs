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

void GroupDetails(AllSubgroups<WElt> subgroups, ANameElt[] names, bool rename = false, int maxLt = -1)
{
    var nbSharp = 16;
    var (g, infos) = (subgroups.Parent, subgroups.Infos);
    if (rename)
        g.Name = names[0].Name;

    var name = g.Name;
    maxLt = int.Max(name.Length, maxLt);
    var diff = (maxLt - name.Length) / 2;
    var space = Enumerable.Repeat(' ', diff).Glue();
    var lt = Enumerable.Repeat('#', maxLt + 4).Glue();
    var sharp = Enumerable.Repeat('#', nbSharp).Glue();
    var line = $"{sharp}{lt}{sharp}";
    var fmt = $"{sharp}{space}  {{0,{-maxLt + diff * 2}}}  {space}{sharp}";
    Console.WriteLine(line);
    Console.WriteLine(fmt, g.Name);
    Console.WriteLine(line);
    
    Console.WriteLine(g.ShortName);
    Console.WriteLine($"Type      {g.GroupType}");
    DisplayGroup.Orders(g);
    Console.WriteLine(infos);
    Console.WriteLine();
    
    var frattini = subgroups.FrattiniSubGroup;
    var fratName = NamesTree.BuildName(subgroups.Restriction(frattini))[0].Name;
    Console.WriteLine($"Frattini Φ(G) = {fratName}");
    
    var fitting = subgroups.FittingSubGroup;
    var fitName = NamesTree.BuildName(subgroups.Restriction(fitting))[0].Name;
    Console.WriteLine($"Fitting  F(G) = {fitName}");
    Console.WriteLine();

    var rels = Graph.DefiningRelatorsOfGroup(g, details: false);
    var gens = rels.Where(c => char.IsLetter(c)).Distinct().Select(c => char.ToLower(c)).Order().ToArray();
    var def = $"< {gens.Glue(",")} | {rels.Replace(" ", "").Replace(",", ", ").Replace("=", " = ")} >";

    Console.WriteLine("Word Group");
    Console.WriteLine(def);
    Console.WriteLine();
    
    var gapInfos = FG.FindIdGroup(g, infos);
    var s = gapInfos.Length > 1 ? " (TODO)" : "";
    foreach (var e in gapInfos)
        Console.WriteLine($"{$"Gap SmallGroup({e.Order},{e.No})",-24} Name:{e.Name}{s}");
    
    Console.WriteLine();
    names.Println("Group names");
    Console.WriteLine();

    Console.WriteLine("Characters Table");
    Console.WriteLine();
    if (g.Name == "SL(2,3)")
    {
        var gl23 = FG.GL2p(3);
        var ctGL23 = FG.CharacterTable(gl23);

        var a = gl23[1, 1, 0, 1];
        var b = gl23[0, 1, 2, 0];

        var sl23 = Group.Generate("SL(2,3)", gl23, a, b);
        var ctSL23 = FG.CharacterTableEmpty(sl23);
        ctSL23.DerivedSubGroupLift();
        ctSL23.RestrictionFromSuperGroup(ctGL23);
        ctSL23.DisplayCells(tableOnly: true);
        Console.WriteLine();
    }
    else
    {
        var ctG = FG.CharacterTableEmpty(g);
        ctG.DerivedSubGroupLift();
        ctG.InductionFromSubGroups(subgroups);
        ctG.DisplayCells(tableOnly: true);
        Console.WriteLine();
    }
}

void DisplayGroupDetails(IEnumerable<(AllSubgroups<WElt>subsg, ANameElt[] names)> seq, bool rename = false)
{
    var lt = seq.OrderBy(e => e.subsg.Parent.Count())
        .ThenBy(e => e.subsg.Parent.GroupType)
        .ThenBy(e => e.names[0])
        .ThenByDescending(e => e.subsg.Parent.ElementsOrders.Values.Max())
        .ThenBy(e => e.subsg.Infos)
        .ToArray();
    var dicOrd = lt.Select(e => e.subsg.Parent.Count()).Distinct().ToDictionary(k => k, _ => 0);
    var maxLt = rename ? lt.Max(e => e.names[0].Name.Length) : lt.Max(e => e.subsg.Parent.Name.Length);
    var nb = lt.Length;
    foreach (var (subsg, names) in lt)
    {
        var o = subsg.Parent.Count();
        Console.WriteLine($"Group{o}[{++dicOrd[o]}]");
        GroupDetails(subsg, names, rename, maxLt);
    }

    Console.WriteLine($"Total Groups:{nb}");
}

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

    var allOrderLess32 = allAb.AppendIsomorphic(allMtCycSdp, ext8, ext12, ext16, ext18, ext27, ext24, ext32)
        .Naming()
        .ToArray();

    DisplayGroupDetails(allOrderLess32);

    GlobalStopWatch.Show();
    Console.Beep();
    
    // Total Groups:144
    // #  Time:25.414s
    // 
}

{
    AllGroupsOrderLess32();
}
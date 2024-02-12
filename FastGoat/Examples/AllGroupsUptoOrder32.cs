using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.Naming;
using FastGoat.Structures.Subgroups;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Words.Tools;

namespace FastGoat.Examples;

public static class AllGroupsUptoOrder32
{
    static void GroupDetails(AllSubgroups<WElt> subgroups, ANameElt[] names, int maxLt = -1)
    {
        var nbSharp = 16;
        var (g, infos) = (subgroups.Parent, subgroups.Infos);
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

        var comChain = Group.CommutatorsChain(g);
        comChain.ForEach(sg =>
            sg.Name = sg.Count() == g.Count() ? g.Name : subgroups.First(cj => cj.Contains(sg)).Representative.Name);
        var comChainComplete = comChain.Last().Count() == 1;

        var zentrumsChain = Group.ZentrumsChainFast(g);
        zentrumsChain.ForEach(
            sg => sg.Name = sg.Count() == g.Count() ? g.Name : subgroups.First(cj => cj.Contains(sg)).Representative.Name);
        var isNilpotent = zentrumsChain.Last().Count() == g.Count();

        var derivedChain = Group.DerivedChain(g);
        derivedChain.ForEach(sg =>
            sg.Name = sg.Count() == g.Count() ? g.Name : subgroups.First(cj => cj.Contains(sg)).Representative.Name);
        var isSolvable = derivedChain.Last().Count() == 1;

        Console.WriteLine(g.ShortName);
        var simplicity = subgroups.IsSimple() ? "Simple" : "NotSimple";
        var nilpotency = isNilpotent ? "Nilpotent" : "NotNilpotent";
        var solubility = isSolvable ? "Solvable" : "NotSolvable";
        Console.WriteLine($"{simplicity}, {g.GroupType}, {nilpotency}, {solubility}");
        DisplayGroup.Orders(g);
        Console.WriteLine();

        Console.WriteLine(infos);
        Console.WriteLine($"Lower Central {(comChainComplete ? "Serie" : "Chain")}");
        Console.WriteLine(comChain.Glue(" ---> "));
        Console.WriteLine($"Upper Central {(isNilpotent ? "Serie" : "Chain")}");
        Console.WriteLine(zentrumsChain.Glue(" ---> "));
        Console.WriteLine($"Derived {(isSolvable ? "Serie" : "Chain")}");
        Console.WriteLine(derivedChain.Glue(" ---> "));

        var zentrum = Group.Zentrum(g);
        var zentrumName = zentrum.Count() == g.Count() ? g.Name : subgroups.First(cj => cj.Contains(zentrum)).Representative.Name;
        Console.WriteLine($"Zentrum  Z(G) = {zentrumName}");
        
        var frattini = subgroups.FrattiniSubGroup;
        var fratName = subgroups.First(cj => cj.Contains(frattini)).Representative.Name;
        Console.WriteLine($"Frattini Φ(G) = {fratName}");

        var fitting = subgroups.FittingSubGroup;
        var fitName = subgroups.First(cj => cj.Contains(fitting)).Representative.Name;
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

    static void DisplayGroupDetails(IEnumerable<(AllSubgroups<WElt>subsg, ANameElt[] names)> seq, bool html)
    {
        if (html)
        {
            Console.WriteLine(@"<!DOCTYPE html>");
            Console.WriteLine(@"<html>");
            Console.WriteLine(@"<head><style>div { page-break-after: always; font-family: monospace; white-space: pre; }</style></head>");
            Console.WriteLine(@"<body>");
        }
        
        var lt = seq.OrderBy(e => e.subsg.Parent.Count())
            .ThenBy(e => e.subsg.Parent.GroupType)
            .ThenBy(e => e.names[0])
            .ThenByDescending(e => e.subsg.Parent.ElementsOrders.Values.Max())
            .ThenBy(e => e.subsg.Infos)
            .ToArray();
        var dicOrd = lt.Select(e => e.subsg.Parent.Count()).Distinct().ToDictionary(k => k, _ => 0);
        var maxLt = lt.Max(e => e.names[0].Name.Length);
        var nb = lt.Length;
        foreach (var (subsg, names) in lt)
        {
            var o = subsg.Parent.Count();
            if (html)
                Console.WriteLine(@"<div>");
            
            Console.WriteLine($"Group{o}[{++dicOrd[o]}]");
            GroupDetails(subsg, names, maxLt);

            if (html)
                Console.WriteLine(@"</div>");
        }

        Console.WriteLine($"Total Groups:{nb}");
        
        if (html)
        {
            Console.WriteLine(@"</body>");
            Console.WriteLine(@"</html>");
        }
    }

    static void DetailsGroupsUptoOrder32(bool html = false)
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

        DisplayGroupDetails(allOrderLess32, html);

        GlobalStopWatch.Show();
        Console.Beep();
    }

    public static void Run()
    {
        Ring.DisplayPolynomial = MonomDisplay.StarCaret;
        DetailsGroupsUptoOrder32(html: true); // #  Time:15m50s
    }

    public static void A5()
    {
        Ring.DisplayPolynomial = MonomDisplay.StarCaret;
        var g = FG.Alternate(5).ToGroupWrapper();
        var subs = g.AllSubgroups();
        var names = NamesTree.BuildName(subs);
        GroupDetails(subs, names);
    }
    
    public static void Group_81_10()
    {
        Ring.DisplayPolynomial = MonomDisplay.StarCaret;
        var g = FG.WordGroup("a3b3, ababa-2b, ab2ab-1ab-1").ToGroupWrapper();
        var subs = g.AllSubgroups();
        var names = NamesTree.BuildName(subs);
        GroupDetails(subs, names);
    }
}
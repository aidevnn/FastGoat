using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Naming;

namespace FastGoat.Examples;

public static class GroupNaming
{
    static void ShowNames<T>(ConcreteGroup<T> G) where T : struct, IElt<T>
    {
        NamesTree.BuildName(G).Println($"Group:{G.ShortName}");
    }

    public static void Example1()
    {
        ShowNames(FG.Symmetric(3));
        ShowNames(FG.Dihedral(4));
        ShowNames(FG.DihedralSdp(5));
        ShowNames(FG.Dihedral(8));
        ShowNames(FG.Alternate(4));
        ShowNames(FG.Symmetric(4));
        ShowNames(Product.Generate(new Cn(5), Group.SemiDirectProd(new Cn(3), new Cn(4))));
        ShowNames(FG.Quaternion(8));
        ShowNames(FG.DiCyclic(6));
        ShowNames(FG.DiCyclic(7));
        ShowNames(FG.DiCyclic(8));
        ShowNames(FG.SemiDihedral(5));
    }

    public static void Example2()
    {
        var exts24 = FG.AllExtensions((FG.Abelian(12), FG.Abelian(2)))
            .OrderBy(e => e.ext.GroupType)
            .ThenByDescending(e => e.ext.ElementsOrders.Values.Max())
            .ThenBy(e => ((int, int, int))e.allSubs.Infos).ToList();

        foreach (var extInfos in exts24)
        {
            var it = NamesTree.BuildName(extInfos.allSubs.ToTable());
            extInfos.ext.Name = it.First().Name;
            CocyclesDFS.DisplayInfosGroups([(extInfos.ext, ((int, int, int))extInfos.allSubs.Infos)], naming: false);
            it.Println("Group Names");
        }
    }

    public static void Example3()
    {
        ShowNames(FG.Abelian(5));
        ShowNames(FG.Alternate(5));
        ShowNames(FG.Symmetric(5));
        ShowNames(FG.Abelian(7, 3, 2));
        ShowNames(FG.Abelian(5, 2));
        ShowNames(FG.GLnp(3, 2));
        ShowNames(FG.SL2p(3));
        ShowNames(FG.SL2p(5));

        // Top Down simple non abelian group construction, name dont changes
        ShowNames(FG.L2p(5));
        ShowNames(FG.L2p(7));
        ShowNames(FG.L2p(11));
    }

    public static void Example4()
    {
        var allExts20 = FG.AllExtensions(FG.AllAbelianGroupsOfOrder(4).Select(e => (FG.Abelian(5), e)).ToArray())
            .OrderBy(e => e.ext.GroupType)
            .ThenByDescending(e => e.ext.ElementsOrders.Values.Max())
            .ThenBy(e => ((int, int, int))e.allSubs.Infos).ToList();

        var allExts40 = FG.AllExtensions(FG.AllAbelianGroupsOfOrder(4).Select(e => (FG.Abelian(2, 5), e)).ToArray())
            .OrderBy(e => e.ext.GroupType)
            .ThenByDescending(e => e.ext.ElementsOrders.Values.Max())
            .ThenBy(e => ((int, int, int))e.allSubs.Infos).ToList();

        foreach (var extInfos in allExts20)
        {
            var it = NamesTree.BuildName(extInfos.allSubs.ToTable());
            extInfos.ext.Name = it.First().Name;
            CocyclesDFS.DisplayInfosGroups([(extInfos.ext, ((int, int, int))extInfos.allSubs.Infos)], naming: false);
            it.Println("Group Names");
        }

        Console.WriteLine();

        foreach (var extInfos in allExts40)
        {
            var it = NamesTree.BuildName(extInfos.allSubs.ToTable());
            extInfos.ext.Name = it.First().Name;
            CocyclesDFS.DisplayInfosGroups([(extInfos.ext, ((int, int, int))extInfos.allSubs.Infos)], naming: false);
            it.Println("Group Names");
        }
    }

    public static void Example5()
    {
        var allExts = FG.AllExtensions((FG.Abelian(4, 4), FG.Abelian(2)))
            .OrderBy(e => e.ext.GroupType)
            .ThenByDescending(e => e.ext.ElementsOrders.Values.Max())
            .ThenBy(e => ((int, int, int))e.allSubs.Infos).ToList();

        foreach (var extInfos in allExts)
        {
            var it = NamesTree.BuildName(extInfos.allSubs.ToTable());
            extInfos.ext.Name = it.First().Name;
            CocyclesDFS.DisplayInfosGroups([(extInfos.ext, ((int, int, int))extInfos.allSubs.Infos)], naming: false);
            it.Println("Group Names");
        }
    }

    public static void Example5_1()
    {
        var allExts = FG.AllExtensions((FG.ModularMaxSdp(4), FG.Abelian(2)))
            .OrderBy(e => e.ext.GroupType)
            .ThenByDescending(e => e.ext.ElementsOrders.Values.Max())
            .ThenBy(e => ((int, int, int))e.allSubs.Infos).ToList();

        foreach (var extInfos in allExts)
        {
            var it = NamesTree.BuildName(extInfos.allSubs.ToTable());
            extInfos.ext.Name = it.First().Name;
            CocyclesDFS.DisplayInfosGroups([(extInfos.ext, ((int, int, int))extInfos.allSubs.Infos)], naming: false);
            it.Println("Group Names");
        }
    }

    public static void Example6()
    {
        GlobalStopWatch.Restart();
        var nbOpsMax = 10000;
        var allExts = FG.AllExtensions(
                (nbOpsMax, FG.Abelian(16), FG.Abelian(2)),
                (nbOpsMax, FG.Abelian(2, 8), FG.Abelian(2)),
                (nbOpsMax, FG.Abelian(4, 4), FG.Abelian(2)),
                (nbOpsMax, FG.Abelian(2, 2, 4), FG.Abelian(2)),
                (nbOpsMax, FG.Abelian(2, 4), FG.Abelian(4)),
                (16, FG.Abelian(2, 4), FG.Abelian(2, 2)),
                (1, FG.Abelian(2, 2, 2), FG.Abelian(2, 2)))
            .Take(51)
            .OrderBy(e => e.ext.GroupType)
            .ThenByDescending(e => e.ext.ElementsOrders.Values.Max())
            .ThenBy(e => ((int, int, int))e.allSubs.Infos).ToList();

        foreach (var extInfos in allExts)
        {
            var it = NamesTree.BuildName(extInfos.allSubs.ToTable());
            extInfos.ext.Name = it.First().Name;
            CocyclesDFS.DisplayInfosGroups([(extInfos.ext, ((int, int, int))extInfos.allSubs.Infos)], naming: false);
            it.Println("Group Names");
        }

        GlobalStopWatch.Show($"Ord32: {allExts.Count}"); // Time:199294 ms ~ 3min
        Console.Beep();
        Console.Write("Checking that all extensions are valid groups...");
        if (allExts.Any(e => !Group.IsGroup(e.ext)))
            throw new GroupException(GroupExceptionType.GroupDef);
        Console.WriteLine(" Done.");
    }

    public static void Example6_1()
    {
        GlobalStopWatch.Restart();
        var allExts = FG.AllExtensions(FG.AllAbelianGroupsOfOrder(27).Select(e => (e, FG.Abelian(3))).ToArray())
            .Take(15)
            .OrderBy(e => e.ext.GroupType)
            .ThenByDescending(e => e.ext.ElementsOrders.Values.Max())
            .ThenBy(e => ((int, int, int))e.allSubs.Infos).ToList();

        foreach (var extInfos in allExts)
        {
            var it = NamesTree.BuildName(extInfos.allSubs.ToTable());
            extInfos.ext.Name = it.First().Name;
            CocyclesDFS.DisplayInfosGroups([(extInfos.ext, ((int, int, int))extInfos.allSubs.Infos)], naming: false);
            it.Println("Group Names");
        }

        GlobalStopWatch.Show($"Ord81: {allExts.Count}"); // Time:130697 ms ~ 2min
        Console.Beep();
        Console.Write("Checking that all extensions are valid groups...");
        if (allExts.Any(e => !Group.IsGroup(e.ext)))
            throw new GroupException(GroupExceptionType.GroupDef);
        Console.WriteLine(" Done.");
    }

    public static void Example6_2()
    {
        GlobalStopWatch.Restart();
        var allExts = FG.AllExtensions((FG.Abelian(2), FG.Abelian(4, 4)))
            .Take(15)
            .OrderBy(e => e.ext.GroupType)
            .ThenByDescending(e => e.ext.ElementsOrders.Values.Max())
            .ThenBy(e => ((int, int, int))e.allSubs.Infos).ToList();

        foreach (var extInfos in allExts)
        {
            var it = NamesTree.BuildName(extInfos.allSubs.ToTable());
            extInfos.ext.Name = it.First().Name;
            CocyclesDFS.DisplayInfosGroups([(extInfos.ext, ((int, int, int))extInfos.allSubs.Infos)], naming: false);
            it.Println("Group Names");
        }

        GlobalStopWatch.Show($"C2 . (C4 x C4): {allExts.Count}"); // Time:595836 ms ~ 10min
        Console.Beep();
        Console.Write("Checking that all extensions are valid groups...");
        if (allExts.Any(e => !Group.IsGroup(e.ext)))
            throw new GroupException(GroupExceptionType.GroupDef);
        Console.WriteLine(" Done.");
    }

    public static void Example7()
    {
        GlobalStopWatch.Restart();
        var allExts = FG.AllExtensions((FG.SL2p(3), FG.Abelian(2)))
            .OrderBy(e => e.ext.GroupType)
            .ThenByDescending(e => e.ext.ElementsOrders.Values.Max())
            .ThenBy(e => ((int, int, int))e.allSubs.Infos).ToList();

        foreach (var extInfos in allExts)
        {
            var it = NamesTree.BuildName(extInfos.allSubs.ToTable());
            extInfos.ext.Name = it.First().Name;
            CocyclesDFS.DisplayInfosGroups([(extInfos.ext, ((int, int, int))extInfos.allSubs.Infos)], naming: false);
            it.Println("Group Names");
        }

        GlobalStopWatch.Show($"SL(2,3) . C2: {allExts.Count}"); // Time:5497 ms
        Console.Beep();
        Console.Write("Checking that all extensions are valid groups...");
        if (allExts.Any(e => !Group.IsGroup(e.ext)))
            throw new GroupException(GroupExceptionType.GroupDef);
        Console.WriteLine(" Done.");
    }
    /* 
       ###################################################
       #################  SL(2,3) x: C2  #################
       ###################################################
       |SL(2,3) x: C2| = 48
       Type        NonAbelianGroup
       BaseGroup   SL(2,3) . C2
       
       Elements Orders : [1]:1, [2]:7, [3]:8, [4]:8, [6]:8, [12]:16
       
       AllSubGr:37 AllConjsCl:15 AllNorms:7
       
       Group Names
           SL(2,3) x: C2
           (D8 x: C2) x: C3
           C4 . A4
           Q8 . C6
           C2 . (C2 x A4)
       #############################################
       #################  C2 . S4  #################
       #############################################
       |C2 . S4| = 48
       Type        NonAbelianGroup
       BaseGroup   SL(2,3) . C2
       
       Elements Orders : [1]:1, [2]:1, [3]:8, [4]:18, [6]:8, [8]:12
       
       AllSubGr:35 AllConjsCl:13 AllNorms:5
       
       Group Names
           C2 . S4
           Q8 . S3
           SL(2,3) . C2
       ##############################################
       #################  Q8 x: S3  #################
       ##############################################
       |Q8 x: S3| = 48
       Type        NonAbelianGroup
       BaseGroup   SL(2,3) . C2
       
       Elements Orders : [1]:1, [2]:13, [3]:8, [4]:6, [6]:8, [8]:12
       
       AllSubGr:55 AllConjsCl:16 AllNorms:5
       
       Group Names
           Q8 x: S3
           SL(2,3) x: C2
           C2 . S4
       ##################################################
       #################  C2 x SL(2,3)  #################
       ##################################################
       |C2 x SL(2,3)| = 48
       Type        NonAbelianGroup
       BaseGroup   SL(2,3) . C2
       
       Elements Orders : [1]:1, [2]:3, [3]:8, [4]:12, [6]:24
       
       AllSubGr:41 AllConjsCl:18 AllNorms:9
       
       Group Names
           C2 x SL(2,3)
           Q8 x: C6
           (C2 x Q8) x: C3
           (C2 x C2) . A4
           C2 . (C2 x A4)
       # SL(2,3) . C2: 4 Time:5326 ms
       Checking that all extensions are valid groups... Done.
     */
}
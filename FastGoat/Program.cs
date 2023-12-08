using System.Collections;
using System.ComponentModel;
using System.Diagnostics;
using System.Numerics;
using FastGoat.Commons;
using FastGoat.Examples;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
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

IEnumerable<GroupTable> AllProducts(ConcreteGroup<TableElt> a, ConcreteGroup<TableElt> b)
{
    yield return Product.Generate(a, b).ToTable().gt;
    foreach (var e in Group.AllSemiDirectProd(a, b))
        yield return e.ToTable().gt;
    foreach (var e in Group.AllSemiDirectProd(b, a))
        yield return e.ToTable().gt;
}

(ConcreteGroup<TableElt> g, AllSubgroups<TableElt> subGroups, GroupNaming.ITreeElt<TableElt>[] names)[]
    AllGroupNames(IEnumerable<ConcreteGroup<TableElt>> elts)
{
    var dicExts = new Dictionary<SubGroupsInfos, HashSet<ConcreteGroup<TableElt>>>();
    var list = new List<(ConcreteGroup<TableElt> g, AllSubgroups<TableElt> subGroups, GroupNaming.ITreeElt<TableElt>[] names)>();
    foreach (var g in elts)
    {
        var subGroups = new AllSubgroups<TableElt>(g);
        var it = GroupNaming.BuildName(subGroups);
        g.Name = it[0].Name;
        var gn = (g, subGroups, it);
        var infos = gn.subGroups.Infos;
        if (dicExts.ContainsKey(infos))
        {
            if (dicExts[infos].Add(g))
                list.Add(gn);
        }
        else
        {
            dicExts[infos] = new(new IsomorphEquality<TableElt>()) { g };
            list.Add(gn);
        }
    }

    return list.OrderBy(e => e.g.GroupType)
        .ThenByDescending(e => e.g.ElementsOrders.Values.Max())
        .ThenBy(e => ((int, int, int))e.subGroups.Infos)
        .ToArray();
}


void DisplayGroupNames((ConcreteGroup<TableElt> g, AllSubgroups<TableElt> subGroups, GroupNaming.ITreeElt<TableElt>[] names)[] elts)
{
    var maxLt = elts.Select(e => e.g.Name).Max(e => e.Length);
    var lt = Enumerable.Repeat('#', maxLt + 4).Glue();
    var line = $"#################{lt}#################";
    var fmt = $"#################  {{0,{-maxLt}}}  #################";

    foreach (var (g, subgroups, it) in elts)
    {
        Console.WriteLine(line);
        Console.WriteLine(fmt, g.Name);
        Console.WriteLine(line);
        DisplayGroup.HeadOrders(g);
        Console.WriteLine(subgroups.Infos);
        it.Println("Group Names");
        Console.WriteLine();
    }

    Console.WriteLine($"Count:{elts.Length}");
    Console.WriteLine();
}

void ProductOrd48()
{
    var (c2, c3, c4) = 3.Range(2).Select(i => new Cn(i).ToTable().gt).Deconstruct();

    var allOrd8 = AllGroupNames([
            ..FG.AllAbelianGroupsOfOrder(8).Select(e => e.ToTable().gt),
            FG.Dihedral(4).ToTable().gt,
            FG.Quaternion(8).ToTable().gt
        ]
    );
    DisplayGroupNames(allOrd8);

    var allOrd12 = AllGroupNames([
            ..FG.AllAbelianGroupsOfOrder(12).Select(e => e.ToTable().gt),
            FG.Dihedral(6).ToTable().gt,
            FG.Alternate(4).ToTable().gt,
            FG.Frobenius(12)[0].ToTable().gt
        ]
    );
    DisplayGroupNames(allOrd12);

    var ord16a = FG.AllAbelianGroupsOfOrder(4).ToArray().Grid2D().SelectMany(e => AllProducts(e.t1.ToTable().gt, e.t2.ToTable().gt));
    var ord16b = allOrd8.SelectMany(e => AllProducts(e.g, c2)).Append(FG.Quaternion(16).ToTable().gt);
    var allOrd16 = AllGroupNames([..ord16a, ..ord16b, new Cn(16).ToTable().gt]);
    DisplayGroupNames(allOrd16);

    var ord24a = allOrd12.SelectMany(e => AllProducts(e.g, c2));
    var ord24b = allOrd8.SelectMany(e => AllProducts(e.g, c3));
    var allOrd24 = AllGroupNames([..ord24a, ..ord24b]);
    DisplayGroupNames(allOrd24);

    GlobalStopWatch.Restart();
    var ord48a = allOrd24.SelectMany(e => AllProducts(e.g, c2));
    var ord48b = allOrd16.SelectMany(e => AllProducts(e.g, c3));
    var ord48c = allOrd12.SelectMany(e => AllProducts(e.g, c4));
    var allOrd48 = AllGroupNames([..ord48a, ..ord48b, ..ord48c]);
    DisplayGroupNames(allOrd48);
    GlobalStopWatch.Show("Ord48");
    Console.Beep();
}

void ProductOrd32()
{
    var c2 = FG.Abelian(2).ToTable().gt;
    var allOrd4 = FG.AllAbelianGroupsOfOrder(4).Select(e => e.ToTable().gt).ToArray();

    var allOrd8 = AllGroupNames([
            ..FG.AllAbelianGroupsOfOrder(8).Select(e => e.ToTable().gt),
            FG.Dihedral(4).ToTable().gt,
            FG.Quaternion(8).ToTable().gt
        ]
    );

    var ord16a = allOrd4.Grid2D().SelectMany(e => AllProducts(e.t1, e.t2));
    var ord16b = allOrd8.SelectMany(e => AllProducts(e.g, c2)).Append(FG.Quaternion(16).ToTable().gt);
    var allOrd16 = AllGroupNames([..ord16a, ..ord16b, new Cn(16).ToTable().gt]);
    DisplayGroupNames(allOrd16);

    GlobalStopWatch.Restart();
    var ord32a = allOrd16.SelectMany(e => AllProducts(e.g, c2)).Append(FG.Quaternion(32).ToTable().gt);
    var ord32b = allOrd8.Grid2D(allOrd4).SelectMany(e => AllProducts(e.t1.g, e.t2));
    var allOrd32 = AllGroupNames([..ord32a, ..ord32b, new Cn(32).ToTable().gt]);
    DisplayGroupNames(allOrd32);
    GlobalStopWatch.Show("Ord32");
    Console.Beep();
}

{
    // ProductOrd32(); // missing 3 extensions
    // ProductOrd48(); // missing 1 extension
}

{

    var exts = FG.AllExtensions((Product.Generate(FG.Abelian(2, 2), FG.Symmetric(3)), FG.Abelian(2)))
        .OrderBy(e => e.ext.GroupType)
        .ThenByDescending(e => e.ext.ElementsOrders.Values.Max())
        .ThenBy(e => ((int, int, int))e.allSubs.Infos).ToList();
    
    foreach (var extInfos in exts)
    {
        var it = GroupNaming.BuildName(extInfos.allSubs);
        extInfos.ext.Name = it.First().Name;
        CocyclesDFS.DisplayInfosGroups([(extInfos.ext, ((int, int, int))extInfos.allSubs.Infos)], naming: false);
        it.Println("Group Names");
    }

    // ~ 10minutes
}
/* https://people.maths.bris.ac.uk/~matyd/GroupNames/1/e15/C2byC2%5E2xS3.html#s2
   6 Extensions : S3 x C2^3, C2 x D24, S3 x D8, C2 x (C3 x: D8), D12 x: C4, S3 x C2 x C4 

   Total Exts:6
   ###############################################
   #################  D12 x: C4  #################
   ###############################################
   |D12 x: C4| = 48
   Type        NonAbelianGroup
   BaseGroup   ((C2 x C2) x Symm3) . C2
   
   Elements Orders : [1]:1, [2]:15, [3]:2, [4]:16, [6]:6, [12]:8
   
   AllSubGr:76 AllConjsCl:34 AllNorms:17
   
   Group Names
       D12 x: C4
       (C4 x C2) x: S3
       (C12 x C2) x: C2
       M(6x:4)5 x: C2
       C3 x: ((C2 x C2) x: C4)
   ##############################################
   #################  C4 x D12  #################
   ##############################################
   |C4 x D12| = 48
   Type        NonAbelianGroup
   BaseGroup   ((C2 x C2) x Symm3) . C2
   
   Elements Orders : [1]:1, [2]:15, [3]:2, [4]:16, [6]:6, [12]:8
   
   AllSubGr:92 AllConjsCl:54 AllNorms:35
   
   Group Names
       C4 x D12
       (C4 x C2) x S3
       D12 x: C4
       (C4 x S3) x: C2
       C6 x: (C4 x C2)
       S3 x: (C4 x C2)
       (C12 x C2) x: C2
       C12 x: (C2 x C2)
       C3 x: (C4 x C2 x C2)
       M(6x:4)5 x: C2
       F(3x:4)2 x: (C2 x C2)
   #############################################
   #################  D8 x S3  #################
   #############################################
   |D8 x S3| = 48
   Type        NonAbelianGroup
   BaseGroup   ((C2 x C2) x Symm3) . C2
   
   Elements Orders : [1]:1, [2]:23, [3]:2, [4]:8, [6]:10, [12]:4
   
   AllSubGr:120 AllConjsCl:54 AllNorms:25
   
   Group Names
       D8 x S3
       D8 x: S3
       S3 x: D8
       C4 x: D12
       D24 x: C2
       (C3 x D8) x: C2
       (C4 x S3) x: C2
       C3 x: (C2 x D8)
       (C2 x C2) x: D12
       (C2 x D12) x: C2
       C12 x: (C2 x C2)
       D12 x: (C2 x C2)
       (C6 x C2) x: (C2 x C2)
       (C3 x: D8) x: C2
       F(3x:4)2 x: (C2 x C2)
   ##############################################
   #################  C2 x D24  #################
   ##############################################
   |C2 x D24| = 48
   Type        NonAbelianGroup
   BaseGroup   ((C2 x C2) x Symm3) . C2
   
   Elements Orders : [1]:1, [2]:27, [3]:2, [4]:4, [6]:6, [12]:8
   
   AllSubGr:124 AllConjsCl:54 AllNorms:27
   
   Group Names
       C2 x D24
       C6 x: D8
       C4 x: D12
       D24 x: C2
       (C4 x C2) x: S3
       C3 x: (C2 x D8)
       (C12 x C2) x: C2
       (C2 x D12) x: C2
       C12 x: (C2 x C2)
       D12 x: (C2 x C2)
   ##############################################
   #################  C6 x: D8  #################
   ##############################################
   |C6 x: D8| = 48
   Type        NonAbelianGroup
   BaseGroup   ((C2 x C2) x Symm3) . C2
   
   Elements Orders : [1]:1, [2]:19, [3]:2, [4]:12, [6]:14
   
   AllSubGr:108 AllConjsCl:54 AllNorms:27
   
   Group Names
       C6 x: D8
       C2 x (C3 x: D8)
       C3 x: (C2 x D8)
       (C2 x C2) x: D12
       (C2 x D12) x: C2
       D12 x: (C2 x C2)
       (C2 x C2 x C2) x: S3
       (C6 x C2 x C2) x: C2
       (C6 x C2) x: (C2 x C2)
       M(6x:4)5 x: C2
       (C3 x: D8) x: C2
       F(3x:4)2 x: (C2 x C2)
   #####################################################
   #################  (C2 x C2) x D12  #################
   #####################################################
   |(C2 x C2) x D12| = 48
   Type        NonAbelianGroup
   BaseGroup   ((C2 x C2) x Symm3) . C2
   
   Elements Orders : [1]:1, [2]:31, [3]:2, [6]:14
   
   AllSubGr:236 AllConjsCl:134 AllNorms:83
   
   Group Names
       (C2 x C2) x D12
       (C2 x C2 x C2) x S3
       (C2 x D12) x: C2
       D12 x: (C2 x C2)
       (C6 x C2 x C2) x: C2
       C6 x: (C2 x C2 x C2)
       S3 x: (C2 x C2 x C2)
       (C6 x C2) x: (C2 x C2)
       C3 x: (C2 x C2 x C2 x C2)
   
   Process finished with exit code 0.
   
*/
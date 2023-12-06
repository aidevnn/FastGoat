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

IEnumerable<ConcreteGroup<TableElt>> AllProducts(ConcreteGroup<TableElt> a, ConcreteGroup<TableElt> b)
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

    var allOrd8 = AllGroupNames(
        FG.AllAbelianGroupsOfOrder(8).Select(e => e.ToTable().gt)
            .Append(FG.Dihedral(4).ToTable().gt)
            .Append(FG.Quaternion(8).ToTable().gt)
    );
    DisplayGroupNames(allOrd8);

    var allOrd12 = AllGroupNames(
        FG.AllAbelianGroupsOfOrder(12).Select(e => e.ToTable().gt)
            .Append(FG.Dihedral(6).ToTable().gt)
            .Append(FG.Alternate(4).ToTable().gt)
            .Append(FG.Frobenius(12)[0].ToTable().gt)
    );
    DisplayGroupNames(allOrd12);

    var ord16a = FG.AllAbelianGroupsOfOrder(4).ToArray().Grid2D().SelectMany(e => AllProducts(e.t1.ToTable().gt, e.t2.ToTable().gt));
    var ord16b = allOrd8.SelectMany(e => AllProducts(e.g, c2)).Append(FG.Quaternion(16).ToTable().gt);
    var allOrd16 = AllGroupNames(ord16a.Concat(ord16b).Append(new Cn(16).ToTable().gt));
    DisplayGroupNames(allOrd16);
    
    var ord24a = allOrd12.SelectMany(e => AllProducts(e.g, c2));
    var ord24b = allOrd8.SelectMany(e => AllProducts(e.g, c3));
    var allOrd24 = AllGroupNames(ord24a.Concat(ord24b));
    DisplayGroupNames(allOrd24);

    GlobalStopWatch.Restart();
    var ord48a = allOrd24.SelectMany(e => AllProducts(e.g, c2));
    var ord48b = allOrd16.SelectMany(e => AllProducts(e.g, c3));
    var ord48c = allOrd12.SelectMany(e => AllProducts(e.g, c4));
    var allOrd48 = AllGroupNames( ord48a.Concat(ord48b).Concat(ord48c));
    DisplayGroupNames(allOrd48);
    GlobalStopWatch.Show("Ord48");
    Console.Beep();
}

void ProductOrd32()
{
    var c2 = FG.Abelian(2).ToTable().gt;
    var allOrd4 = FG.AllAbelianGroupsOfOrder(4).Select(e => e.ToTable().gt).ToArray();

    var allOrd8 = AllGroupNames(
        FG.AllAbelianGroupsOfOrder(8).Select(e => e.ToTable().gt)
            .Append(FG.Dihedral(4).ToTable().gt)
            .Append(FG.Quaternion(8).ToTable().gt)
    );

    var ord16a = allOrd4.Grid2D().SelectMany(e => AllProducts(e.t1, e.t2));
    var ord16b = allOrd8.SelectMany(e => AllProducts(e.g, c2)).Append(FG.Quaternion(16).ToTable().gt);
    var allOrd16 = AllGroupNames(ord16a.Concat(ord16b).Append(new Cn(16).ToTable().gt));
    DisplayGroupNames(allOrd16);
    
    GlobalStopWatch.Restart();
    var ord32a = allOrd16.SelectMany(e => AllProducts(e.g, c2)).Append(FG.Quaternion(32).ToTable().gt);
    var ord32b = allOrd8.Grid2D(allOrd4).SelectMany(e => AllProducts(e.t1.g, e.t2));
    var allOrd32 = AllGroupNames(ord32a.Concat(ord32b).Append(new Cn(32).ToTable().gt));
    DisplayGroupNames(allOrd32);
    GlobalStopWatch.Show("Ord32");
    Console.Beep();
}

{
    // ProductOrd32(); // missing 3 extensions
    ProductOrd48(); // missing 1 extension
}

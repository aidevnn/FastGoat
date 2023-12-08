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

{
    var exts = FG.AllExtensions((FG.Abelian(4), FG.Abelian(2, 2)))
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

    // Found 9 extensions expected 8 extensions according to
    // https://people.maths.bris.ac.uk/~matyd/GroupNames/1/e3/C2%5E2byC4.html#d6
    // C8 x C2 is extra ???
}

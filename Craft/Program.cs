using System.Numerics;
using System.Reflection;
using System.Text;
using Craft;
using Craft.Craft;
using Examples;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.Naming;
using FastGoat.Structures.Subgroups;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.DatabaseSmallGroups;
using FastGoat.UserGroup.EllCurve;
using FastGoat.UserGroup.Floats;
using FastGoat.UserGroup.GModuleN;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Matrix;
using FastGoat.UserGroup.Perms;
using FastGoat.UserGroup.Polynoms;
using FastGoat.UserGroup.Words;
using FastGoat.UserGroup.Words.Tools;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");

Perm.Style = DisplayPerm.CyclesComma;

string GapExport(Perm[] gens)
{
    var old = Perm.Style;
    Perm.Style = DisplayPerm.Gap;
    var export = $"Group([{gens.Glue(", ")}]);";
    Perm.Style = old;
    return export;
}

{
    for (int m = 3; m <= 64; m++)
    {
        var (dicm, autDicm) = FG.AutomorphismDicyclicPg(m);
        DisplayGroup.HeadOrders(dicm);
        DisplayGroup.HeadOrdersGenerators(autDicm); // for m>=3, Aut[Dic[m]] ~ Aut[D4m] ~ Hol[C2m]
        if (m >= 3)
        {
            var holC2m = UGCraft.HolCn(FG.Cycles(2 * m)).ToHashSet();
            if (!autDicm.SetEquals(holC2m))
                throw new();
        }

        Console.WriteLine($"autDic{m}:={GapExport(autDicm.GetGenerators().ToArray())}");
        if (m < 32)
            Console.WriteLine($"IdGroup(autDic{m});IdGroup(AutomorphismGroup(DicyclicGroup({4 * m})));");
        Console.WriteLine();
    }
}

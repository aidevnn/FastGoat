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

{
    GlobalStopWatch.Restart();

    for (int ord = 1; ord <= 32; ++ord)
    {
        foreach (var G in FG.AllGroupsOfOrder(ord).Where(g => g.GroupType == GroupType.NonAbelianGroup))
        {
            if (!Group.IsSoluble(G))
                continue;
            
            var pols = PolyCyclicSequence<Word>.Decomp(G).ToArray();
            Console.WriteLine($"##################### {G.ShortName}");
            foreach (var polyCyclic in pols)
            {
                polyCyclic.Generators.Select(e => (e.e, e.K.Order, e.Split))
                    .Println($"{polyCyclic.Normals.Select((N, i) => $"G{i}[{N.Order}]").Glue(" ≤ ")} = G");
                
                var rels = polyCyclic.Relators();
                var G0 = FG.WordGroup("G", rels);
                if (!GroupCraft.AreIsomorphic(G, G0))
                    throw new();
            }
        }
    }

    GlobalStopWatch.Show();
    
    // ##################### |Dic6| = 24
    // G0[1] ≤ G1[12] ≤ G2[24] = G
    //     (a-1, 12, True)
    //     (b-1, 2, False)
    //     ((), 1, True)
    // 
    // PCGroup Relators: a2 = b6, b12, a-1ba = b11
    // 
    // G0[1] ≤ G1[3] ≤ G2[12] ≤ G3[24] = G
    //     (a4, 3, True)
    //     (b-1, 4, True)
    //     (a-1, 2, False)
    //     ((), 1, True)
    // 
    // PCGroup Relators: a2 = b2c, b4, c3, a-1ba = b3c, a-1ca = c, b-1cb = c2
    // 
    // G0[1] ≤ G1[3] ≤ G2[12] ≤ G3[24] = G
    //     (a4, 3, True)
    //     (a-1b-1, 4, True)
    //     (a-1, 2, False)
    //     ((), 1, True)
    // 
    // PCGroup Relators: a2 = b2c, b4, c3, a-1ba = b3c, a-1ca = c, b-1cb = c2
    // 
}

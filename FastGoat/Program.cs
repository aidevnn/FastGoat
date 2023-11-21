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

IEnumerable<(ConcreteGroup<Ep2<Tn, Tg>> ext, Dictionary<ConcreteGroup<Ep2<Tn, Tg>>, List<ConcreteGroup<Ep2<Tn, Tg>>>> allSubs, 
    (int, int, int) infos)> AllExtensions<Tn, Tg>(ConcreteGroup<Tn> N, ConcreteGroup<Tg> G)
    where Tg : struct, IElt<Tg>
    where Tn : struct, IElt<Tn>
{
    var CN = Group.Zentrum(N);
    var autN = Group.AutomorphismGroup(N);
    var ops = Group.AllHomomorphisms(G, autN);
    var dicExts = new Dictionary<(int, int, int), HashSet<ConcreteGroup<Ep2<Tn, Tg>>>>();
    foreach (var (op, i) in ops.Select((l, i) => (l, i + 1)))
    {
        var L = op.ToMapElt(autN);
        var lbl = $"Lbl{i}/{ops.Count}";
        var (cohs, cobs, cocs) = ZNSolver.ReduceCohomologies(CN, G, L, lbl: lbl);
        foreach (var c in cohs)
        {
            var c0 = c.ToMapElt;
            var ext = Group.ExtensionGroup(N, L, c0, G);
            var extBase = ((ExtensionGroupBase<Tn, Tg>)ext.BaseGroup);
            if (extBase.IsGroup)
            {
                var allSubs = Group.AllSubGroups(ext);
                var infos = CocyclesDFS.SubGroupsDetails(allSubs);
                if (dicExts.ContainsKey(infos))
                {
                    if (dicExts[infos].Add(ext))
                        yield return (ext, allSubs, infos);
                }
                else
                {
                    dicExts[infos] = new(new IsomorphEquality<Ep2<Tn, Tg>>()) { ext };
                    yield return (ext, allSubs, infos);
                }
            }
            else
            {
                Console.WriteLine("????????????????????? Extension isnt a group");
            }
        }

        Console.WriteLine($"Nb Exts:{dicExts.Values.Sum(v => v.Count)}");
    }
}

{
    GlobalStopWatch.Restart();
    var allExts = AllExtensions(N: FG.Abelian(2), G: FG.Abelian(4, 4)).ToArray();
    CocyclesDFS.DisplayInfosGroups(allExts.Select(e => (e.ext, e.infos)).ToArray(), naming: true, prefix: "Ext");
    GlobalStopWatch.Show();
    Console.Beep();
}

// 
// ############# Lbl1/1       #############
// H2(G, N) with N:|Z(C2)| = 2 and G:|C4 x C4| = 16
// #### |B2|:8192 |Z2|:65536 |H2|:8
// 
// Step:12 Gens:3/16 Dim:65536/65536
// B2(G,N):8192=2x2x2x2x2x2x2x2x2x2x2x2x2
// Z2(G,N):65536=2x2x2x2x2x2x2x2x2x2x2x2x2x2x2x2
// Cosets:8/8
// H2(G,N):8=1x2x2x2 Expected:8
// Nb Exts:4
// ##########################################################
// #################    Ext32[1] found   ####################
// ##########################################################
// |Ext32[1]| = 32
// Type        AbelianGroup
// BaseGroup   C2 . (C4 x C4)
// 
// Elements Orders : [1]:1, [2]:3, [4]:12, [8]:16
// 
// AllSubGr:22 AllConjsCl:22 AllNorms:22
// 
// ##########################################################
// #################    Ext32[2] found   ####################
// ##########################################################
// |Ext32[2]| = 32
// Type        AbelianGroup
// BaseGroup   C2 . (C4 x C4)
// 
// Elements Orders : [1]:1, [2]:7, [4]:24
// 
// AllSubGr:54 AllConjsCl:54 AllNorms:54
// 
// ##########################################################
// #################    Ext32[3] found   ####################
// ##########################################################
// |Ext32[3]| = 32
// Type        NonAbelianGroup
// BaseGroup   C2 . (C4 x C4)
// 
// Elements Orders : [1]:1, [2]:3, [4]:12, [8]:16
// 
// AllSubGr:22 AllConjsCl:20 AllNorms:18
// 
// ##########################################################
// #################    Ext32[4] found   ####################
// ##########################################################
// |Ext32[4]| = 32
// Type        NonAbelianGroup
// BaseGroup   C2 . (C4 x C4)
// 
// Elements Orders : [1]:1, [2]:7, [4]:24
// 
// AllSubGr:50 AllConjsCl:38 AllNorms:26
// 
// #  Time:779182 ms
// 
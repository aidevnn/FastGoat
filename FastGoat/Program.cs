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
    var g0 = FG.Abelian(2, 4);
    var allExts = AllExtensions(N: g0, G: g0);
    var nonSplits = allExts.Where(extInfos =>
    {
        if (extInfos.ext.GroupType == GroupType.AbelianGroup)
            return false;

        // Normal Subgroup C2 x C4
        var lt0 = extInfos.allSubs.Where(e => 
                e.Value.Count == 1 && 
                e.Key.GroupType == GroupType.AbelianGroup && 
                e.Key.IsIsomorphicTo(g0))
            .ToArray();
        if (lt0.Length == 0)
            return false;

        // Ext/(C2 x C4) is isomorphic to C2 x C4
        var lt1 = lt0.Select(e => extInfos.ext.Over(e.Key)).Where(e => e.GroupType == GroupType.AbelianGroup && e.IsIsomorphicTo(g0))
            .ToArray();
        if (lt1.Length == 0)
            return false;

        return NonSplitExtension.AllSplittingGroups(g0, extInfos.ext, g0, details: false).Any(e => e.s.IsNull);
    }).Take(10).ToArray();
    CocyclesDFS.DisplayInfosGroups(nonSplits.Select(e => (e.ext, e.infos)).ToArray(), naming: true, prefix: "Ext");
    GlobalStopWatch.Show();
    Console.Beep();
}

//
// ############# Lbl1/32      #############
// H2(G, N) with N:|Z(C2 x C4)| = 8 and G:|C2 x C4| = 8
// #### |B2|:65536 |Z2|:8388608 |H2|:128
// 
// Step:9 Gens:6/16 Dim:8388608/8388608
// B2(G,N):65536=4x4x4x4x4x2x2x2x2x2x2
// Z2(G,N):8388608=4x4x4x4x4x4x4x2x2x2x2x2x2x2x2x2
// Cosets:128/128
// H2(G,N):128=1x4x2x2x2x2x2 Expected:128
// ##########################################################
// #################    Ext64[1] found   ####################
// ##########################################################
// |Ext64[1]| = 64
// Type        NonAbelianGroup
// BaseGroup   (C2 x C4) . (C2 x C4)
// 
// Elements Orders : [1]:1, [2]:15, [4]:48
// 
// AllSubGr:185 AllConjsCl:129 AllNorms:73
// 
// ##########################################################
// #################    Ext64[2] found   ####################
// ##########################################################
// |Ext64[2]| = 64
// Type        NonAbelianGroup
// BaseGroup   (C2 x C4) . (C2 x C4)
// 
// Elements Orders : [1]:1, [2]:7, [4]:56
// 
// AllSubGr:121 AllConjsCl:97 AllNorms:73
// 
// ##########################################################
// #################    Ext64[3] found   ####################
// ##########################################################
// |Ext64[3]| = 64
// Type        NonAbelianGroup
// BaseGroup   (C2 x C4) . (C2 x C4)
// 
// Elements Orders : [1]:1, [2]:7, [4]:24, [8]:32
// 
// AllSubGr:81 AllConjsCl:71 AllNorms:61
// 
// ##########################################################
// #################    Ext64[4] found   ####################
// ##########################################################
// |Ext64[4]| = 64
// Type        NonAbelianGroup
// BaseGroup   (C2 x C4) . (C2 x C4)
// 
// Elements Orders : [1]:1, [2]:7, [4]:24, [8]:32
// 
// AllSubGr:77 AllConjsCl:59 AllNorms:41
// 
// ##########################################################
// #################    Ext64[5] found   ####################
// ##########################################################
// |Ext64[5]| = 64
// Type        NonAbelianGroup
// BaseGroup   (C2 x C4) . (C2 x C4)
// 
// Elements Orders : [1]:1, [2]:3, [4]:12, [8]:48
// 
// AllSubGr:37 AllConjsCl:33 AllNorms:29
// 
// ##########################################################
// #################    Ext64[6] found   ####################
// ##########################################################
// |Ext64[6]| = 64
// Type        NonAbelianGroup
// BaseGroup   (C2 x C4) . (C2 x C4)
// 
// Elements Orders : [1]:1, [2]:15, [4]:48
// 
// AllSubGr:201 AllConjsCl:165 AllNorms:129
// 
// ##########################################################
// #################    Ext64[7] found   ####################
// ##########################################################
// |Ext64[7]| = 64
// Type        NonAbelianGroup
// BaseGroup   (C2 x C4) . (C2 x C4)
// 
// Elements Orders : [1]:1, [2]:7, [4]:56
// 
// AllSubGr:113 AllConjsCl:89 AllNorms:65
// 
// ##########################################################
// #################    Ext64[8] found   ####################
// ##########################################################
// |Ext64[8]| = 64
// Type        NonAbelianGroup
// BaseGroup   (C2 x C4) . (C2 x C4)
// 
// Elements Orders : [1]:1, [2]:7, [4]:24, [8]:32
// 
// AllSubGr:73 AllConjsCl:59 AllNorms:45
// 
// ##########################################################
// #################    Ext64[9] found   ####################
// ##########################################################
// |Ext64[9]| = 64
// Type        NonAbelianGroup
// BaseGroup   (C2 x C4) . (C2 x C4)
// 
// Elements Orders : [1]:1, [2]:7, [4]:24, [8]:32
// 
// AllSubGr:81 AllConjsCl:73 AllNorms:65
// 
// ##########################################################
// #################   Ext64[10] found   ####################
// ##########################################################
// |Ext64[10]| = 64
// Type        NonAbelianGroup
// BaseGroup   (C2 x C4) . (C2 x C4)
// 
// Elements Orders : [1]:1, [2]:7, [4]:8, [8]:16, [16]:32
// 
// AllSubGr:45 AllConjsCl:33 AllNorms:21
// 
// #  Time:203557 ms
// 

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
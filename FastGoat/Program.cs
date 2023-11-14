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

IEnumerable<ConcreteGroup<Ep2<Tn, Tg>>> AllExtensions<Tn, Tg>(ConcreteGroup<Tn> N, ConcreteGroup<Tg> G, bool bugDetails = false)
    where Tg : struct, IElt<Tg>
    where Tn : struct, IElt<Tn>
{
    var CN = Group.Zentrum(N);
    var autN = Group.AutomorphismGroup(N);
    var ops = Group.AllHomomorphisms(G, autN);
    var allExts = new HashSet<ConcreteGroup<Ep2<Tn, Tg>>>(new IsomorphEquality<Ep2<Tn, Tg>>());
    foreach (var (op, i) in ops.Select((l, i) => (l, i + 1)))
    {
        var L = op.ToMapElt(autN);
        var lbl = $"Lbl{i}/{ops.Count}";
        var (cohs, cobs, cocs) = ZNSolver.Reduce2Cohomologies(CN, G, L, lbl);
        foreach (var c in cohs)
        {
            var c0 = c.ToMapElt;
            var ext = Group.ExtensionGroup(N, L, c0, G);
            var extBase = ((ExtensionGroupBase<Tn, Tg>)ext.BaseGroup);
            if (extBase.IsGroup)
            {
                if (allExts.Add(ext))
                    yield return ext;
            }
            else
            {
                Console.WriteLine("????????????????????? Extension isnt a group");
                if (bugDetails)
                {
                    var gCobs = cobs.allMaps.ToGroupMapElt("B2");
                    var gCocs = cocs.allMaps.ToGroupMapElt("Z2");
                    var gCohs = gCocs.Over(gCobs, "H2");
                    var cs0 = gCohs.First(cs => cs.Contains(c0));
                    var h2dfs = CocyclesDFS.TwoCocycles(N, G, L, lbl);
                    CocyclesDFS.DisplayMapElt("c x", c0, cs0.X);
                    CocyclesDFS.DisplayMapElt("CS0", cs0.OrderMaps(G).ToArray());
                    
                    h2dfs.AllTwoCocycles();
                    var cs1 = h2dfs.AllCosets.Values.First(cs => cs.Contains(c0)).OrderMaps(G).ToArray();
                    CocyclesDFS.DisplayMapElt("CS1", cs1);
                    var ext0 = Group.ExtensionGroup(N, L, cs1.First(), G);
                    var extBase0 = ((ExtensionGroupBase<Tn, Tg>)ext0.BaseGroup);
                    Console.WriteLine($"extBase IsGroup {extBase.IsGroup}/{extBase0.IsGroup}");
                }
            }
        }

        Console.WriteLine($"Nb Exts:{allExts.Count}");
    }
}

{
    var allExts = AllExtensions(N: FG.Abelian(4), G: FG.Abelian(2, 2), bugDetails: true).ToArray(); // BUG
}

void All16Order()
{
    GlobalStopWatch.Restart();
    var allExts = new HashSet<ConcreteGroup<Ep2<Ep<ZnInt>, Ep<ZnInt>>>>(new IsomorphEquality<Ep2<Ep<ZnInt>, Ep<ZnInt>>>());

    allExts.UnionWith(AllExtensions(N: FG.Abelian(8), G: FG.Abelian(2)));
    Console.WriteLine($"Total Exts:{allExts.Count}");
    allExts.UnionWith(AllExtensions(N: FG.Abelian(2, 4), G: FG.Abelian(2)));
    Console.WriteLine($"Total Exts:{allExts.Count}");
    allExts.UnionWith(AllExtensions(N: FG.Abelian(2, 2, 2), G: FG.Abelian(2)));
    Console.WriteLine($"Total Exts:{allExts.Count}");

    CocyclesDFS.DisplayInfosGroups(allExts, naming: true, prefix: "Ext");
    GlobalStopWatch.Show();
    Console.Beep();
}

void All24Order()
{
    GlobalStopWatch.Restart();
    var allExts8 = new HashSet<ConcreteGroup<Ep2<Ep<ZnInt>, Ep<ZnInt>>>>(new IsomorphEquality<Ep2<Ep<ZnInt>, Ep<ZnInt>>>());
    var allExts12 = new HashSet<ConcreteGroup<Ep2<Ep<ZnInt>, Ep<ZnInt>>>>(new IsomorphEquality<Ep2<Ep<ZnInt>, Ep<ZnInt>>>());
    var allExts24 =
        new HashSet<ConcreteGroup<Ep2<Ep2<Ep<ZnInt>, Ep<ZnInt>>, ZnInt>>>(
            new IsomorphEquality<Ep2<Ep2<Ep<ZnInt>, Ep<ZnInt>>, ZnInt>>());

    GlobalStopWatch.AddLap();
    allExts8.UnionWith(AllExtensions(N: FG.Abelian(4), G: FG.Abelian(2)));
    allExts8.UnionWith(AllExtensions(N: FG.Abelian(2, 2), G: FG.Abelian(2)));
    GlobalStopWatch.Show($"allExts8:{allExts8.Count}");

    GlobalStopWatch.AddLap();
    allExts12.UnionWith(AllExtensions(N: FG.Abelian(6), G: FG.Abelian(2)));
    allExts12.UnionWith(AllExtensions(N: FG.Abelian(2, 2), G: FG.Abelian(3)));
    GlobalStopWatch.Show($"allExts12:{allExts12.Count}");

    GlobalStopWatch.AddLap();
    allExts24.UnionWith(allExts8.SelectMany(n => AllExtensions(n, new Cn(3))));
    allExts24.UnionWith(allExts12.SelectMany(n => AllExtensions(n, new Cn(2))));
    GlobalStopWatch.Show($"allExts24:{allExts24.Count}");

    CocyclesDFS.DisplayInfosGroups(allExts8, naming: true, prefix: "Ext");
    CocyclesDFS.DisplayInfosGroups(allExts12, naming: true, prefix: "Ext");
    CocyclesDFS.DisplayInfosGroups(allExts24, naming: true, prefix: "Ext");
    GlobalStopWatch.Show();

    Console.Beep();
}

void All32Order()
{
    GlobalStopWatch.Restart();
    var allExts = new HashSet<ConcreteGroup<Ep2<Ep<ZnInt>, Ep<ZnInt>>>>(new IsomorphEquality<Ep2<Ep<ZnInt>, Ep<ZnInt>>>());

    allExts.UnionWith(AllExtensions(N: FG.Abelian(16), G: FG.Abelian(2)));
    Console.WriteLine($"Total Exts:{allExts.Count}");
    allExts.UnionWith(AllExtensions(N: FG.Abelian(2, 8), G: FG.Abelian(2)));
    Console.WriteLine($"Total Exts:{allExts.Count}");
    allExts.UnionWith(AllExtensions(N: FG.Abelian(4, 4), G: FG.Abelian(2)));
    Console.WriteLine($"Total Exts:{allExts.Count}");
    allExts.UnionWith(AllExtensions(N: FG.Abelian(2, 2, 4), G: FG.Abelian(2)));
    Console.WriteLine($"Total Exts:{allExts.Count}");

    allExts.UnionWith(AllExtensions(N: FG.Abelian(2, 4), G: FG.Abelian(4)));
    Console.WriteLine($"Total Exts:{allExts.Count}");
    allExts.UnionWith(AllExtensions(N: FG.Abelian(2, 2, 2), G: FG.Abelian(4)));
    Console.WriteLine($"Total Exts:{allExts.Count}");

    allExts.UnionWith(AllExtensions(N: FG.Abelian(8), G: FG.Abelian(2, 2)));
    Console.WriteLine($"Total Exts:{allExts.Count}");

    foreach (var ext in AllExtensions(N: FG.Abelian(2, 4), G: FG.Abelian(2, 2)))
    {
        if (allExts.Add(ext) && allExts.Count == 50)
            break;
    }

    foreach (var ext in AllExtensions(N: FG.Abelian(2, 2, 2), G: FG.Abelian(2, 2)))
    {
        if (allExts.Add(ext) && allExts.Count == 51)
            break;
    }

    CocyclesDFS.DisplayInfosGroups(allExts, naming: true, prefix: "Ext");
    GlobalStopWatch.Show();
    Console.Beep();
}

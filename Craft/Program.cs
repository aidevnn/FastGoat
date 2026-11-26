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

XArray<T>[] Permutations<T>(T[] start) where T : struct, IElt<T>
{
    XArray<T>[] AllPermutations = [new()];
    for (int k = 1; k <= start.Length; ++k)
    {
        var k1 = start[k - 1];
        AllPermutations = Enumerable.Range(0, k)
            .SelectMany(c => AllPermutations.Select(l => l.Take(c).Append(k1).Concat(l.Skip(c)).ToXArray())).Distinct()
            .OrderBy(a => a).ToArray();
    }

    return AllPermutations;
}

void RunPolyCyclic(int maxOrd = 24)
{
    GlobalStopWatch.Restart();

    for (int ord = 1; ord <= maxOrd; ++ord)
    {
        if (ord == 64)
            continue;
        foreach (var G in FG.AllGroupsOfOrder(ord).Where(g => g.GroupType == GroupType.NonAbelianGroup))
        {
            if (!Group.IsSoluble(G))
                continue;

            var pols = PolyCyclicSequence<Word>.Decomp(G).DistinctBy(pol => pol.RX.ToXArray())
                .OrderBy(pol => pol.Generators.Count).ToList();

            Console.WriteLine($"########################################## {G.ShortName,-50} nb:{pols.Count}");
            foreach (var pcSeq in pols)
            {
                Console.WriteLine($"##################### {G.ShortName}");
                Console.WriteLine(G.WGbase.Definition);
                Console.WriteLine($"{G} gens:{pcSeq.Generators.Select(e => e.e).ToXArray()}");
                pcSeq.Display();

                foreach (var split in pcSeq.FullSplit().Where(e => e.Count != 1))
                {
                    Console.WriteLine($"##################### {G.ShortName}");
                    pcSeq.Display();
                    Console.WriteLine("---------------------");
                    foreach (var (K, isDirectProduct) in split)
                    {
                        K.Display();
                        Console.WriteLine($"isDirectProduct:{isDirectProduct}");
                        Console.WriteLine();
                    }
                }
            }
        }
    }

    GlobalStopWatch.Show();
}

{
    // RunPolyCyclic();
    // RunPolyCyclic();
    // RunPolyCyclic();
}

void NewRunPolyCyclic(int maxOrd = 24)
{
    GlobalStopWatch.Restart();

    for (int ord = 1; ord <= maxOrd; ord++)
    {
        if (ord == 64)
            continue;
        foreach (var G in FG.AllGroupsOfOrder(ord).Where(g => g.GroupType == GroupType.NonAbelianGroup))
        {
            if (!Group.IsSoluble(G))
                continue;
            Console.WriteLine($"{G.ShortName} {G.WGbase.Definition}");
            // var allSubs = G.AllSubgroups();
            // allSubs.Naming();
            foreach (var pCseq in PCseq<Word>.Decomp(G))
            {
                Console.WriteLine($"        PCseq:{pCseq}");
                // pCseq.Relators();
                // foreach (var sub in allSubs.Where(e => !e.IsTrivial && e.IsProper))
                // {
                //     Console.WriteLine($"    {sub,-20}");
                //     foreach (var H in sub.Conjugates)
                //     {
                //         var indPCseq = pCseq.InducedPCseq(H);
                //         Console.WriteLine(
                //             $"        conjGens:{H.GetGenerators().ToXArray(),-30} IndPCseq:{indPCseq}");
                //         if(!H.IsIsomorphicTo(indPCseq.PCGroup))
                //             throw new();
                //
                //         indPCseq.NormalGroup();
                //     }
                // }
            }
        }
    }

    GlobalStopWatch.Show();
}

{
    NewRunPolyCyclic();
    NewRunPolyCyclic();
    NewRunPolyCyclic();
}

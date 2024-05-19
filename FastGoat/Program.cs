using System.CodeDom;
using System.Collections;
using System.Diagnostics;
using System.Numerics;
using FastGoat.Commons;
using FastGoat.Examples;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.Naming;
using FastGoat.Structures.Subgroups;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Characters;
using FastGoat.UserGroup.GModuleN;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Matrix;
using FastGoat.UserGroup.Perms;
using FastGoat.UserGroup.Polynoms;
using FastGoat.UserGroup.Words;
using static FastGoat.Commons.IntExt;
using static FastGoat.Commons.EnumerableExt;
using FastGoat.UserGroup.Padic;
using FastGoat.UserGroup.Words.Tools;
using GroupRegX = System.Text.RegularExpressions;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");

(int m, int n, int r)[] MetaCyclicSdp(int order)
{
    return IntExt.Dividors(order).Where(d => d > 1)
        .SelectMany(m => FG.MetaCyclicSdpGetR(m, order / m).Select(r => (m, n: order / m, r)))
        .ToArray();
}

void AllGensOfMtCycSdpUpToOrder(int maxOrd, int maxDim = 6)
{
    GlobalStopWatch.Restart();
    var missing = new List<(int, int, int)>();
    var allMtCycSdp = (maxOrd - 5).Range(6).SelectMany(ord => MetaCyclicSdp(ord)).ToArray();

    foreach (var e in allMtCycSdp)
    {
        var mtGL = FG.MetaCyclicSdpMat(e.m, e.n, e.r, maxDim);
        if (mtGL.Count() != 1)
        {
            DisplayGroup.HeadOrders(mtGL);
            var ct = FG.CharacterTable(mtGL);
            foreach (var (mat, k) in mtGL.GetGenerators().Select((mat,k)=>(mat, k+1)))
            {
                Console.WriteLine($"gen{k} of order {ct.Classes.GetClassName(mat)}");
                Console.WriteLine(mat);
            }

            ct.DisplayCells(tableOnly: true);
            var dimCt = ct.AllCharacters.Max(chi => chi[mtGL.Neutral()])!.Value.E[0].Num;
            if (dimCt < mtGL.Neutral().GL.N)
                Console.WriteLine($"Differences {dimCt} < {mtGL.Neutral().GL.N}\n");
            
            continue;
        }

        missing.Add(e);
    }

    var total = allMtCycSdp.Length;
    missing.Println(e => $"M({e.Item1}x:{e.Item2}){e.Item3}", $"Missing:{missing.Count} Found:{total - missing.Count}/{total}");
    GlobalStopWatch.Show("END");
    Console.Beep();
}

{
    AllGensOfMtCycSdpUpToOrder(32);
}
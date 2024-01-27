using System.Collections;
using System.ComponentModel;
using System.Diagnostics;
using System.Numerics;
using FastGoat.Commons;
using FastGoat.Examples;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.Naming;
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

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");

void KTransitivitySn(int n, int k)
{
    GlobalStopWatch.Restart();
    Perm.Style = DisplayPerm.CyclesComma;
    var sn = new Symm(n);
    Console.WriteLine($"Transitives groups of order {n}");
    var Xn = new Cn(n);
    var Xnk = FG.Abelian(Enumerable.Repeat(n,k).ToArray()).ToArray();
    var setXnk = Xnk.Where(xnk => xnk.Ei.Distinct().Count() == k).ToArray();
    Ep<ZnInt> Image(Perm g, Ep<ZnInt> x) => new(x.Ei.Select(e => Xn[g.Table[e.K]]).ToArray());
    var allSubgSn = sn.AllSubgroups();
    var nb = 0;
    foreach (var g1 in allSubgSn.AllRepresentatives)
    {
        var allOrbits = Group.AllOrbits(g1, setXnk, Image);
        if (allOrbits.Count == 1)
        {
            ++nb;
            var subg = allSubgSn.Restriction(g1).ToGroupWrapper();
            g1.Name = NamesTree.BuildName(subg)[0].Name;
            Console.WriteLine(g1.ShortName);
            Console.WriteLine($"Generators:{{ {g1.GetGenerators().Glue(", ")} }}");
        }
    }
    
    GlobalStopWatch.Show($"Total:{nb} transitives groups of order {n}");
    Console.WriteLine();
}

void TransitivitySn(int n) => KTransitivitySn(n, 1);

{
    TransitivitySn(2);
    TransitivitySn(3);
    TransitivitySn(4);
    TransitivitySn(5);
    TransitivitySn(6);
    TransitivitySn(7);
    Console.Beep();
}

/* Transitives groups of order 2
   |C2| = 2
   Generators:{ [(1, 2)] }
   # Total:1 transitives groups of order 2 Time:992ms
   
   Transitives groups of order 3
   |C3| = 3
   Generators:{ [(1, 2, 3)] }
   |S3| = 6
   Generators:{ [(2, 3)], [(1, 2)] }
   # Total:2 transitives groups of order 3 Time:116ms
   
   Transitives groups of order 4
   |C2 x C2| = 4
   Generators:{ [(1, 2), (3, 4)], [(1, 3), (2, 4)] }
   |C4| = 4
   Generators:{ [(1, 2, 3, 4)] }
   |D8| = 8
   Generators:{ [(3, 4)], [(1, 2)], [(1, 3), (2, 4)] }
   |A4| = 12
   Generators:{ [(2, 3, 4)], [(1, 2), (3, 4)] }
   |S4| = 24
   Generators:{ [(3, 4)], [(2, 3)], [(1, 2)] }
   # Total:5 transitives groups of order 4 Time:158ms
   
   Transitives groups of order 5
   |C5| = 5
   Generators:{ [(1, 2, 3, 4, 5)] }
   |D10| = 10
   Generators:{ [(2, 4), (3, 5)], [(1, 2), (3, 4)] }
   |F(5x:4)3| = 20
   Generators:{ [(2, 3, 4, 5)], [(1, 2), (3, 4)] }
   |A5| = 60
   Generators:{ [(3, 4, 5)], [(2, 3), (4, 5)], [(1, 2), (4, 5)] }
   |S5| = 120
   Generators:{ [(4, 5)], [(3, 4)], [(2, 3)], [(1, 2)] }
   # Total:5 transitives groups of order 5 Time:523ms
   
   Transitives groups of order 6
   |S3| = 6
   Generators:{ [(1, 2, 3), (4, 6, 5)], [(1, 4), (2, 5), (3, 6)] }
   |C6| = 6
   Generators:{ [(1, 2, 3, 4, 5, 6)] }
   |A4| = 12
   Generators:{ [(3, 4), (5, 6)], [(1, 2), (5, 6)], [(1, 3, 5), (2, 4, 6)] }
   |D12| = 12
   Generators:{ [(2, 6), (3, 5)], [(1, 2), (3, 4)], [(1, 3), (2, 4), (5, 6)] }
   |M(3x:6)2| = 18
   Generators:{ [(2, 3, 4)], [(1, 2), (3, 5), (4, 6)] }
   |C2 x A4| = 24
   Generators:{ [(5, 6)], [(3, 4)], [(1, 2)], [(1, 3, 5), (2, 4, 6)] }
   |S4| = 24
   Generators:{ [(2, 3, 5, 6)], [(1, 2), (3, 6), (4, 5)] }
   |S4| = 24
   Generators:{ [(3, 4), (5, 6)], [(3, 5), (4, 6)], [(1, 2), (4, 5)], [(1, 3), (2, 6)] }
   |(C3 x C3) x: C4| = 36
   Generators:{ [(3, 4, 5)], [(2, 6), (4, 5)], [(1, 2), (4, 5)], [(1, 3), (2, 4, 6, 5)] }
   |S3 x S3| = 36
   Generators:{ [(3, 5), (4, 6)], [(2, 4, 6)], [(1, 2), (3, 4), (5, 6)] }
   |C2 x S4| = 48
   Generators:{ [(3, 6)], [(2, 3), (5, 6)], [(1, 2), (4, 5)] }
   |A5| = 60
   Generators:{ [(3, 4), (5, 6)], [(2, 3), (4, 6)], [(1, 2), (5, 6)] }
   |(C3 x C3) x: D8| = 72
   Generators:{ [(5, 6)], [(4, 5)], [(2, 3)], [(1, 2)], [(1, 4), (2, 5), (3, 6)] }
   |S5| = 120
   Generators:{ [(3, 4, 6, 5)], [(2, 3), (4, 6)], [(1, 2), (4, 5)] }
   |A6| = 360
   Generators:{ [(4, 5, 6)], [(3, 4), (5, 6)], [(2, 3), (5, 6)], [(1, 2), (5, 6)] }
   |S6| = 720
   Generators:{ [(5, 6)], [(4, 5)], [(3, 4)], [(2, 3)], [(1, 2)] }
   # Total:16 transitives groups of order 6 Time:5.185s
   
   Transitives groups of order 7
   |C7| = 7
   Generators:{ [(1, 2, 3, 4, 5, 6, 7)] }
   |D14| = 14
   Generators:{ [(2, 5), (3, 6), (4, 7)], [(1, 2), (3, 5), (4, 6)] }
   |F(7x:3)2| = 21
   Generators:{ [(2, 4, 6), (3, 5, 7)], [(1, 2, 3), (4, 6, 5)] }
   |F(7x:6)5| = 42
   Generators:{ [(2, 3), (4, 5), (6, 7)], [(2, 4, 6), (3, 5, 7)], [(1, 2), (3, 4), (5, 7)] }
   |SL(3,2)| = 168
   Generators:{ [(3, 4), (6, 7)], [(3, 6), (4, 7)], [(2, 3), (5, 6)], [(1, 2), (6, 7)] }
   |A7| = 2520
   Generators:{ [(5, 6, 7)], [(4, 5), (6, 7)], [(3, 4), (6, 7)], [(2, 3), (6, 7)], [(1, 2), (6, 7)] }
   |S7| = 5040
   Generators:{ [(6, 7)], [(5, 6)], [(4, 5)], [(3, 4)], [(2, 3)], [(1, 2)] }
   # Total:7 transitives groups of order 7 Time:4m28s
 */
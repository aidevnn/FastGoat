using System.Diagnostics;
using System.IO.IsolatedStorage;
using System.Numerics;
using FastGoat.Commons;
using FastGoat.Examples;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Characters;
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

void AllSubGrs<T>(ConcreteGroup<T> g) where T : struct, IElt<T>
{
    GlobalStopWatch.Restart();
    var table = Group.AllSubGroups(g);
    var (conjsRepr, allSubGrs) = (table.Keys.ToList(), table.Values.SelectMany(e => e).ToList());
    var frat = Group.FrattiniSubGroup(allSubGrs, g);
    GlobalStopWatch.Show($"{g.ShortName,-15} Frattini {frat.ShortName,-20} NbClasses {conjsRepr.Count,-5} NbSubGrs {allSubGrs.Count,-5}");
    
    var lattice = Group.SubGroupsLattice(conjsRepr).ToArray();
    foreach (var tower in lattice.OrderBy(t => t.Count).ThenBy(t => t[1].Count()))
        Console.WriteLine("    {0}", tower.Select(g0 => $"{g0}[{g0.Count()}]").Glue("-------"));

    Console.WriteLine();
    DisplayGroup.Head(g.Over(frat));
}

void PermsGroups(int maxK = 3, int maxN = 6)
{
    for (int k = 0; k < maxK; ++k)
    {
        GlobalStopWatch.Restart();
        for (int i = 3; i <= maxN; i++)
        {
            AllSubGrs(FG.Alternate(i));
            AllSubGrs(FG.Symmetric(i));
        }

        Console.WriteLine();
    }
}
/*
   # |Alt3| = 3      NbIsos 2     NbSubGrs 2     Time:0 ms
   # |Symm3| = 6     NbIsos 4     NbSubGrs 6     Time:0 ms
   # |Alt4| = 12     NbIsos 5     NbSubGrs 10    Time:1 ms
   # |Symm4| = 24    NbIsos 9     NbSubGrs 30    Time:19 ms
   # |Alt5| = 60     NbIsos 9     NbSubGrs 59    Time:82 ms
   # |Symm5| = 120   NbIsos 16    NbSubGrs 156   Time:420 ms
   # |Alt6| = 360    NbIsos 16    NbSubGrs 501   Time:8971 ms
   # |Symm6| = 720   NbIsos 29    NbSubGrs 1455  Time:90667 ms
*/

void OthersGroups()
{
    Console.WriteLine();
    for (int i = 3; i <= 10; i++)
        AllSubGrs(FG.Dihedral(i));
    /*
       # |D6| = 6        NbClasses 4     NbSubGrs 6     Time:34 ms
       # |D8| = 8        NbClasses 8     NbSubGrs 10    Time:47 ms
       # |D10| = 10      NbClasses 4     NbSubGrs 8     Time:2 ms
       # |D12| = 12      NbClasses 10    NbSubGrs 16    Time:8 ms
       # |D14| = 14      NbClasses 4     NbSubGrs 10    Time:6 ms
       # |D16| = 16      NbClasses 11    NbSubGrs 19    Time:19 ms
       # |D18| = 18      NbClasses 6     NbSubGrs 16    Time:12 ms
       # |D20| = 20      NbClasses 10    NbSubGrs 22    Time:24 ms
     */

    Console.WriteLine();
    for (int i = 3; i <= 10; i++)
        AllSubGrs(FG.DiCyclic(i));
    /*
       # |Dic3| = 12     NbClasses 6     NbSubGrs 8     Time:42 ms
       # |Dic4| = 16     NbClasses 9     NbSubGrs 11    Time:109 ms
       # |Dic5| = 20     NbClasses 6     NbSubGrs 10    Time:79 ms
       # |Dic6| = 24     NbClasses 12    NbSubGrs 18    Time:202 ms
       # |Dic7| = 28     NbClasses 6     NbSubGrs 12    Time:138 ms
       # |Dic8| = 32     NbClasses 12    NbSubGrs 20    Time:260 ms
       # |Dic9| = 36     NbClasses 9     NbSubGrs 19    Time:284 ms
       # |Dic10| = 40    NbClasses 12    NbSubGrs 24    Time:490 ms
     */

    Console.WriteLine();
    for (int i = 3; i <= 7; i++)
        AllSubGrs(FG.Quaternion(2.Pow(i)));
    /*
       # |Q8| = 8        NbClasses 6     NbSubGrs 6     Time:48 ms
       # |Q16| = 16      NbClasses 9     NbSubGrs 11    Time:11 ms
       # |Q32| = 32      NbClasses 12    NbSubGrs 20    Time:42 ms
       # |Q64| = 64      NbClasses 15    NbSubGrs 37    Time:185 ms
       # |Q128| = 128    NbClasses 18    NbSubGrs 70    Time:649 ms
     */

    var gl23 = new GL(2, 3);
    var GL23mat = Group.Generate("GL2(3)", gl23, gl23[2, 1, 0, 1], gl23[1, 0, 1, 1]);
    AllSubGrs(GL23mat);
    // # |GL2(3)| = 48   NbClasses 16    NbSubGrs 55    Time:62 ms
    var SL23mat = Group.Generate("SL2(3)", gl23, gl23[1, 1, 0, 1], gl23[0, 1, 2, 0]);
    AllSubGrs(SL23mat);
}

{
    // PermsGroups();
    OthersGroups();
    // PermsGroups(1, 5);

    // E = a, b, c | a 4 = b2 = c3 = e, bab = a 3 , bcb = c, aca 3 = c2 .
    // var E = FG.WordGroup("E", "a4, b2, c3, bab = a3, bcb = c, aca3 = c2");
    // DisplayGroup.HeadOrders(E);
    // SemiDirectProduct<Ep2<ZnInt, ZnInt>, ZnInt> E1 = Group
    //     .AllSemiDirectProd("(C6 x C2) x: C2", Product.Generate(new Cn(6), new Cn(2)), new Cn(2)).First(g0 => g0.IsIsomorphicTo(E));
    //
    // AllSubGrs(E);
    // AllSubGrs(E1);
}
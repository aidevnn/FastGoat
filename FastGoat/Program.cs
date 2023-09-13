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
    var (isos, subgs) = (table.Keys.ToList(), table.Values.SelectMany(e => e).ToList());
    GlobalStopWatch.Show($"{g.ShortName,-15} NbClasses {isos.Count,-5} NbSubGrs {subgs.Count,-5}");
    
    var lattice = Group.SubGroupsLattice(isos.ToList()).ToArray();
    foreach (var tower in lattice.OrderBy(t => t.Count).ThenBy(t => t[1].Count()))
    {
        string name(ConcreteGroup<T> g0) =>
            table[g0].Count == 1
                ? $"{g0}[{g0.Count()}]"
                : $"{g0}[{g0.Count()}](x{table[g0].Count})";
        
        Console.WriteLine("    {0}", tower.Select(g0 => name(g0)).Glue("-------"));
    }

    Console.WriteLine();
    var frat = Group.FrattiniSubGroup(subgs, g);
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
       # |D6| = 6        NbClasses 4     NbSubGrs 6     Time:74 ms
       # |D8| = 8        NbClasses 5     NbSubGrs 10    Time:3 ms
       # |D10| = 10      NbClasses 4     NbSubGrs 8     Time:2 ms
       # |D12| = 12      NbClasses 7     NbSubGrs 16    Time:8 ms
       # |D14| = 14      NbClasses 4     NbSubGrs 10    Time:8 ms
       # |D16| = 16      NbClasses 7     NbSubGrs 19    Time:20 ms
       # |D18| = 18      NbClasses 6     NbSubGrs 16    Time:19 ms
       # |D20| = 20      NbClasses 7     NbSubGrs 22    Time:32 ms
     */
    
    Console.WriteLine();
    for (int i = 3; i <= 10; i++)
        AllSubGrs(FG.DiCyclic(i));
    /*
       # |Dic3| = 12     NbClasses 6     NbSubGrs 8     Time:66 ms
       # |Dic4| = 16     NbClasses 6     NbSubGrs 11    Time:60 ms
       # |Dic5| = 20     NbClasses 6     NbSubGrs 10    Time:110 ms
       # |Dic6| = 24     NbClasses 9     NbSubGrs 18    Time:172 ms
       # |Dic7| = 28     NbClasses 6     NbSubGrs 12    Time:169 ms
       # |Dic8| = 32     NbClasses 8     NbSubGrs 20    Time:374 ms
       # |Dic9| = 36     NbClasses 9     NbSubGrs 19    Time:366 ms
       # |Dic10| = 40    NbClasses 9     NbSubGrs 24    Time:807 ms
     */

    Console.WriteLine();
    for (int i = 3; i <= 7; i++)
        AllSubGrs(FG.Quaternion(2.Pow(i)));
    /*
       # |Q8| = 8        NbClasses 4     NbSubGrs 6     Time:51 ms
       # |Q16| = 16      NbClasses 6     NbSubGrs 11    Time:7 ms
       # |Q32| = 32      NbClasses 8     NbSubGrs 20    Time:36 ms
       # |Q64| = 64      NbClasses 10    NbSubGrs 37    Time:205 ms
       # |Q128| = 128    NbClasses 12    NbSubGrs 70    Time:1064 ms
     */
    
    var gl23 = new GL(2, 3);
    var r0 = gl23[2, 1, 0, 1];
    var r1 = gl23[1, 0, 1, 1];
    var GL23mat = Group.Generate("GL2(3)", gl23, r0, r1);
    AllSubGrs(GL23mat);
    // # |GL2(3)| = 48   NbClasses 14    NbSubGrs 55    Time:62 ms
}

{
    // PermsGroups();
    // OthersGroups();
    PermsGroups(1, 5);
}
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
    var (isos, subgs) = Group.AllSubGroups(g);
    GlobalStopWatch.Show($"{g.ShortName,-15} NbIsos {isos.Count,-5} NbSubGrs {subgs.Count,-5}");
}

void PermsGroups()
{
    for (int k = 0; k < 4; ++k)
    {
        GlobalStopWatch.Restart();
        for (int i = 3; i < 7; i++)
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

{
    Console.WriteLine();
    for (int i = 3; i <= 10; i++)
        AllSubGrs(FG.Dihedral(i));
    
    Console.WriteLine();
    for (int i = 3; i <= 10; i++)
        AllSubGrs(FG.DiCyclic(i));

    Console.WriteLine();
    for (int i = 3; i <= 7; i++)
        AllSubGrs(FG.Quaternion(2.Pow(i)));
}

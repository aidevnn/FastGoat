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
    
    GlobalStopWatch.Show(
        $"{g.ShortName,-15} Frattini {frat.ShortName,-20} NbClasses {conjsRepr.Count,-5} NbSubGrs {allSubGrs.Count,-5}");
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

void OthersGroups()
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

    var gl23 = new GL(2, 3);
    var GL23mat = Group.Generate("GL2(3)", gl23, gl23[2, 1, 0, 1], gl23[1, 0, 1, 1]);
    AllSubGrs(GL23mat);
    
    var SL23mat = Group.Generate("SL2(3)", gl23, gl23[1, 1, 0, 1], gl23[0, 1, 2, 0]);
    AllSubGrs(SL23mat);

    Console.WriteLine();
}

{
    OthersGroups();
    PermsGroups(1, 7);
}

/*
   # |D6| = 6        Frattini |Φ(D6)| = 1          NbClasses 4     NbSubGrs 6     Time:42 ms
   # |D8| = 8        Frattini |Φ(D8)| = 2          NbClasses 8     NbSubGrs 10    Time:3 ms
   # |D10| = 10      Frattini |Φ(D10)| = 1         NbClasses 4     NbSubGrs 8     Time:1 ms
   # |D12| = 12      Frattini |Φ(D12)| = 1         NbClasses 10    NbSubGrs 16    Time:6 ms
   # |D14| = 14      Frattini |Φ(D14)| = 1         NbClasses 4     NbSubGrs 10    Time:3 ms
   # |D16| = 16      Frattini |Φ(D16)| = 4         NbClasses 11    NbSubGrs 19    Time:7 ms
   # |D18| = 18      Frattini |Φ(D18)| = 3         NbClasses 6     NbSubGrs 16    Time:7 ms
   # |D20| = 20      Frattini |Φ(D20)| = 1         NbClasses 10    NbSubGrs 22    Time:13 ms
   
   # |Dic3| = 12     Frattini |Φ(Dic3)| = 2        NbClasses 6     NbSubGrs 8     Time:37 ms
   # |Dic4| = 16     Frattini |Φ(Dic4)| = 4        NbClasses 9     NbSubGrs 11    Time:37 ms
   # |Dic5| = 20     Frattini |Φ(Dic5)| = 2        NbClasses 6     NbSubGrs 10    Time:44 ms
   # |Dic6| = 24     Frattini |Φ(Dic6)| = 2        NbClasses 12    NbSubGrs 18    Time:134 ms
   # |Dic7| = 28     Frattini |Φ(Dic7)| = 2        NbClasses 6     NbSubGrs 12    Time:85 ms
   # |Dic8| = 32     Frattini |Φ(Dic8)| = 8        NbClasses 12    NbSubGrs 20    Time:159 ms
   # |Dic9| = 36     Frattini |Φ(Dic9)| = 6        NbClasses 9     NbSubGrs 19    Time:164 ms
   # |Dic10| = 40    Frattini |Φ(Dic10)| = 2       NbClasses 12    NbSubGrs 24    Time:264 ms
   
   # |Q8| = 8        Frattini |Φ(Q8)| = 2          NbClasses 6     NbSubGrs 6     Time:21 ms
   # |Q16| = 16      Frattini |Φ(Q16)| = 4         NbClasses 9     NbSubGrs 11    Time:4 ms
   # |Q32| = 32      Frattini |Φ(Q32)| = 8         NbClasses 12    NbSubGrs 20    Time:28 ms
   # |Q64| = 64      Frattini |Φ(Q64)| = 16        NbClasses 15    NbSubGrs 37    Time:106 ms
   # |Q128| = 128    Frattini |Φ(Q128)| = 32       NbClasses 18    NbSubGrs 70    Time:349 ms
   # |GL2(3)| = 48   Frattini |Φ(GL2(3))| = 2      NbClasses 16    NbSubGrs 55    Time:28 ms
   # |SL2(3)| = 24   Frattini |Φ(SL2(3))| = 2      NbClasses 7     NbSubGrs 15    Time:5 ms
   
   # |Alt3| = 3      Frattini |Φ(Alt3)| = 1        NbClasses 2     NbSubGrs 2     Time:0 ms
   # |Symm3| = 6     Frattini |Φ(Symm3)| = 1       NbClasses 4     NbSubGrs 6     Time:1 ms
   # |Alt4| = 12     Frattini |Φ(Alt4)| = 1        NbClasses 5     NbSubGrs 10    Time:4 ms
   # |Symm4| = 24    Frattini |Φ(Symm4)| = 1       NbClasses 11    NbSubGrs 30    Time:19 ms
   # |Alt5| = 60     Frattini |Φ(Alt5)| = 1        NbClasses 9     NbSubGrs 59    Time:49 ms
   # |Symm5| = 120   Frattini |Φ(Symm5)| = 1       NbClasses 19    NbSubGrs 156   Time:216 ms
   # |Alt6| = 360    Frattini |Φ(Alt6)| = 1        NbClasses 22    NbSubGrs 501   Time:866 ms
   # |Symm6| = 720   Frattini |Φ(Symm6)| = 1       NbClasses 56    NbSubGrs 1455  Time:6123 ms
   # |Alt7| = 2520   Frattini |Φ(Alt7)| = 1        NbClasses 40    NbSubGrs 3786  Time:69673 ms
   # |Symm7| = 5040  Frattini |Φ(Symm7)| = 1       NbClasses 96    NbSubGrs 11300 Time:493585 ms
*/
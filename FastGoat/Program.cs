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
using FastGoat.UserGroup.Words.ToddCoxeter;
using static FastGoat.Commons.IntExt;
using static FastGoat.Commons.EnumerableExt;
using FastGoat.UserGroup.Padic;
using GroupRegX = System.Text.RegularExpressions.Group;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");

void Coincidences1()
{
    // Algebre Tome 1, Daniel Guin – Thomas Hausberger
    // p107 step by step algorithm
    ToddCoxeterAlgo.Run("a", "a3, b3, aba2b", details: true);
    ToddCoxeterAlgo.Run("a3, b3, aba2b", details: true);
}

void Coincidences2()
{
    // Ken Brown paper toddcox.pdf
    ToddCoxeterAlgo.Run("aba-1 = b2, bab-1 = a2", details: true);
}

void Group_576_8282()
{
    // Ken Brown paper toddcox.pdf
    ToddCoxeterAlgo.Run("a,b", "a3,b2,c2,abababab,acac,bcbcbc", details: true);
    Console.ReadLine();
    
    GlobalStopWatch.Restart();
    var g = FG.WordGroup("a3,b2,c2,abababab,acac,bcbcbc");
    DisplayGroup.HeadOrders(g);
    DisplayGroup.HeadOrdersNames(g);
    GlobalStopWatch.Show();
}

/*
   |SL(2,3) x: S4| = 576
   Type        NonAbelianGroup
   BaseGroup   WG[a,b,c]
   
   Elements Orders : [1]:1, [2]:91, [3]:80, [4]:84, [6]:80, [8]:144, [12]:96
   SubGroupsInfos { AllSubGr = 1731, AllConjsCl = 127, AllNorms = 13 }
   Group names
       SL(2,3) x: S4
       (SL(2,3) x: A4) x: C2
       (SL(2,3) x: (C2 x C2)) x: S3
       (Q8 x: A4) x: S3
       (D8 x: (C2 x C2)) x: ((C3 x C3) x: C2)
       C2 . (A4 x: S4)
       Q8 . (C3 x: S4)
   
   #  Time:1m38s
   
 */

void Symm6a()
{
    GlobalStopWatch.Restart();
    ToddCoxeterAlgo.Run("a", "a2, b6, ababababab, ab2ab-2ab2ab-2, abab-1abab-1abab-1", details: true);
    GlobalStopWatch.Show();
}

void Symm6b()
{
    GlobalStopWatch.Restart();
    ToddCoxeterAlgo.Run("b", "a2, b6, ababababab, ab2ab-2ab2ab-2, abab-1abab-1abab-1", details: true);
    GlobalStopWatch.Show();
}

void Symm6()
{
    GlobalStopWatch.Restart();
    ToddCoxeterAlgo.Run("a2, b6, ababababab, ab2ab-2ab2ab-2, abab-1abab-1abab-1", details: true);
    GlobalStopWatch.Show();
}

void Symm6Orders()
{
    GlobalStopWatch.Restart();
    var g = FG.WordGroup("a2, b6, ababababab, ab2ab-2ab2ab-2, abab-1abab-1abab-1");
    DisplayGroup.HeadOrders(g);
    GlobalStopWatch.Show();
}

{
    // Symm6Orders();
}
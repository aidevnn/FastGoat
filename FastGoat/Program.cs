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
using FastGoat.UserGroup.Floats;
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

{
    foreach (var g in FG.AllGroupsOfOrder(8))
        DisplayGroup.HeadCayleyGraph(g);
    
    foreach (var g in FG.AllGroupsOfOrder(12))
        DisplayGroup.HeadCayleyGraph(g);
    
    foreach (var g in FG.AllGroupsOfOrder(16))
        DisplayGroup.HeadElementsCayleyGraph(g);
}

/*
   |Q8| = 8
   Type        NonAbelianGroup
   BaseGroup   WG[a,b]
   
   Cycles with Gen: (5)[4] = a
       (1) --> (5) --> (2) --> (3) --> (1)
       (4) --> (8) --> (6) --> (7) --> (4)
   Cycles with Gen: (6)[4] = b
       (1) --> (6) --> (2) --> (4) --> (1)
       (5) --> (8) --> (3) --> (7) --> (5)
   Nb arrows:16
   
   |D8| = 8
   Type        NonAbelianGroup
   BaseGroup   WG[a,b]
   
   Cycles with Gen: (8)[4] = a
       (1) --> (8) --> (3) --> (7) --> (1)
       (2) --> (5) --> (6) --> (4) --> (2)
   Cycles with Gen: (2)[2] = b
       (1) --> (2) --> (1)
       (7) --> (5) --> (7)
       (8) --> (4) --> (8)
       (3) --> (6) --> (3)
   Nb arrows:16
*/
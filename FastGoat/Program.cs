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

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");

{
    DisplayGroup.HeadOrdersNames(FG.DiCyclicSdp(20));
    DisplayGroup.HeadOrdersNames(FG.GLnp(2, 5));
}
/*
   |Dic20| = 80
   Type        NonAbelianGroup
   BaseGroup   C5 x Q16
   
   Elements Orders : [1]:1, [2]:1, [4]:42, [5]:4, [8]:4, [10]:4, [20]:8, [40]:16
   
   SubGroupsInfos { AllSubGr = 50, AllConjsCl = 18, AllNorms = 11 }
   Group names
       Dic20
       C5 x: Q16
       C10 . D8
       C2 . D40
       C4 . D20
       C40 . C2
       C8 . D10
       Dic10 . C2
       C20 . (C2 x C2)
   
   |GL(2,5)| = 480
   Type        NonAbelianGroup
   BaseGroup   GL(2,5)
   
   Elements Orders : [1]:1, [2]:31, [3]:20, [4]:152, [5]:24, [6]:20, [8]:40, [10]:24, [12]:40, [20]:48, [24]:80
   
   SubGroupsInfos { AllSubGr = 466, AllConjsCl = 48, AllNorms = 6 }
   Group names
       GL(2,5)
       SL(2,5) x: C4
       C4 . S5
       C2 . (A5 x: C4)
       (SL(2,5) x: C2) . C2
*/

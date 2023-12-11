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
    var gl23 = FG.GLnp(2, 3);
    DisplayGroup.HeadOrders(gl23);
    var names = NamesTree.BuildName(gl23);
    names.Println("Group names");
}
/*
   |GL(2,3)| = 48
   Type        NonAbelianGroup
   BaseGroup   GL(2,3)
   
   Elements Orders : [1]:1, [2]:13, [3]:8, [4]:6, [6]:8, [8]:12
   
   Group names
       GL(2,3)
       Q8 x: S3
       SL(2,3) x: C2
       C2 . S4
*/
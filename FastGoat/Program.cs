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
    var g21_1 = FG.MetaCyclicSdpWg(7, 3, 2);
    var (_, mtGL, matSubgrs, names) = GroupMatrixFormPart2.MatrixFormOfGroup(g21_1);
    FG.DisplayName(mtGL, matSubgrs, names, false, false, true, 20);
    GroupMatrixFormPart2.GetCharacter(mtGL, matSubgrs);
    
    var (a, bi) = (g21_1["a"], g21_1["b-1"]);
    DisplayGroup.ElementsCayleyGraph(g21_1, SortBy.Order, a, bi);
}
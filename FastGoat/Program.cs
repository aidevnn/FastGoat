using System.Collections;
using System.ComponentModel;
using System.Diagnostics;
using System.Numerics;
using FastGoat.Commons;
using FastGoat.Examples;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
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
    var ab = FG.Abelian(8, 10, 25);
    DisplayGroup.Head(ab);
    DisplayGroup.Head(FG.AbelianDirectSum(ab).AbCanonic);

    var abt = ab.ToTable();
    DisplayGroup.Head(abt);
    GlobalStopWatch.Bench(5, "Ab", () => Group.AbelianInvariants(ab));
    GlobalStopWatch.Bench(5, "AbTable", () => Group.AbelianInvariants(abt));
    GlobalStopWatch.Bench(5, "Ab", () => Group.AbelianInvariants(ab));
    GlobalStopWatch.Bench(5, "AbTable", () => Group.AbelianInvariants(abt));
    GlobalStopWatch.Bench(5, "Ab", () => Group.AbelianInvariants(ab));
    GlobalStopWatch.Bench(5, "AbTable", () => Group.AbelianInvariants(abt));
    DisplayGroup.AreIsomorphics(ab, abt);
}

{
    DisplayGroup.HeadElements(FG.DihedralWg(4));
}
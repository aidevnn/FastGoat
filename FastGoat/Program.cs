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
    var s3 = ZnInt.ZnZero(3);
    var ma3 = new[] { 1, 1, 0, 1 }.ToKMatrix(2, s3);
    var mb3 = new[] { 0, 1, -1, 0 }.ToKMatrix(2, s3);
    var sl23 = Group.MulGroup("SL(2,3)", ma3, mb3);
    DisplayGroup.Head(sl23);

    var sl23wg = FG.WordGroup("SL(2,3)wg", "a3, abab-1a-1b-1");
    DisplayGroup.Head(sl23wg);

    DisplayGroup.AreIsomorphics(sl23wg, sl23);
    Console.WriteLine();

    var s5 = ZnInt.ZnZero(5);
    var ma5 = new[] { 2, 0, 0, 3 }.ToKMatrix(2, s5);
    var mb5 = new[] { -1, 1, -1, 0 }.ToKMatrix(2, s5);
    var sl25 = Group.MulGroup("SL(2,5)", ma5, mb5);
    DisplayGroup.Head(sl25);
    var sl25wg = FG.WordGroup("SL(2,5)wg", "a5, b3, abababab, ababa-1b-1a-1b-1");
    DisplayGroup.Head(sl25wg);

    DisplayGroup.AreIsomorphics(sl25wg, sl25);
    Console.WriteLine();

    var mc7 = new[] { 1, 3, 2, 5 }.ToKMatrix(2, ZnInt.ZnZero(7));
    DisplayGroup.HeadElements(Group.MulGroup("T", mc7));
}
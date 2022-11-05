using System.Collections;
using System.Diagnostics;
using FastGoat.Commons;
using FastGoat.Examples;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Perms;
using FastGoat.UserGroup.Polynoms;
using FastGoat.UserGroup.Words;
using static FastGoat.Commons.IntExt;
using static FastGoat.Commons.EnumerableExt;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");

{
    DisplayGroup.HeadElements(Group.Galois(16));
    DisplayGroup.HeadElements(Group.Galois(32));
    DisplayGroup.HeadElements(Group.Galois(9));
    DisplayGroup.HeadElements(Group.Galois(7));
}
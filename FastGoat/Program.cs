using System.CodeDom;
using System.Collections;
using System.Diagnostics;
using System.Numerics;
using System.Reflection.Emit;
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
using FastGoat.UserGroup.DatabaseSmallGroups;
using FastGoat.UserGroup.EllCurve;
using FastGoat.UserGroup.Floats;
using FastGoat.UserGroup.GModuleN;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Matrix;
using FastGoat.UserGroup.GModuleN;
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
    var E = new EllGroup<Rational>("-36", "0");
    var O = E.O;
    EllPt<Rational> P = ("-3", "9");
    EllPt<Rational> Q = ("-2", "8");
    Console.WriteLine(new { O, P, Q });
    Console.WriteLine($"-P = {E.Invert(P)}");
    Console.WriteLine($"P + Q = {E.Op(P, Q)}");
    Console.WriteLine($"2P = {E.Times(P, 2)}");
    Console.WriteLine($"2Q = {E.Times(Q, 2)}");
    Console.WriteLine($"2P + 2Q) = {E.Times(E.Op(P, Q), 2)}");
    Console.WriteLine($"2(P + Q) = {E.Op(E.Op(P, Q), E.Op(P, Q))}");
}
using System.Numerics;
using System.Reflection;
using FastGoat.Commons;
using FastGoat.Examples;
using FastGoat.Structures;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.EllCurve;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.LWE;
using FastGoat.UserGroup.Matrix;
using FastGoat.UserGroup.Padic;
using FastGoat.UserGroup.Polynoms;
using Rq = FastGoat.Structures.VecSpace.KPoly<FastGoat.UserGroup.Integers.Rational>;
using static FastGoat.Commons.IntExt;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");
Ring.DisplayPolynomial = MonomDisplay.StarCaret;
// RecomputeAllPrimesUpTo(5000000);
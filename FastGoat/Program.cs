using System.Numerics;
using FastGoat.Commons;
using FastGoat.Examples;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.LWE;
using FastGoat.UserGroup.Polynoms;
using Rq = FastGoat.Structures.VecSpace.KPoly<FastGoat.UserGroup.Integers.Rational>;
using static FastGoat.Commons.IntExt;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");
// RngSeed(259128);
Ring.DisplayPolynomial = MonomDisplay.StarCaret;

{
    var staticVoidMeths = typeof(LearningWithErrors).GetMethods()
        .Where(m => m.IsPublic && m.IsStatic && 
                    m.GetParameters().Length == 0 &&
                    m.ReturnType.Equals(typeof(void)))
        .ToArray();
    
    foreach (var meth in staticVoidMeths)
        meth.Invoke(null, null);
}
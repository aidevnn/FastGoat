using System.Numerics;
using FastGoat.Commons;
using FastGoat.Examples;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Floats;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Polynoms;
using static FastGoat.Commons.IntExt;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");

{
    GlobalStopWatch.Restart();
    Ring.DisplayPolynomial = MonomDisplay.StarCaret;
    
    var staticVoidMeths = typeof(AlgebraicIntegerRelationPSLQ).GetMethods()
        .Where(m => m.IsPublic && m.IsStatic && 
                    m.GetParameters().Length == 0 &&
                    m.ReturnType.Equals(typeof(void)))
        .ToArray();
    
    foreach (var meth in staticVoidMeths)
        meth.Invoke(null, null);
    
    Console.Beep();
}
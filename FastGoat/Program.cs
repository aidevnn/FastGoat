using System.Numerics;
using System.Reflection;
using FastGoat.Commons;
using FastGoat.Examples;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.LWE;
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

{
    var staticVoidMeths = typeof(EllipticCurves).GetMethods()
        .Where(m => m.IsPublic && m.IsStatic &&
                    m.GetParameters().Length == 0 &&
                    m.ReturnType.Equals(typeof(void)))
        .OrderBy(m => m.Name)
        .ToArray();
    
    foreach (var meth in staticVoidMeths)
    {
        Console.WriteLine($"@@@@@@ EllipticCurves.{meth.Name}()");
        meth.Invoke(null, null);
    }
}
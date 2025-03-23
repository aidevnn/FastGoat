using System.Numerics;
using System.Reflection;
using FastGoat.Commons;
using FastGoat.Examples;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
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

{
    var (sk, tk) = (5, 5);
    // 0->4   done ~10min 0 error
    // 5->9   done ~25min 0 error
    // 10->14 done ~15min 0 error
    // 15->19 done ~12min 0 error
    // 20->24 done ~12min 2 planned errors
    // 25->29 done ~2min  0 error
    // 30->34 done ~20min 0 error
    // 35->39 done ~16min 0 error
    // 40->44 done ~40min 0 error
    // 45->49 done ~8min  0 error
    // 50->54 done ~2min  0 error
    // 55->57 done ~8min  0 error
    var examplesList = Assembly.GetExecutingAssembly().GetTypes()
        .Where(t => t.Namespace == "FastGoat.Examples" && t.IsPublic)
        .OrderBy(t => t.FullName)
        .ToList();

    foreach (var (idxType, examplesClass) in examplesList.Index().Skip(sk).Take(tk))
    {
        if (examplesClass.FullName != null && examplesClass.FullName.Contains("CocyclesDFS"))
            continue;

        Console.WriteLine($"@@@@@@[{idxType}] {examplesClass.FullName} @@@@@@");
        Console.WriteLine();

        var staticVoidMeths = examplesClass.GetMethods()
            .Where(m => m.IsPublic && m.IsStatic &&
                        m.GetParameters().Length == 0 &&
                        m.ReturnType.Equals(typeof(void)))
            .OrderBy(m => m.Name)
            .ToArray();

        foreach (var (idxMeth, meth) in staticVoidMeths.Index())
        {
            Console.WriteLine($"@@@@@@[{idxType}/{idxMeth}] {examplesClass.FullName}.{meth.Name}()");
            try
            {
                meth.Invoke(null, null);
            }
            catch (Exception e)
            {
                Console.WriteLine(e);
                Console.WriteLine($"Error at {examplesClass.FullName}.{meth.Name}()");
                Console.WriteLine();
            }
        }
    }

    Console.Beep();
}
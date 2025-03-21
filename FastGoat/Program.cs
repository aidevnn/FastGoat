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
    var (sk, tk) = (0, 1000);
    var examplesList = Assembly.GetExecutingAssembly().GetTypes()
        .Where(t => t.Namespace == "FastGoat.Examples" && t.IsPublic && t.IsClass)
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
}
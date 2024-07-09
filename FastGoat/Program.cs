using System.CodeDom;
using System.Collections;
using System.Diagnostics;
using System.Numerics;
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
using FastGoat.UserGroup.Floats;
using FastGoat.UserGroup.GModuleN;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Matrix;
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
    foreach (var n in 7.Range(2))
    {
        foreach (var p in Primes10000.Take(40))
        {
            var ogl = FG.DPGLnpOrder(n, p);
            if (ogl > FG.MatrixGroupMaxOrder)
                break;

            var dpgl = FG.DPGLnp(n, p);
            Console.WriteLine(dpgl.ShortName);

            if (ogl != dpgl.Count())
                throw new($"{dpgl.ShortName} and ord:{ogl}");
        }

        Console.WriteLine();
    }

    foreach (var n in 7.Range(2))
    {
        foreach (var p in Primes10000.Take(40))
        {
            var osl = FG.DPSLnpOrder(n, p);
            if (osl > FG.MatrixGroupMaxOrder)
                break;

            var dpsl = FG.DPSLnp(n, p);
            Console.WriteLine(dpsl.ShortName);

            if (osl != dpsl.Count())
                throw new($"{dpsl.ShortName} and ord:{osl}");

            if (dpsl.Any(m => m.Det != 1))
            {
                Console.WriteLine(dpsl.First(m => m.Det != 1));
                throw new("det != 1");
            }
        }

        Console.WriteLine();
    }
}

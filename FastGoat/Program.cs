using System.ComponentModel.DataAnnotations;
using System.Diagnostics;
using System.Numerics;
using System.Security.Cryptography;
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
    for (int n = 3; n < 200; n++)
    {
        var un = FG.UnInt(n);
        var seq = AbelianInvariantsFactors.Reduce(un).ToArray();
        var ab = FG.Abelian(seq);
        DisplayGroup.AreIsomorphics(un, ab);
    }
}

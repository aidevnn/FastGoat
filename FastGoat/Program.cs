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
using GroupRegX = System.Text.RegularExpressions.Group;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");

{
    foreach (var q in new[] { 4, 5, 7, 8, 9, 11, 13, 17, 19, 23, 25 })
    {
        Console.WriteLine($"############# Start F{q} #############");
        var x = FG.FqX(q);
        var FqMul = Group.MulGroup($"F{q}*", x);
        var (p, r) = PrimesDec(q).First();
        var gens = r.Range(1).Select(i => x.Pow(i)).ToArray();
        var FqAdd = Group.AddGroup($"F{q}+", gens);
        DisplayGroup.HeadOrdersGenerators(FqAdd);
        DisplayGroup.HeadOrdersGenerators(FqMul);

        var sdp1 = Group.SemiDirectProd(FqAdd, FqMul);
        DisplayGroup.HeadOrders(sdp1);

        var cpr = FG.Abelian(Enumerable.Repeat(p, r).ToArray());
        var cqx = FG.Abelian(q - 1);
        var sdp2 = Group.SemiDirectProd(cpr, cqx);
        DisplayGroup.HeadOrders(sdp2);
        DisplayGroup.AreIsomorphics(sdp1, sdp2);
        Console.WriteLine();
        Console.WriteLine($"############# End   F{q} #############");
        Console.WriteLine();
    }
}
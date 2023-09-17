using System.Diagnostics;
using System.IO.IsolatedStorage;
using System.Numerics;
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
    var gl = new GL(2, 3);
    var a0 = gl[2, 0, 0, 1];
    var b0 = gl[2, 1, 2, 0];

    var a1 = gl[1, 1, 0, 1];
    var b1 = gl[0, 1, 2, 0];

    var gl23 = Group.Generate(gl, a0, b0);
    var sl23 = Group.Generate("SL(2,3)", gl, a1, b1);
    var q8 = Group.IsomorphicSubgroup(sl23, FG.Quaternion(8), "Q8");
    DisplayGroup.Head(sl23.Over(q8, "C3"));

    Group.AllSylowPSubgroups(gl23).Select(e => $"{e.Key.ShortName} NbConjs {e.Value.Count}").Println();
    Group.AllSylowPSubgroups(sl23).Select(e => $"{e.Key.ShortName} NbConjs {e.Value.Count}").Println();
    Group.AllSylowPSubgroups(Group.SemiDirectProd(new Cn(5), new Cn(4))).Select(e => $"{e.Key.ShortName} NbConjs {e.Value.Count}")
        .Println();
    Group.AllSylowPSubgroups(Product.Generate(new Cn(2), FG.Dihedral(4))).Select(e => $"{e.Key.ShortName} NbConjs {e.Value.Count}")
        .Println();
    Group.AllSylowPSubgroups(FG.Symmetric(5)).Select(e => $"{e.Key.ShortName} NbConjs {e.Value.Count}").Println();
    
    Group.AllSylowPSubgroups(Product.Generate(new Cn(3), FG.Dihedral(5))).Select(e => $"{e.Key.ShortName} NbConjs {e.Value.Count}")
        .Println();
    Group.AllSylowPSubgroups(Group.SemiDirectProd(FG.Abelian(3, 3), FG.Abelian(8))).Select(e => $"{e.Key.ShortName} NbConjs {e.Value.Count}")
        .Println();
    Group.AllSylowPSubgroups(FG.Abelian(6, 8, 18)).Select(e => $"{e.Key.ShortName} NbConjs {e.Value.Count}")
        .Println();
}

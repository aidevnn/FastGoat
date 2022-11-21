using System.Diagnostics;
using System.Runtime.Intrinsics.X86;
using System.Threading.Channels;
using FastGoat.Commons;
using FastGoat.Examples;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Matrix;
using FastGoat.UserGroup.Perms;
using FastGoat.UserGroup.Polynoms;
using FastGoat.UserGroup.Words;
using FastGoat.UserGroup.Words.ToddCoxeter;
using static FastGoat.Commons.IntExt;
using static FastGoat.Commons.EnumerableExt;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");

{
    var V = Product.Generate(new Cn(2), new Cn(2));
    var autV = Group.AutomorphismGroup(V);
    DisplayGroup.AreIsomorphics(autV, new Symm(3));
    var sdp = Group.SemiDirectProd(V, autV);
    DisplayGroup.AreIsomorphics(new Symm(4), sdp);
}
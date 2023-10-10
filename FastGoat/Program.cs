using System.Collections;
using System.ComponentModel;
using System.Diagnostics;
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
    var G = FG.Abelian(4, 6, 18);
    var abInvs = Group.AbelianInvariants(G);
    var Gs = abInvs.Select(e => Group.Generate($"C{e.o}", G, e.g)).ToArray();
    Gs.Select(g => $"{g.ShortName,-40} generator {g.GetGenerators().First()}").Println($"{G.ShortName}");
}

{
    var G = FG.Abelian(20, 30);
    var abInvs = Group.AbelianInvariants(G);
    var Gs = abInvs.Select(e => Group.Generate($"C{e.o}", G, e.g)).ToArray();
    Gs.Select(g => $"{g.ShortName,-40} generator {g.GetGenerators().First()}").Println($"{G.ShortName}");
}

{
    var G = FG.Abelian(14, 21);
    var abInvs = Group.AbelianInvariants(G);
    var Gs = abInvs.Select(e => Group.Generate($"C{e.o}", G, e.g)).ToArray();
    Gs.Select(g => $"{g.ShortName,-40} generator {g.GetGenerators().First()}").Println($"{G.ShortName}");
}

{
    var G = FG.Abelian(8, 18, 30);
    var abInvs = Group.AbelianInvariants(G);
    var Gs = abInvs.Select(e => Group.Generate($"C{e.o}", G, e.g)).ToArray();
    Gs.Select(g => $"{g.ShortName,-40} generator {g.GetGenerators().First()}").Println($"{G.ShortName}");
}

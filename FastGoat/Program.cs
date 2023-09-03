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
    var x = FG.ZPoly(3);
    var Q = x.Pow(3) + 2 * x.Pow(2) - x - 1;
    var f0 = IntFactorisation.Firr(Q, x.KOne * 2).ToArray();
    f0.Println($"Fact({Q}) irreductible in F3");
    var (xa, a) = FG.EPolyXc(Q, 'a');
    var gr = Group.MulGroup("GF(27)", a, a - 1, a - 2);
    var c0 = gr.GetGenerators().First();
    DisplayGroup.Head(gr);
    var facts = IntFactorisation.Firr(Q.Substitute(xa), c0).ToArray();
    var roots = facts.Select(f => -f[0]).ToList();
    GaloisTheory.GaloisGroup(roots, details: true);
}

{
    var x = FG.ZPoly(3);
    var Q = x.Pow(4) + 2 * x.Pow(3) + 2 * x.Pow(2) + x + 2;
    var f0 = IntFactorisation.Firr(Q, x.KOne * 2).ToArray();
    f0.Println($"Fact({Q}) irreductible in F3");

    var (xc, c) = FG.EPolyXc(Q, 'c');
    var gr = Group.MulGroup("GF(81)", c, c - 1, c - 2);
    var c0 = gr.GetGenerators().First();
    DisplayGroup.Head(gr);
    var facts = IntFactorisation.Firr(Q.Substitute(xc), c0).ToArray();
    var roots = facts.Select(f => -f[0]).ToList();
    GaloisTheory.GaloisGroup(roots, details: true);
}

{
    var (x, a) = FG.FqX_Poly(81);
    var Q = x.Pow(4) + 2 * x.Pow(3) + 2 * x.Pow(2) + x + 2;
    var facts = IntFactorisation.Firr(Q, a).ToArray();
    facts.Println($"Fact({Q}) in F81 ~ F3(a) with {a.F} = 0");
}

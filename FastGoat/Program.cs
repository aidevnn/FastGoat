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
Ring.DisplayPolynomial = MonomDisplay.StarPowFct;
void SymbolicSolve()
{
    var Fp = FG.UnInt(19);
    var gen = Fp.GetGenerators().First();
    var xs = Ring.Polynomial(gen, MonomOrder.Lex, (6, "x")).Deconstruct();
    var bs = new PolynomialBasis<ZnInt, Xi>(xs.SelectMany(P => new[] { P.Pow(9) - 1, P.Pow(19) - 1 }).ToArray());
    var Xs = xs.Select(P => new EPolynomial<ZnInt>(P, P.One, bs)).Deconstruct();
    var (z, o) = (Xs[0].Zero, Xs[0].One);

    var M0 = 6.Range().Grid2D().Select(e => e.t1 == e.t2 ? Xs[e.t1] : z).ToKMatrix(6);
    Console.WriteLine("M0");
    Console.WriteLine(M0);

    Console.WriteLine($"M0^9 = I6 : {M0.Pow(9).Equals(M0.One)}");
    Console.WriteLine();
    var M0_2 = M0.Pow(2);

    var gl6 = new GL(6, 2);
    var a0 = gl6[
        0, 1, 0, 0, 0, 0,
        1, 0, 0, 0, 0, 0,
        0, 0, 1, 0, 0, 0,
        0, 0, 0, 1, 0, 0,
        0, 0, 0, 0, 1, 0,
        0, 0, 0, 0, 0, 1
    ];
    var a1 = gl6[
        0, 1, 0, 0, 0, 0,
        0, 0, 1, 0, 0, 0,
        0, 0, 0, 1, 0, 0,
        0, 0, 0, 0, 1, 0,
        0, 0, 0, 0, 0, 1,
        1, 0, 0, 0, 0, 0
    ];
    
    var s6GL = Group.Generate("S6", gl6, a0, a1);
    var ord6 = s6GL.Where(e => s6GL.ElementsOrders[e] == 6)
        .Select(e => Group.Generate("C6", gl6, e))
        .ToHashSet(new GroupSetEquality<Mat>());

    var sys0 = (M0.Pow(9) - M0.One).Where(e => !e.IsZero()).Select(e => e.Num).ToArray();
    // foreach (var e in ord6.Shuffle().SelectMany(c6 => c6.GetGenerators()))
    foreach (var e in ord6.SelectMany(c6 => c6.GetGenerators()))
    {
        var M1 = e.Table.Select(c => c * o).ToKMatrix(6);
        var Sys = (M1.Inv() * M0 * M1 - M0_2).Where(P => !P.IsZero()).Select(P => P.Num).Concat(sys0).ToArray();
        Console.WriteLine("M1");
        Console.WriteLine(M1);
        Sys.Println("Sys in GL(6, 19), M0^9=Id, M1^6=Id, and M1^-1 * M0 * M1 = M1^2");
        var Sols = Ring.ReducedGrobnerBasis(Sys);
        Sols.Select(P => new EPolynomial<ZnInt>(P, P.One, bs)).Println("Sols");
        Console.WriteLine();
        if (Sols.Any(P => P.Degree > 1))
            break;
    }
}

// M0
// [x0,  0,  0,  0,  0,  0]
// [ 0, x1,  0,  0,  0,  0]
// [ 0,  0, x2,  0,  0,  0]
// [ 0,  0,  0, x3,  0,  0]
// [ 0,  0,  0,  0, x4,  0]
// [ 0,  0,  0,  0,  0, x5]
// M0^9 = I6 : True
// 
// M1
// [ 0,  1,  0,  0,  0,  0]
// [ 0,  0,  1,  0,  0,  0]
// [ 0,  0,  0,  1,  0,  0]
// [ 0,  0,  0,  0,  1,  0]
// [ 0,  0,  0,  0,  0,  1]
// [ 1,  0,  0,  0,  0,  0]
// Sys in GL(6, 19), M0^9=Id, M1^6=Id, and M1^-1 * M0 * M1 = M1^2
//     18*x0.Pow(2) + x5
//     x0 + 18*x1.Pow(2)
//     x1 + 18*x2.Pow(2)
//     x2 + 18*x3.Pow(2)
//     x3 + 18*x4.Pow(2)
//     x4 + 18*x5.Pow(2)
// Sols
//      0
//     (x4 + 18*x5.Pow(2))
//     (x3 + 18*x5.Pow(4))
//     (x2 + 18*x5.Pow(8))
//     (x1 + 18*x5.Pow(7))
//     (x0 + 18*x5.Pow(5))
// 

{
    SymbolicSolve();
    
    var Fp = FG.UnInt(19);
    var x5 = Fp.First(e => Fp.ElementsOrders[e] == 9);
    var x0 = x5.Pow(5);
    var x1 = x5.Pow(7);
    var x2 = x5.Pow(8);
    var x3 = x5.Pow(4);
    var x4 = x5.Pow(2);
    
    var gl6 = new GL(6, 19);
    var a0 = gl6[
        x0.K, 0, 0, 0, 0, 0,
        0, x1.K, 0, 0, 0, 0,
        0, 0, x2.K, 0, 0, 0,
        0, 0, 0, x3.K, 0, 0,
        0, 0, 0, 0, x4.K, 0,
        0, 0, 0, 0, 0, x5.K
    ];
    
    var a1 = gl6[
        0, 1, 0, 0, 0, 0,
        0, 0, 1, 0, 0, 0,
        0, 0, 0, 1, 0, 0,
        0, 0, 0, 0, 1, 0,
        0, 0, 0, 0, 0, 1,
        1, 0, 0, 0, 0, 0
    ];

    var mtSdp = FG.MetaCyclicSdp(9, 6, 2);
    var mtGL = Group.Generate("M(9x:6)2", gl6, a0, a1);
    DisplayGroup.HeadGenerators(mtGL);
    DisplayGroup.AreIsomorphics(mtSdp, mtGL);
}

// |M(9x:6)2| = 54
// Type        NonAbelianGroup
// BaseGroup   GL(6,19)
// 
// Generators of M(9x:6)2
// gen1 of order 6
// [ 0,  1,  0,  0,  0,  0]
// [ 0,  0,  1,  0,  0,  0]
// [ 0,  0,  0,  1,  0,  0]
// [ 0,  0,  0,  0,  1,  0]
// [ 0,  0,  0,  0,  0,  1]
// [ 1,  0,  0,  0,  0,  0]
// gen2 of order 9
// [17,  0,  0,  0,  0,  0]
// [ 0,  6,  0,  0,  0,  0]
// [ 0,  0,  5,  0,  0,  0]
// [ 0,  0,  0,  9,  0,  0]
// [ 0,  0,  0,  0, 16,  0]
// [ 0,  0,  0,  0,  0,  4]
// 
// MtCyc(9,6,2) IsIsomorphicTo M(9x:6)2 : True
// 
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

// {
//     Ring.DisplayPolynomial = MonomDisplay.Caret;
//     var x = FG.QPoly('X');
//     var (a0, _, b0, _) = IntFactorisation.SplittingField(x.Pow(4) - 2 * x.Pow(2) - 2).Deconstruct();
//     var minPoly = a0.F;
//
//     var subFields = GaloisTheory.SubFields(minPoly).ToArray();
//     var extTowers = GaloisCorrespondenceExamples.ExtensionsTower(subFields);
//     GaloisCorrespondenceExamples.GaloisCorrespondence(extTowers);
//
//     var y = subFields.First().primElt.X;
//     var (a, b) = (a0.Substitute(y), b0.Substitute(y));
//     GaloisCorrespondenceExamples.FindExtension(subFields, a.One, "Q");
//     
//     GaloisCorrespondenceExamples.FindExtension(subFields, a.Pow(2), "Q(a^2)");
//     GaloisCorrespondenceExamples.FindExtension(subFields, a.Pow(3) * b - a * b, "Q(a^3*b - a*b)");
//     GaloisCorrespondenceExamples.FindExtension(subFields, a * b, "Q(a*b)");
//     
//     GaloisCorrespondenceExamples.FindExtension(subFields, a, "Q(a)");
//     GaloisCorrespondenceExamples.FindExtension(subFields, b, "Q(b)");
//     GaloisCorrespondenceExamples.FindExtension(subFields, a + b, "Q(a+b)");
//     GaloisCorrespondenceExamples.FindExtension(subFields, a - b, "Q(a-b)");
//
//     var primElt_a2_ab = GaloisTheory.PrimitiveEltComb(a.Pow(2), a * b).W;
//     GaloisCorrespondenceExamples.FindExtension(subFields, primElt_a2_ab, "Q(a^2, a*b)");
//     
//     var primElt_a_b = GaloisTheory.PrimitiveEltComb(a, b).W;
//     GaloisCorrespondenceExamples.FindExtension(subFields, primElt_a_b, "Q(a, b)");
// }
//
// {
//     Ring.DisplayPolynomial = MonomDisplay.Caret;
//     var x = FG.QPoly('X');
//     var P = x.Pow(10) + 10 * x.Pow(8) + 125 * x.Pow(6) + 500 * x.Pow(4) + 2500 * x.Pow(2) + 4000;
//     // var P = x.Pow(12) + 96 * x.Pow(8) + 1664 * x.Pow(6) - 16128 * x.Pow(4) + 165888 * x.Pow(2) + 331776;
//     // var P = x.Pow(12) + 6 * x.Pow(8) + 26 * x.Pow(6) - 63 * x.Pow(4) + 162 * x.Pow(2) + 81;
//     // var P = x.Pow(18) + 171 * x.Pow(12) + 5130 * x.Pow(6) + 27;
//     // var P = x.Pow(12) - 572 * x.Pow(6) + 470596;
//     var (X, y) = FG.EPolyXc(P, 'y');
//     var (P0, X0) = IntFactorisation.Deflate(P, 3);
//     var roots = IntFactorisation.AlgebraicRoots(P0.Substitute(X), true)
//         .SelectMany(r => IntFactorisation.AlgebraicRoots(X0.Substitute(X) - r, true))
//         .ToList();
//     
//     var subFields = GaloisTheory.SubFields(roots, nbGens: 3).ToArray();
//     var extTowers = GaloisCorrespondenceExamples.ExtensionsTower(subFields);
//     GaloisCorrespondenceExamples.GaloisCorrespondence(extTowers);
// }

{
    Ring.DisplayPolynomial = MonomDisplay.Caret;
    var x = FG.QPoly('X');
    
    var P = x.Pow(6) + 243;
    var (X, y) = FG.EPolyXc(P, 'y');
    var rootsK = IntFactorisation.AlgebraicRoots(P);
    rootsK.Println($"P = {P}");
    Console.WriteLine(rootsK.Aggregate(X.One, (prod, r) => prod * (X - r)));
    Console.WriteLine();
    
    var xc = FG.CplxPoly('X');
    var Pc = P.Substitute(xc);
    var rootsC = FG.NRoots(Pc);
    rootsC.Select(e => xc - e).Println($"P = {Pc}");
    Console.WriteLine(rootsC.Aggregate(xc.One, (prod, r) => prod * (xc - r)));
    Console.WriteLine();

    foreach (var yc in rootsK.Order())
    {
        var map = rootsC.Order().ToDictionary(e => e, e => yc.Poly.Substitute(e));
        map.Println($"With yc = {yc}");
        Console.WriteLine($"Bijection : {map.Keys.ToHashSet().SetEquals(map.Values)}");
        Console.WriteLine();
    }
}

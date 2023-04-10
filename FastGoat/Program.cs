using System.ComponentModel.DataAnnotations;
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
    Ring.DisplayPolynomial = MonomDisplay.Caret;
    var x = FG.QPoly('X');
    var (a0, _, b0, _) = IntFactorisation.SplittingField(x.Pow(4) - 2 * x.Pow(2) - 2).Deconstruct();
    var minPoly = a0.F;

    var subFields = GaloisTheory.SubFields(minPoly).ToArray();
    var extTowers = GaloisCorrespondenceExamples.ExtensionsTower(subFields);
    GaloisCorrespondenceExamples.GaloisCorrespondence(extTowers);

    var y = subFields.First().primElt.X;
    var (a, b) = (a0.Substitute(y), b0.Substitute(y));
    GaloisCorrespondenceExamples.FindExtension(subFields, a.One, "Q");
    
    GaloisCorrespondenceExamples.FindExtension(subFields, a.Pow(2), "Q(a^2)");
    GaloisCorrespondenceExamples.FindExtension(subFields, a.Pow(3) * b - a * b, "Q(a^3*b - a*b)");
    GaloisCorrespondenceExamples.FindExtension(subFields, a * b, "Q(a*b)");
    
    GaloisCorrespondenceExamples.FindExtension(subFields, a, "Q(a)");
    GaloisCorrespondenceExamples.FindExtension(subFields, b, "Q(b)");
    GaloisCorrespondenceExamples.FindExtension(subFields, a + b, "Q(a+b)");
    GaloisCorrespondenceExamples.FindExtension(subFields, a - b, "Q(a-b)");

    var primElt_a2_ab = GaloisTheory.PrimitiveEltComb(a.Pow(2), a * b).W;
    GaloisCorrespondenceExamples.FindExtension(subFields, primElt_a2_ab, "Q(a^2, a*b)");
    
    var primElt_a_b = GaloisTheory.PrimitiveEltComb(a, b).W;
    GaloisCorrespondenceExamples.FindExtension(subFields, primElt_a_b, "Q(a, b)");
}

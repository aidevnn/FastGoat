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

void CheckChebotarev(KPoly<Rational> P, int nb, ConcreteGroup<Perm> gal, bool detail = false)
{
    var lt = IntFactorisation.ChebotarevTypes(P, nb, detail)
        .GroupBy(e => e.Deconstruct())
        .ToDictionary(e => e.Key, e => e.Count());

    var types = gal.Select(perm => PermutationToCycles(perm.Sn.N, perm.Table).Select(l => l.Length).Order().Deconstruct())
        .GroupBy(e => e).ToDictionary(e => e.Key, e => e.Count());
    types.Keys.OrderBy(e => e.ToString()).ToDictionary(e => e, e => types[e]).Println($"Expected types");

    var d = gal.Count();
    lt.Keys.OrderBy(e => e.ToString()).ToDictionary(e => e, e => (1.0 * d * lt[e]) / nb).Println("Actual types");

    Console.WriteLine();
}

{
    Ring.DisplayPolynomial = MonomDisplay.Caret;
    var x = FG.QPoly('X');
    
    {
        var P = x.Pow(3) - 3 * x - 1;
        var rootsK = IntFactorisation.AlgebraicRoots(P);
        var gal = GaloisTheory.GaloisGroup(rootsK, details: true);
    
        CheckChebotarev(P, 20, gal);
    }
    
    {
        var P = x.Pow(4) + x.Pow(3) + x.Pow(2) + x + 1;
        var rootsK = IntFactorisation.AlgebraicRoots(P);
        var gal = GaloisTheory.GaloisGroup(rootsK, details: true);
    
        CheckChebotarev(P, 30, gal);
    }
    
    {
        var P = x.Pow(5) + x.Pow(4) - 4 * x.Pow(3) - 3 * x.Pow(2) + 3 * x + 1;
        var rootsK = IntFactorisation.AlgebraicRoots(P);
        var gal = GaloisTheory.GaloisGroup(rootsK, details: true);
    
        CheckChebotarev(P, 20, gal);
    }
    
    {
        var P = x.Pow(6) + 243;
        var rootsK = IntFactorisation.AlgebraicRoots(P);
        var gal = GaloisTheory.GaloisGroup(rootsK, details: true);
    
        CheckChebotarev(P, 30, gal);
    }
    
    {
        var P = x.Pow(8) + 1;
        var rootsK = IntFactorisation.AlgebraicRoots(P);
        var gal = GaloisTheory.GaloisGroup(rootsK, details: true);
    
        CheckChebotarev(P, 30, gal);
    }

    {
        var P = x.Pow(4) + x + 1;
        CheckChebotarev(P, 60, new Symm(4), true); // slow
    }
    
    {
        var P = x.Pow(5) + 2;
        var s5 = new Sn(5);
        var gal = Group.Generate("C5 x: C4", s5, s5[(2, 3, 5, 4)], s5[(1, 2, 3, 4, 5)]);
        DisplayGroup.HeadElements(gal);
        CheckChebotarev(P, 60, gal, true); // slow
    }
}

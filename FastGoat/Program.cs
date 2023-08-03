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
Ring.DisplayPolynomial = MonomDisplay.StarCaret;

// JS Milne, FT v5.10, page 50-51
void GalGrQuarticPol(KPoly<Rational> P)
{
    if (P.Degree != 4 || !P.LT.Equals(P.KOne))
        throw new();

    var Pfacts = IntFactorisation.FirrZ2(P);
    if (Pfacts.Length != 1)
        throw new();

    Console.WriteLine($"P = {P}");
    Console.WriteLine("Chebotarev Gal(P) = {0}", GaloisApplicationsPart2.GaloisGroupChebotarev(P).Name);
    
    var (_, b, c, d, e) = P.Coefs.Reverse().Deconstruct();
    var x = P.X;
    var Q = x.Pow(3) - c * x.Pow(2) + (b * d - 4 * e) * x - b.Pow(2) * e + 4 * c * e - d.Pow(2);
    var Qfacts = IntFactorisation.FirrZ2(Q);
    if (Qfacts.Length != 1)
        Console.WriteLine($"Res(P) = {Q} = {Qfacts.Glue("*", "({0})")}");
    else
        Console.WriteLine($"Res(P) = {Q} is irreductible");
    
    var QfactsDeg2or3 = Qfacts.Where(e => e.Degree > 1).ToArray();
    if (QfactsDeg2or3.Length == 0)
    {
        Console.WriteLine("Gal(P) = V4");
        Console.WriteLine();
        return;
    }

    var Q0 = QfactsDeg2or3[0];
    var Q0facts = IntFactorisation.AlgebraicFactors(Q0);
    Console.WriteLine($"{Q0} = {Q0facts.Glue("*", "({0})")} in Q(y)");
    
    if (Q0facts.Count != Q0.Degree)
        Console.WriteLine("Gal(P) = S4");
    else if (Q0facts.Count == 3)
        Console.WriteLine("Gal(P) = A4");
    else if (Q0facts.Count == 2)
    {
        var P0 = P.Substitute(Q0facts[0].X);
        var P0facts = IntFactorisation.AlgebraicFactors(P0);
        if (P0facts.Count != 1)
            Console.WriteLine($"{P0} = {P0facts.Glue("*", "({0})")} in Q(y)");
        else
            Console.WriteLine($"{P0} is irreductible in Q(y)");
        
        if (P0facts.Count == 1)
            Console.WriteLine("Gal(P) = D4");
        else
            Console.WriteLine("Gal(P) = C4");
    }

    Console.WriteLine();
}

{
    var x = FG.QPoly();
    GalGrQuarticPol(x.Pow(4) - 4 * x.Pow(2) + 1);
    GalGrQuarticPol(x.Pow(4) - 10 * x.Pow(2) + 4);
    GalGrQuarticPol(x.Pow(4) + 1);
    
    GalGrQuarticPol(x.Pow(4) + x.Pow(3) + x.Pow(2) + x + 1);
    GalGrQuarticPol(x.Pow(4) + 4 * x.Pow(2) + 2);
    GalGrQuarticPol(x.Pow(4) - 5 * x.Pow(2) + 5);

    GalGrQuarticPol(x.Pow(4) - 2);
    GalGrQuarticPol(x.Pow(4) - 2 * x.Pow(2) + 3);
    GalGrQuarticPol(x.Pow(4) + x.Pow(2) - 1);
    
    GalGrQuarticPol(x.Pow(4) + 8 * x + 12);

    GalGrQuarticPol(x.Pow(4) - 4 * x + 2);
    GalGrQuarticPol(x.Pow(4) + x + 1);
    GalGrQuarticPol(x.Pow(4) + x.Pow(3) - 5);
}
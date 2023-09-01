using System.ComponentModel.DataAnnotations;
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
//
// {
//     var x = Ring.Polynomial(Rational.KZero(), "x")[0];
//     Console.WriteLine(Ring.LcmPolynomial((x + 1) * (x + 5).Pow(3) * (x - 6), (x + 2) * (x + 5).Pow(2)));
//     Console.WriteLine((x + 1) * (x + 5).Pow(3) * (x + 2) * (x - 6));
//     Console.WriteLine();
//     Console.WriteLine(Ring.GcdPolynomial((x + 1) * (x + 5).Pow(3) * (x - 6), (x + 2) * (x + 5).Pow(2)));
//     Console.WriteLine((x + 5).Pow(2));
//     Console.WriteLine();
// }
//
// {
//     var x = FG.QPoly();
//     Console.WriteLine(Ring.Lcm((x + 1) * (x + 5).Pow(3) * (x - 6), (x + 2) * (x + 5).Pow(2)).Monic);
//     Console.WriteLine((x + 1) * (x + 5).Pow(3) * (x + 2) * (x - 6));
//     Console.WriteLine();
//     Console.WriteLine(Ring.Gcd((x + 1) * (x + 5).Pow(3) * (x - 6), (x + 2) * (x + 5).Pow(2)).Monic);
//     Console.WriteLine((x + 5).Pow(2));
//     Console.WriteLine();
// }
//
// {
//     var (x, y) = Ring.Polynomial(Rational.KZero(), "x", "y").Deconstruct();
//     Console.WriteLine(Ring.LcmPolynomial((x + 1) * (y + 5).Pow(3) * (x - 6), (x + 2) * (y + 5).Pow(2)));
//     Console.WriteLine(Ring.LcmPolynomial((x + 2) * (y + 5).Pow(2), (x + 1) * (y + 5).Pow(3) * (x - 6)));
//     Console.WriteLine((x + 1) * (y + 5).Pow(3) * (x + 2) * (x - 6));
//     Console.WriteLine();
//     Console.WriteLine(Ring.GcdPolynomial((x + 1) * (y + 5).Pow(3) * (x - 6), (x + 2) * (y + 5).Pow(2)));
//     Console.WriteLine(Ring.GcdPolynomial((x + 2) * (y + 5).Pow(2), (x + 1) * (y + 5).Pow(3) * (x - 6)));
//     Console.WriteLine((y + 5).Pow(2));
//     Console.WriteLine();
// }

{

    Ring.DisplayPolynomial = MonomDisplay.StarCaret;
    Console.WriteLine("Swinnerton–Dyer polynomials");
    for (int m = 1; m <= 4; m++)
    {
        var (X0, xis) = Ring.Polynomial(Rational.KZero(), MonomOrder.Graded, (m, "a"), "X");
        var basis = Primes10000.Take(m).Select((i, k) => xis[k].Pow(2) - i).ToArray();
        var B = new PolynomialBasis<Rational, Xi>(basis);
        var Xis = xis.Select(xi => new EPolynomial<Rational>(xi, B)).ToArray();
        var X = new EPolynomial<Rational>(X0, B);
        var Sm = YieldAllCombinations(m).Select(i => i.Select((j, k) => j ? -Xis[k] : Xis[k]).Aggregate(X, (acc, a1) => acc + a1))
            .Aggregate((a0, a1) => a0 * a1);

        var minPol = Sm.Num.ToKPoly(X.Num);
        if (minPol.Degree != 2.Pow(m))
            throw new();
        
        Console.WriteLine($"S{m} = {minPol}");
    }

    Console.WriteLine();

    var x = FG.QPoly('X');
    for (int m = 1; m <= 6; m++)
    {
        var polys = Primes10000.Take(m).Select(i => x.Pow(2) - i).ToArray();
        var minPol = polys[0];
        foreach (var h in polys.Skip(1))
        {
            var (X0, _) = FG.EPolyXc(minPol, 'y');
            minPol = IntFactorisation.PrimitiveElt(h.Substitute(X0)).r.SubstituteChar('X');
        }

        if (minPol.Degree != 2.Pow(m))
            throw new();
        
        Console.WriteLine($"S{m} = {minPol}");
    }
}
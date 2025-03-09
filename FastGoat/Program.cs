using System.Numerics;
using FastGoat.Commons;
using FastGoat.Examples;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.LWE;
using FastGoat.UserGroup.Padic;
using FastGoat.UserGroup.Polynoms;
using Rq = FastGoat.Structures.VecSpace.KPoly<FastGoat.UserGroup.Integers.Rational>;
using static FastGoat.Commons.IntExt;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");
// RngSeed(259128);
Ring.DisplayPolynomial = MonomDisplay.StarCaret;
// RecomputeAllPrimesUpTo(1000000);

int NthRootsANTV1Composite(int a, int r, int p)
{
    var ai = AmodP(a, p);
    if (r == 1 || ai <= 1)
        return ai;

    return PrimesDecomposition(r).Aggregate(ai, (acc, ri) => NthRootsANTV1Step(acc, ri, p));
}

int NthRootsANTV1Step(int a, int r, int p)
{
    if (Gcd(r, p - 1) == 1)
    {
        var invr = InvModPbez(r, p - 1);
        return PowMod(a, invr, p);
    }

    return NumberTheory.NthRootANTV1(a, r, p);
}

{
    var p = 97;
    for (int r = 2; r < p / 2; r++)
    {
        for (int a = 0; a < p; a++)
        {
            try
            {
                var g = NthRootsANTV1Composite(a, r, p);
                var gr = PowMod(g, r, p);
                Console.WriteLine($"p={p} r={r,3} a={a,3} g={g,3} g^r={gr,3} check:{gr == a}");
            }
            catch (Exception e)
            {
                Console.WriteLine($"p={p} r={r,3} a={a,3} error msg:{e.Message}");
            }
        }

        Console.WriteLine();
    }
}
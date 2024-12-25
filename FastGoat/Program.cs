using System.Numerics;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Lattice;
using Rq = FastGoat.Structures.VecSpace.KPoly<FastGoat.UserGroup.Integers.Rational>;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");
// IntExt.RngSeed(7532159);
Ring.DisplayPolynomial = MonomDisplay.StarCaret;

// ZnBInt-PadicSequence()
IEnumerable<BigInteger> DecompBase(BigInteger a, int b = 2, int size = -1, bool sizeStrict = false)
{
    var a0 = a;
    var k = 0;
    do
    {
        (a0, var r) = BigInteger.DivRem(a0, b);
        yield return r;
    } while (!a0.IsZero || (sizeStrict && ++k < size));
}
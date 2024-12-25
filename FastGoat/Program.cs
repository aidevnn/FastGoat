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
IntExt.RngSeed(7532159);
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
        ++k;
    } while (!a0.IsZero || (sizeStrict && k < size));
}

Vec<Rational> VecBase(int k, int B = 2)
{
    return new(k.SeqLazy().Select(i => new Rational(BigInteger.Pow(B, i))).ToArray());
}

Vec<Rq> DecomposeRq(Rq poly, int degree, int k, int B = 2)
{
    var extr = degree.SeqLazy().Select(i => DecompBase(poly[i].Num, B, k, sizeStrict: true).ToArray()).ToArray();
    // extr.Println(l => $"[{l.Glue(", ")}]", "Decompose);
    return k.SeqLazy().Select(j => degree.SeqLazy().Select(i => new Rational(extr[i][j])).ToKPoly()).ToVec();
}

{
    var (k, B) = (6, 2);
    var vBase = VecBase(k, B);
    var (n, t0, q0) = (16, 8, B.Pow(k));

    var m = FHE.GenUnif(n, q0);
    Console.WriteLine($"m:{m}");
    var gm = DecomposeRq(m, n, k, B);
    vBase.Zip(gm).Println("Decompose");
    var prod = vBase.MulA(gm);
    prod.Println("Decompose");
    Console.WriteLine($"m:{m}");
    Console.WriteLine($" :{prod.Sum()}");
}

// m:45*x^15 + 23*x^14 + 55*x^13 + 20*x^12 + 51*x^11 + 48*x^10 + 11*x^9 + 53*x^8 + 8*x^7 + 12*x^6 + 50*x^5 + 37*x^4 + 41*x^3 + 28*x^2 + 31*x + 16
// Decompose
//     [0, 0, 0, 0, 1, 0]
//     [1, 1, 1, 1, 1, 0]
//     [0, 0, 1, 1, 1, 0]
//     [1, 0, 0, 1, 0, 1]
//     [1, 0, 1, 0, 0, 1]
//     [0, 1, 0, 0, 1, 1]
//     [0, 0, 1, 1, 0, 0]
//     [0, 0, 0, 1, 0, 0]
//     [1, 0, 1, 0, 1, 1]
//     [1, 1, 0, 1, 0, 0]
//     [0, 0, 0, 0, 1, 1]
//     [1, 1, 0, 0, 1, 1]
//     [0, 0, 1, 0, 1, 0]
//     [1, 1, 1, 0, 1, 1]
//     [1, 1, 1, 0, 1, 0]
//     [1, 0, 1, 1, 0, 1]
// Decompose
//     (1, x^15 + x^14 + x^13 + x^11 + x^9 + x^8 + x^4 + x^3 + x)
//     (2, x^14 + x^13 + x^11 + x^9 + x^5 + x)
//     (4, x^15 + x^14 + x^13 + x^12 + x^8 + x^6 + x^4 + x^2 + x)
//     (8, x^15 + x^9 + x^7 + x^6 + x^3 + x^2 + x)
//     (16, x^14 + x^13 + x^12 + x^11 + x^10 + x^8 + x^5 + x^2 + x + 1)
//     (32, x^15 + x^13 + x^11 + x^10 + x^8 + x^5 + x^4 + x^3)
// Decompose
//     x^15 + x^14 + x^13 + x^11 + x^9 + x^8 + x^4 + x^3 + x
//     2*x^14 + 2*x^13 + 2*x^11 + 2*x^9 + 2*x^5 + 2*x
//     4*x^15 + 4*x^14 + 4*x^13 + 4*x^12 + 4*x^8 + 4*x^6 + 4*x^4 + 4*x^2 + 4*x
//     8*x^15 + 8*x^9 + 8*x^7 + 8*x^6 + 8*x^3 + 8*x^2 + 8*x
//     16*x^14 + 16*x^13 + 16*x^12 + 16*x^11 + 16*x^10 + 16*x^8 + 16*x^5 + 16*x^2 + 16*x + 16
//     32*x^15 + 32*x^13 + 32*x^11 + 32*x^10 + 32*x^8 + 32*x^5 + 32*x^4 + 32*x^3
// m:45*x^15 + 23*x^14 + 55*x^13 + 20*x^12 + 51*x^11 + 48*x^10 + 11*x^9 + 53*x^8 + 8*x^7 + 12*x^6 + 50*x^5 + 37*x^4 + 41*x^3 + 28*x^2 + 31*x + 16
//  :45*x^15 + 23*x^14 + 55*x^13 + 20*x^12 + 51*x^11 + 48*x^10 + 11*x^9 + 53*x^8 + 8*x^7 + 12*x^6 + 50*x^5 + 37*x^4 + 41*x^3 + 28*x^2 + 31*x + 16
// 
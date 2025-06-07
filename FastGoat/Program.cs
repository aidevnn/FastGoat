using System.Numerics;
using System.Text;
using FastGoat.Commons;
using FastGoat.Examples;
using FastGoat.Structures;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.EllCurve;
using FastGoat.UserGroup.Integers;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");
Ring.DisplayPolynomial = MonomDisplay.StarCaret;
IntExt.RecomputeAllPrimesUpTo(15000000);

{
    EllipticCurvesPart2.Example10EllApSchoof(); // Time:4.486s
    EllipticCurvesPart2.Example11EllApSchoof(); // Time:3m0s
    
    EC.EllApSchoofBigInt([1, 1], 10.Pow(6) + 3); // Time:3.484s
    EC.EllApSchoofBigInt([1, 1], BigInteger.Pow(10, 10) + 19); // Time:25.543s
    EC.EllApSchoofBigInt([1, 1], BigInteger.Pow(10, 14) + 31);
    // Frob Traces
    //     [2, 0]
    //     [3, 2]
    //     [5, 1]
    //     [7, 0]
    //     [11, 7]
    //     [13, 10]
    //     [17, 4]
    //     [19, 13]
    //     [23, 7]
    // p = 100000000000031 crt = 208037606 ap = -15055264 L = 223092870
    // # EllApSchoof(Ell[1,1](Z/100000000000031Z)) Time:4m1s
    // 
    // 
    // pari
    // ? ellap(ellinit([1,1]),10^14+31)
    // %42 = -15055264
    // 
}
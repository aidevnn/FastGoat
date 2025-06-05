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
IntExt.RecomputeAllPrimesUpTo(200000);

{
    EC.CheckDivOfZero = false;
    EllipticCurvesPart2.Example10EllApSchoof(); // Time:7.632s
    // EC.EllApSchoofBigInt([1, 1], 10.Pow(6) + 3); // slow Time:34.572s
    
    EC.CheckDivOfZero = true;
    EllipticCurvesPart2.Example10EllApSchoof(); // Time:7.081s
    // EllipticCurvesPart2.Example11EllApSchoof(); // Time:5m20s
    
    EC.EllApSchoofBigInt([1, 1], 10.Pow(6) + 3); // Time:4.470s
    // EC.EllApSchoofBigInt([1, 1], BigInteger.Pow(10, 10) + 19); // Time:25.600s

}
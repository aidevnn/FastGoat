using System.Numerics;
using FastGoat.Commons;
using FastGoat.Examples;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.LWE;
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

{
    RecomputeAllPrimesUpTo(500000);
    for (int i = 7; i < 12; i++)
    {
        var n = 1 << i;
        Console.WriteLine($"fg{n}a = {new Regev(n).ExportParams}");
        Console.WriteLine($"fg{n}b = {new RLWE(2 * n).ExportParams}");
        Console.WriteLine();
    }
    // from estimator import *
    // 
    // D = ND.DiscreteGaussian
    // T = ND.SparseTernary
    // U = ND.UniformMod
    // 
    // fg128a = LWEParameters(n=128, q=12547, Xs=U(12547), Xe=D(9.0292), m=+oo, tag='fg128')
    // fg128b = LWEParameters(n=128, q=13313, Xs=T(58, n=128), Xe=D(9.5804), m=+oo, tag='fg128')
    // 
    // fg256a = LWEParameters(n=256, q=32771, Xs=U(32771), Xe=D(12.7673), m=+oo, tag='fg256')
    // fg256b = LWEParameters(n=256, q=36353, Xs=T(120, n=256), Xe=D(14.1628), m=+oo, tag='fg256')
    // 
    // fg512a = LWEParameters(n=512, q=82963, Xs=U(82963), Xe=D(18.0582), m=+oo, tag='fg512')
    // fg512b = LWEParameters(n=512, q=83969, Xs=T(244, n=512), Xe=D(18.2772), m=+oo, tag='fg512')
    // 
    // fg1024a = LWEParameters(n=1024, q=204803, Xs=U(204803), Xe=D(25.5327), m=+oo, tag='fg1024')
    // fg1024b = LWEParameters(n=1024, q=249857, Xs=T(496, n=1024), Xe=D(31.1495), m=+oo, tag='fg1024')
    // 
    // fg2048a = LWEParameters(n=2048, q=495617, Xs=U(495617), Xe=D(36.1082), m=+oo, tag='fg2048')
    // fg2048b = LWEParameters(n=2048, q=520193, Xs=T(1001, n=2048), Xe=D(37.8987), m=+oo, tag='fg2048')
    // 
    // LWE.estimate(fg-xyz)
}
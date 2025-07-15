using System.Numerics;
using System.Text;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.EllCurve;
using FastGoat.UserGroup.Floats;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Polynoms;
using static FastGoat.Commons.IntExt;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");

// 11.1. Use the double-and-add algorithm (XI.1.1) to compute [n]P in E(Fp) for each of the
// following curves and points.
// (a) E : Y 2 = X 3 + 143X + 367, p = 613,
// P = (195, 9),n = 23.
// (b) E : Y 2 = X 3 + 1541X + 1335, p = 3221,
// P = (2898, 439),n = 3211.

{
    var (p, n) = (613, 23);
    var E = EC.EllGroup([143, 367]);
    var Ep = E.ToZnInt(p);
    var P = Ep[195, 9];
    var nP = Ep.Times(P, n);
    var ordP = EC.BSGSlong(Ep, P, Ep.O, p + (int)double.Sqrt(16 * p));
    Console.WriteLine(new { E, p, P, nP, ordP });
}

{
    var (p, n) = (3221, 3211);
    var E = EC.EllGroup([1541, 1335]);
    var Ep = E.ToZnInt(p);
    var P = Ep[2898, 439];
    var nP = Ep.Times(P, n);
    var ordP = EC.BSGSlong(Ep, P, Ep.O, p + (int)double.Sqrt(16 * p));
    Console.WriteLine(new { E, p, P, n, nP, ordP });
}
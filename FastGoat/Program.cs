using System.CodeDom;
using System.Collections;
using System.Diagnostics;
using System.Numerics;
using FastGoat.Commons;
using FastGoat.Examples;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.Naming;
using FastGoat.Structures.Subgroups;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Characters;
using FastGoat.UserGroup.Floats;
using FastGoat.UserGroup.GModuleN;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Matrix;
using FastGoat.UserGroup.Perms;
using FastGoat.UserGroup.Polynoms;
using FastGoat.UserGroup.Words;
using static FastGoat.Commons.IntExt;
using static FastGoat.Commons.EnumerableExt;
using FastGoat.UserGroup.Padic;
using FastGoat.UserGroup.Words.Tools;
using GroupRegX = System.Text.RegularExpressions;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");

KPoly<Rational> RecNewtonInverse(KPoly<Rational> F, int N)
{
    if (F.Degree >= N)
        throw new();

    if (N == 1)
        return F[0].Inv() * F.One;

    var mid = (N / 2) + (N % 2);
    var F0 = F.Div(F.X.Pow(mid)).rem;
    var G = RecNewtonInverse(F0, mid);
    var G0 = (G + (1 - F * G) * G).Div(F.X.Pow(N)).rem;
    return G0;
}

// MolienSum(D8) = 1/(t0^6 - t0^4 - t0^2 + 1)
// MolienSerie(D8) = 3*t^8 + 2*t^6 + 2*t^4 + t^2 + 1
{
    var x = FG.QPoly();
    var P = 1 - x.Pow(2) - x.Pow(4) + x.Pow(6);
    var Pi = RecNewtonInverse(P, 9);
    Console.WriteLine(Pi);
}

// MolienSum(Q8) = (t0^4 - t0^2 + 1)/(t0^6 - t0^4 - t0^2 + 1)
// MolienSerie(Q8) = 3*t^8 + t^6 + 2*t^4 + 1
{
    var x = FG.QPoly();
    var P = 1 - x.Pow(2) - x.Pow(4) + x.Pow(6);
    var Q = 1 - x.Pow(2) + x.Pow(4);
    var Pi = RecNewtonInverse(P, 9);
    Console.WriteLine((Q * Pi).Div(x.Pow(9)).rem);
}

// MolienSum(D12) = 1/(t0^8 - t0^6 - t0^2 + 1)
// MolienSerie(D12) = 3*t^12 + 2*t^10 + 2*t^8 + 2*t^6 + t^4 + t^2 + 1
{
    var x = FG.QPoly();
    var P = 1 - x.Pow(2) - x.Pow(6) + x.Pow(8);
    var Pi = RecNewtonInverse(P, 13);
    Console.WriteLine(Pi);
}

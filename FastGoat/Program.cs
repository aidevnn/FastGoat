using System.CodeDom;
using System.Collections;
using System.Diagnostics;
using System.Numerics;
using System.Reflection.Emit;
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
using FastGoat.UserGroup.DatabaseSmallGroups;
using FastGoat.UserGroup.Floats;
using FastGoat.UserGroup.GModuleN;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Matrix;
using FastGoat.UserGroup.GModuleN;
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

// wikipedia algorithme de chudnovski
// def binary_split(a, b):
//     if b == a + 1:
//         Pab = -(6*a - 5)*(2*a - 1)*(6*a - 1)
//         Qab = 10939058860032000 * a**3
//         Rab = Pab * (545140134*a + 13591409)
//     else:
//         m = (a + b) // 2
//         Pam, Qam, Ram = binary_split(a, m)
//         Pmb, Qmb, Rmb = binary_split(m, b)
//         
//         Pab = Pam * Pmb
//         Qab = Qam * Qmb
//         Rab = Qmb * Ram + Pam * Rmb
//     return Pab, Qab, Rab
// 
// 
// def chudnovsky(n):
//     """Chudnovsky algorithm."""
//     P1n, Q1n, R1n = binary_split(1, n)
//     return (426880 * decimal.Decimal(10005).sqrt() * Q1n) / (13591409*Q1n + R1n)
// 
// 
// print(f"1 = {chudnovsky(2)}")  # 3.141592653589793238462643384
// 
// decimal.getcontext().prec = 100 # number of digits of decimal precision
// for n in range(2,10):
//     print(f"{n} = {chudnovsky(n)}")  # 3.14159265358979323846264338...

(BigReal Pab, BigReal Qab, BigReal Rab) BinarySplit(int a, int b, int O = 20)
{
    if (b == a + 1)
    {
        var ba = BigReal.BrOne(O) * a;
        var Pab = -(6 * ba - 5) * (2 * ba - 1) * (6 * ba - 1);
        var Qab = BigReal.FromBigInteger(new(10939058860032000), O) * ba.Pow(3);
        var Rab = Pab * (545140134 * ba + 13591409);
        return (Pab, Qab, Rab);
    }
    else
    {
        var m = (a + b) / 2;
        var (Pam, Qam, Ram) = BinarySplit(a, m, O);
        var (Pmb, Qmb, Rmb) = BinarySplit(m, b, O);

        var Pab = Pam * Pmb;
        var Qab = Qam * Qmb;
        var Rab = Qmb * Ram + Pam * Rmb;
        return (Pab, Qab, Rab);
    }
}

BigReal Chudnovsky(int n, int O = 20)
{
    var (P1n, Q1n, R1n) = BinarySplit(1, n, O + 3);
    var pi = (426880 * BigReal.FromBigInteger(10005, O + 3).Sqrt() * Q1n) / (13591409 * Q1n + R1n);
    return BigReal.Round(pi, O - 1).ToBigReal(O);
}

BigReal EulerNumber(int O)
{
    var fact = BigReal.BrOne(O + 3);
    var facti = fact;
    var n = fact.One;
    var sum = fact.One;
    while (!facti.IsZero())
    {
        sum += facti;
        n += 1;
        fact *= n;
        facti = fact.Inv();
    }

    return BigReal.Round(sum, O - 1).ToBigReal(O);
}

{
    foreach (var n in new[] { 2, 10, 20, 80 })
    {
        var O = 25 * n / 2;
        var pi0 = Chudnovsky(n, O);
        Console.WriteLine(pi0.ToFixForm());
        Console.WriteLine((pi0 - BigReal.Pi(O)).ToFixForm());
    }

    Console.WriteLine();

    foreach (var O in new[] { 25, 125, 250, 1000 })
    {
        var e0 = EulerNumber(O);
        var e1 = BigReal.E(O);
        Console.WriteLine(e0.ToFixForm());
        Console.WriteLine((e0 - e1).ToFixForm());
    }
}
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

T Simpson<T>(F<T> f, T a, T b, int n) where T : struct, IElt<T>, IRingElt<T>, IFieldElt<T>, IFloatElt<T>
{
    var h = (b - a) / n;
    var sum = f(a) + f(b);
    var xi = a + h;
    for (int i = 1; i < n; i++)
    {
        var s = i % 2 == 0 ? 2 : 4;
        sum += s * f(xi);
        xi += h;
    }
    
    return sum * h / 3;
}

T SimpsonThreeEighths<T>(F<T> f, T a, T b, int n) where T : struct, IElt<T>, IRingElt<T>, IFieldElt<T>, IFloatElt<T>
{
    n = ((n + 2) / 3) * 3;
    var h = (b - a) / n;
    var sum = f(a) + f(b);
    var xi = a + h;
    for (int i = 1; i < n; i++)
    {
        var s = (i % 3 != 0 ? 3 : 2);
        sum +=  s * f(xi);
        xi += h;
    }
    
    return sum * (3 * h / 8);
}

KPoly<K> Primitive<K>(KPoly<K> P) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
{
    var z0 = P.KZero;
    return new(P.x, z0, P.Coefs.Select((c, i) => c / (i + 1)).Prepend(z0).ToArray());
}

{
    {
        var n = 50;
        var X = FG.BRealPoly(40, 'X');
        var P0 = X.Pow(2) + 3 * X - 5;
        var P1 = Primitive(P0);
        
        var k = Rng.Next(-10, 10);
        var r = Rng.Next(k + 1, k + 10);
        var (a, b) = (P0.KOne * k, P0.KOne * r);
        Console.WriteLine(Simpson(x => P0.Substitute(x), a, b, n).ToFixForm());
        Console.WriteLine(SimpsonThreeEighths(x => P0.Substitute(x), a, b, n).ToFixForm());
        Console.WriteLine((P1.Substitute(b) - P1.Substitute(a)).ToFixForm());
        Console.WriteLine();
    }

    {
        F<Dble> cos = x => double.Cos(x);
        F<Dcml> f = x => 1 / (1 + x * x);

        Dble pi = double.Pi;
        var n = 800;
        var a = pi / 4;
        Console.WriteLine(Simpson(cos, 0, a, n));
        Console.WriteLine(SimpsonThreeEighths(cos, 0, a, n));
        Console.WriteLine(double.Sin(a));
        
        Console.WriteLine(SimpsonThreeEighths(f, 0, 1.0m, n));
        Console.WriteLine(double.Atan(1.0));
        Console.WriteLine(a);
        Console.WriteLine((BigReal.Pi(30) / 4).ToFixForm());
    }
}

delegate T F<T>(T a);
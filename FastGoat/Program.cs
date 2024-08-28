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
        sum += s * f(xi);
        xi += h;
    }

    return sum * (3 * h / 8);
}

KPoly<K> Primitive<K>(KPoly<K> P) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
{
    var z0 = P.KZero;
    return new(P.x, z0, P.Coefs.Select((c, i) => c / (i + 1)).Prepend(z0).ToArray());
}

void run1()
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
        Console.WriteLine();
        
        Console.WriteLine(Simpson(f, 0, 0.15m, 2 * n));
        Console.WriteLine(SimpsonThreeEighths(f, 0, 0.15m, 2 * n));
        Console.WriteLine(ATan<Dcml>(0.15m));
        Console.WriteLine(ATanOld(BigReal.FromBigInteger(15, 124) / 100).ToFixForm()); // errors accumulations
        Console.WriteLine(ATan(BigReal.FromBigInteger(15, 124) / 100).ToFixForm());
        Console.WriteLine(double.Atan(0.15));
    }
}

void CordicDouble(double theta)
{
    var (xi, yi, zi) = (1.0, 0.0, theta);
    var ti = 1.0;
    var i = 0;
    do
    {
        ++i;
        var Ki = 1.0 / double.Sqrt(1 + ti * ti);
        var di = zi > 0 ? 1 : -1;
        (xi, yi) = (Ki * (xi - di * yi * ti), Ki * (yi + di * xi * ti));
        zi = zi - di * double.Atan(ti);
        ti *= 0.5;
    } while (double.Abs(zi) > 1e-17);

    var (cos, sin) = (double.Cos(theta), double.Sin(theta));
    Console.WriteLine(new { i, xi, yi, zi });
    Console.WriteLine(new { cos, sin });
    Console.WriteLine();
}

T ATan<T>(T x) where T : struct, IElt<T>, IRingElt<T>, IFieldElt<T>, IFloatElt<T>
{
    if (x.Equals(x.One))
    {
        if (x is BigReal x0)
            return T.Pi(x0.O) / 4;

        return T.Pi() / 4;
    }

    if (x.CompareTo(x.One) == 1)
    {
        var pi2 = x is BigReal x0 ? T.Pi(x0.O) / 2 : T.Pi() / 2;
        return pi2 - ATan(x.Inv());
    }

    if (x.CompareTo(x.Zero) == -1)
        return -ATan(-x);

    T sum = x.Zero;
    var x2 = -x * x;
    var idx = 1;
    var x2k1 = x;
    while (idx < 10000)
    {
        var c = x2k1 * (idx * x2k1.One).Inv();
        if (c.IsZero())
            break;

        sum += c;
        x2k1 *= x2;
        idx += 2;
    }

    return sum;
}

T ATanOld<T>(T x) where T : struct, IElt<T>, IRingElt<T>, IFieldElt<T>, IFloatElt<T>
{
    T sum = x.Zero;
    var x2 = -x * x;
    var idx = 1;
    var x2k1 = x;
    while (idx < 1000)
    {
        var c = x2k1 / idx;
        if (c.IsZero())
            break;

        sum += c;
        x2k1 *= x2;
        idx += 2;
    }

    return sum;
}

(T xi, T yi) Cordic<T>(T theta) where T : struct, IElt<T>, IRingElt<T>, IFieldElt<T>, IFloatElt<T>
{
    var (xi, yi, zi) = (theta.One, theta.Zero, theta);
    var half = theta.One / 2;
    var ti = theta.One;
    var i = 0;
    do
    {
        var Ki = T.Sqrt(1 + ti * ti).Inv();
        var di = zi.CompareTo(zi.Zero) != -1 ? 1 : -1;
        (xi, yi) = (Ki * (xi - di * yi * ti), Ki * (yi + di * xi * ti));
        zi = zi - di * ATan(ti);
        ti *= half;
        ++i;
    } while (!zi.IsZero());

    var (cos, sin) = (double.Cos(theta.ToDouble), double.Sin(theta.ToDouble));
    Console.WriteLine(new { i, xi, yi, zi });
    Console.WriteLine(new { cos, sin });
    Console.WriteLine();

    return (xi, yi);
}

void run2()
{
    Console.WriteLine(ATan(Dcml.DcmlOne() / 2));
    Console.WriteLine(double.Atan(0.5));
    Console.WriteLine(ATan(Dcml.DcmlOne() / 4));
    Console.WriteLine(double.Atan(0.25));

    {
        var pi = double.Pi;
        CordicDouble(pi / 6);
        CordicDouble(pi / 5);
    }

    {
        var pi = Dble.Pi();
        Cordic(pi / 6);
        Cordic(pi / 5);
    }

    {
        var pi = Dcml.Pi();
        Cordic(pi / 6);
        Cordic(pi / 5);
    }

    {
        var pi = BigReal.Pi(110);
        var (cos, sin) = Cordic(pi / 6);
        Console.WriteLine($"{{ cos = {cos.ToFixForm()}, sin = {sin.ToFixForm()} }}");
        (cos, sin) = Cordic(pi / 5);
        Console.WriteLine($"{{ cos = {cos.ToFixForm()}, sin = {sin.ToFixForm()} }}");
    }
}

{
    run1();
    run2();
}

delegate T F<T>(T a);
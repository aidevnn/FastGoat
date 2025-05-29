using System.Numerics;
using FastGoat.Commons;
using FastGoat.Examples;
using FastGoat.Structures;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.EllCurve;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Polynoms;
using static FastGoat.Commons.IntExt;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");
Ring.DisplayPolynomial = MonomDisplay.StarCaret;

EllPt<FracPoly<FracPoly<T>>> 
    Op<T>(EllGroup<T> E, EllPt<FracPoly<FracPoly<T>>> e1, EllPt<FracPoly<FracPoly<T>>> e2, KPoly<FracPoly<T>> s,
        KPoly<KPoly<T>> R0, KPoly<KPoly<T>> R1) where T : struct, IElt<T>, IRingElt<T>, IFieldElt<T>
{
    if (e1.IsO)
        return e2;

    if (e2.IsO)
        return e1;

    var (a1, a2, a3, a4, a6) = E.ArrCoefs.Select(c => c * e1.X.KOne).Deconstruct();
    var ((x1, y1), (x2, y2)) = (e1, e2);
    if (!x1.Equals(x2))
    {
        var alpha = EC.Simplify((y2 - y1) / (x2 - x1), s, R0, R1);
        var x3 = EC.Simplify(alpha.Pow(2) + a1 * alpha - a2 - x2 - x1, s, R0, R1);
        var y3 = EC.Simplify(x3 * alpha - x1 * alpha + y1, s, R0, R1);
        return new(x3, EC.Simplify(-a1 * x3 - a3 - y3, s, R0, R1));
    }
    else
    {
        if (!y1.Equals(-a1 * x2 - a3 - y2))
        {
            var alpha = EC.Simplify((3 * x1.Pow(2) - y1 * a1 + 2 * x1 * a2 + a4) / (x1 * a1 + a3 + 2 * y1), s, R0, R1);
            var x3 = EC.Simplify(alpha.Pow(2) + a1 * alpha - a2 - 2 * x1, s, R0, R1);
            var y3 = EC.Simplify(alpha * (x3 - x1) + y1, s, R0, R1);
            return new(x3, EC.Simplify(-a1 * x3 - a3 - y3, s, R0, R1));
        }
        else
            return new();
    }
}

{
    var (p, m) = (37, 10);
    var E0 = EC.EllGroup([1, 0, 1, -4, 4]).ToZnInt(p);
    Console.WriteLine(E0);
    Console.WriteLine(E0.EqStr);

    var (R0, R1, psi, fdiv) = EC.DivisionPolynomial(E0, m);
    var Y = R0.ToFrac().X;
    var X = Y.KOne.X * Y.One;
    fdiv.Println("fdiv");
    psi.Println("psi");

    var Ob = new EllPt<FracPoly<FracPoly<ZnInt>>>();
    var ptb = new EllPt<FracPoly<FracPoly<ZnInt>>>(X, Y);
    var lista = new[] { Ob, ptb }.Index().ToDictionary(e => e.Index, e => e.Item);
    var listb = new[] { Ob }.Index().ToDictionary(e => e.Index, e => e.Item);
    Console.WriteLine();

    var psi2 = psi[2].ToFrac().Num;
    psi[-1] = -psi[1];
    
    GlobalStopWatch.AddLap();
    for (int n = 1; n <= m - 2; n++)
    {
        Console.WriteLine(new { n });
        var npta = lista[n] = Op(E0, lista[n / 2], lista[n - n / 2], psi2, R0, R1);
        Console.WriteLine(new { npta });
        Console.WriteLine();
    }
    
    GlobalStopWatch.Show();

    GlobalStopWatch.AddLap();
    for (int n = 1; n <= m - 2; n++)
    {
        Console.WriteLine(new { n });
        var nptb = listb[n] = EC.NPt(n, psi, R0, R1).nP;
        Console.WriteLine(new { nptb });
    }
    
    GlobalStopWatch.Show();
    Console.WriteLine();

    Console.WriteLine(lista.All(e => e.Value.Equals(listb[e.Key])));
}
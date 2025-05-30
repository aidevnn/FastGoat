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

EllPt<EllPoly<T>> Op2<T>(EllGroup<T> E, EllPt<EllPoly<T>> e1, EllPt<EllPoly<T>> e2)
    where T : struct, IElt<T>, IRingElt<T>, IFieldElt<T>
{
    if (e1.IsO)
        return e2;

    if (e2.IsO)
        return e1;

    var (a1, a2, a3, a4, a6) = E.Coefs;
    var ((x1, y1), (x2, y2)) = (e1, e2);
    if (!x1.Equals(x2))
    {
        var alpha = (y2 - y1) / (x2 - x1);
        var x3 = alpha.Pow(2) + a1 * alpha - a2 - x2 - x1;
        var y3 = -(alpha * (x3 - x1) + y1) - a1 * x3 - a3;
        return new(x3, y3);
    }
    else
    {
        if (!y1.Equals(-a1 * x2 - a3 - y2))
        {
            var alpha = (3 * x1.Pow(2) - y1 * a1 + 2 * x1 * a2 + a4) / (x1 * a1 + a3 + 2 * y1);
            var x3 = alpha.Pow(2) + a1 * alpha - a2 - 2 * x1;
            var y3 = -(alpha * (x3 - x1) + y1) - a1 * x3 - a3;
            return new(x3, y3);
        }
        else
            return new();
    }
}

void testNPT()
{
    var (p, m) = (37, 10);
    var E0 = EC.EllGroup([1, 0, 1, -4, 4]).ToZnInt(p);
    Console.WriteLine(E0);
    Console.WriteLine(E0.EqStr);

    var (R0, R1, psi, fdiv) = EC.DivisionPolynomial(E0, m);
    fdiv.Println("fdiv");
    psi.Println("psi");

    var Y0 = R0.ToFrac().X;
    var X0 = Y0.KOne.X * Y0.One;
    var ptb = new EllPt<FracPoly<FracPoly<ZnInt>>>(X0, Y0);
    var Ob = new EllPt<FracPoly<FracPoly<ZnInt>>>();
    var lista = new[] { Ob, ptb }.Index().ToDictionary(e => e.Index, e => e.Item);
    Console.WriteLine();

    psi[-1] = -psi[1];

    GlobalStopWatch.AddLap();
    var psi2 = psi[2].ToFrac().Num;
    for (int n = 1; n <= m - 2; n++)
    {
        Console.WriteLine(new { n });
        var npta = lista[n] = Op(E0, lista[n / 2], lista[n - n / 2], psi2, R0, R1);
        Console.WriteLine(new { npta });
    }

    GlobalStopWatch.Show();
    Console.WriteLine();

    var (_, _, eq, sd) = E0.GetPolynomials();
    var (X2, Y2) = EllPoly<ZnInt>.GetXY(eq, sd, eq.Zero);
    var (y0, x0) = eq.Indeterminates.Deconstruct();
    var listb = new[] { new EllPt<EllPoly<ZnInt>>(), new(X2, Y2) }.Index().ToDictionary(e => e.Index, e => e.Item);
    GlobalStopWatch.AddLap();
    for (int n = 1; n <= m - 2; n++)
    {
        Console.WriteLine(new { n });
        var nptc = listb[n] = Op2(E0, listb[n / 2], listb[n - n / 2]);
        Console.WriteLine(new { nptc });
        Console.WriteLine();
    }

    GlobalStopWatch.Show();
    Console.WriteLine();
    
    var listc = new[] { Ob }.Index().ToDictionary(e => e.Index, e => e.Item);
    GlobalStopWatch.AddLap();
    for (int n = 1; n <= m - 2; n++)
    {
        Console.WriteLine(new { n });
        var nptc = listc[n] = EC.NPt(n, psi, R0, R1).nP;
        Console.WriteLine(new { nptc });
    }

    GlobalStopWatch.Show();
    Console.WriteLine();

    Console.WriteLine(lista.All(e => e.Value.Equals(listb[e.Key].ToFracPoly(x0, y0))));
    Console.WriteLine(lista.All(e => e.Value.Equals(listc[e.Key])));
    Console.WriteLine();
}

{
    GlobalStopWatch.Restart();

    testNPT();
    testNPT();
}
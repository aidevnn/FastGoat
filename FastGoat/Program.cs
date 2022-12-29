using System.Globalization;
using FastGoat.Commons;
using FastGoat.Examples;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Matrix;
using FastGoat.UserGroup.Perms;
using FastGoat.UserGroup.Polynoms;
using FastGoat.UserGroup.Words;
using FastGoat.UserGroup.Words.ToddCoxeter;
using static FastGoat.Commons.IntExt;
using static FastGoat.Commons.EnumerableExt;
using System.Numerics;
using System.Runtime.InteropServices;
using System.Xml;
using FastGoat.UserGroup.Padic;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");

Polynomial<Rational, Xi> Simplify(Polynomial<Rational, Xi> f)
{
    var gcdNum = IntExt.GcdBigInt(f.Coefs.Where(c => !c.Value.IsZero()).Select(c => c.Value.Num).ToArray());
    var gcdDenom = IntExt.GcdBigInt(f.Coefs.Where(c => !c.Value.IsZero()).Select(c => c.Value.Denom).ToArray());
    var r = new Rational(gcdDenom, gcdNum);
    return f * r;
}

//
// {
//     var (X, Y, Z, T) = Ring.Polynomial("X", "Y", "Z", "T", Rational.KZero());
//     X.Indeterminates.SetOrder(MonomOrder.Lex);
//     var g1 = 2 * X * Y - X * Z;
//     var g2 = X.Pow(2);
//     var g3 = X.Pow(2) * Z;
//     Console.WriteLine(new { g1, g2, g3 });
//     Console.WriteLine("X * g1 - 2 * Y * g = -g3 , X * Z * g1 - 2 * Y * g3 = -Z * g3 , Z * g2 - g3 = 0");
//     Console.WriteLine($"{X * g1 - 2 * Y * g2} = {-g3} , {X * Z * g1 - 2 * Y * g3} = {-Z * g3} , {Z * g2 - g3}");
//
//     SyzygieLcm(g1, g2);
//     SyzygieLcm(g1, g3);
//     SyzygieLcm(g2, g3);
// }

{
    Monom.Display = MonomDisplay.StarCaret;
    var xis = Ring.Polynomial(Rational.KZero(), "x", "y", "z");
    var x = xis[0];
    x.Indeterminates.SetOrder(MonomOrder.Lex);
    Console.WriteLine(x.Indeterminates);
    var all = 3.Range().Grid3D().Select(e => xis[0].Pow(e.t1) * xis[1].Pow(e.t2) * xis[2].Pow(e.t3)).Ascending().ToArray();
    Console.WriteLine(all.Glue(", "));
}

{
    Monom.Display = MonomDisplay.StarCaret;
    var xis = Ring.Polynomial(Rational.KZero(), "x", "y", "z");
    var x = xis[0];
    x.Indeterminates.SetOrder(MonomOrder.GrLex);
    Console.WriteLine(x.Indeterminates);
    var all = 3.Range().Grid3D().Select(e => xis[0].Pow(e.t1) * xis[1].Pow(e.t2) * xis[2].Pow(e.t3)).Ascending().ToArray();
    Console.WriteLine(all.Glue(", "));
}

{
    Monom.Display = MonomDisplay.StarCaret;
    var xis = Ring.Polynomial(Rational.KZero(), "x", "y", "z");
    var x = xis[0];
    x.Indeterminates.SetOrder(MonomOrder.GrevLex);
    Console.WriteLine(x.Indeterminates);
    var all = 3.Range().Grid3D().Select(e => xis[0].Pow(e.t1) * xis[1].Pow(e.t2) * xis[2].Pow(e.t3)).Ascending().ToArray();
    Console.WriteLine(all.Glue(", "));
}

{
    var (X, Y) = Ring.Polynomial("X", "Y", Rational.KZero());
    var f = X.Pow(4) + 3 * X - 7 * X * Y.Pow(2) - 3;
    var g = Y.Pow(3) + X * Y - Y + 2;
    var h = f * g + 4 * X * Y.Pow(2) - X * X * Y + 5;
    Console.WriteLine(h);
    var (x1, x2) = (X.Indeterminates[0], X.Indeterminates[1]);
    Console.WriteLine(X.Indeterminates);
    Console.WriteLine(h.CoefMax(x1));
    Console.WriteLine(h.CoefMax(x2));
    Console.WriteLine(f.LeadingDetails);
    Console.WriteLine(g.LeadingDetails);
    Console.WriteLine(h.LeadingDetails);
    Console.WriteLine();

    Console.WriteLine(h.Div(f));
    Console.WriteLine(h.Div(g));
}

{
    var (X, Y, Z, T) = Ring.Polynomial("X", "Y", "Z", "T", Rational.KZero());
    X.Indeterminates.SetOrder(MonomOrder.Lex);
    // X + Y - Z, X2 - 2T2 , Y2 - 5T2
    var (p1, p2, p3) = (X + Y - Z, X.Pow(2) - 2 * T.Pow(2), Y.Pow(2) - 5 * T.Pow(2));
    var bs = Ring.GroebnerBase(p1, p2, p3);
    Console.WriteLine(bs.Select(Simplify).Glue("; "));
    Console.WriteLine();
}

{
    var (X, Y, Z, T) = Ring.Polynomial("X", "Y", "Z", "T", Rational.KZero());
    X.Indeterminates.SetOrder(MonomOrder.Lex);
    var (p1, p2, p3) = (X + Y - Z, X.Pow(2) - 2 * T.Pow(2), Y.Pow(2) - 6 * T.Pow(2));
    var bs = Ring.GroebnerBase(p1, p2, p3);
    Console.WriteLine(bs.Select(Simplify).Glue("; "));
    Console.WriteLine();
}

{
    var (X, Y, Z, T) = Ring.Polynomial("X", "Y", "Z", "T", Rational.KZero());
    X.Indeterminates.SetOrder(MonomOrder.Lex);
    var (p1, p2, p3) = (X + Y - Z, X.Pow(2) + T.Pow(2), Y.Pow(2) - 2 * T.Pow(2));
    var bs = Ring.GroebnerBase(p1, p2, p3);
    Console.WriteLine(bs.Select(Simplify).Glue("; "));
    Console.WriteLine();
}

{
    var (X, Y, Z, T) = Ring.Polynomial("X", "Y", "Z", "T", Rational.KZero());
    X.Indeterminates.SetOrder(MonomOrder.Lex);
    var (p1, p2, p3) = (X + Y - Z, X.Pow(2) + 1, Y.Pow(2) - 2);
    var bs = Ring.GroebnerBase(p1, p2, p3);
    Console.WriteLine(bs.Select(Simplify).Glue("; "));
    Console.WriteLine();
}

{
    var (X, Y, Z, T) = Ring.Polynomial("X", "Y", "Z", "T", Rational.KZero());
    X.Indeterminates.SetOrder(MonomOrder.Lex);
    var (p1, p2, p3) = (X + Y - Z, X.Pow(4) + X.Pow(3) + X.Pow(2) + X + 1, Y.Pow(2) - 5);
    var bs = Ring.GroebnerBase(p1, p2, p3);
    Console.WriteLine(bs.Select(Simplify).Glue("\n"));
    Console.WriteLine();
}

{
    Monom.Display = MonomDisplay.StarCaret;
    var (X, Y, Z, T) = Ring.Polynomial("x", "y", "z", "t", Rational.KZero());
    X.Indeterminates.SetOrder(MonomOrder.Lex);
    var (p1, p2, p3) = (X + Y - Z, X.Pow(4) + X.Pow(3) + X.Pow(2) + X + 1, Y.Pow(2) - 5);
    var bs = Ring.MinimalGroebnerBase(p1, p2, p3);
    Console.WriteLine("[{0}]", bs.Order().Select(Simplify).Glue(", "));
    Console.WriteLine("[{0}]", bs.OrderDescending().Select(Simplify).Glue(", "));
    Console.WriteLine();
}

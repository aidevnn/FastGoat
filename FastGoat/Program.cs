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
 
Monom.Display = MonomDisplay.StarCaret;
{
    var (x, y) = Ring.Polynomial("x", "y", Rational.KZero());
    x.Indeterminates.SetOrder(MonomOrder.GrLex);
    var (p1, p2) = (x.Pow(3) - 2 * x * y, x.Pow(2) * y - 2 * y.Pow(2) + x);
    Console.WriteLine($"[{p1}, {p2}]");
    var bs = Ring.GroebnerBasis(p1, p2);
    Console.WriteLine(bs.Select(Simplify).Glue("\n"));
    Console.WriteLine();
}

{
    var (x, y) = Ring.Polynomial("x", "y", Rational.KZero());
    x.Indeterminates.SetOrder(MonomOrder.GrLex);
    var (p1, p2) = (x.Pow(3) - 2 * x * y, x.Pow(2) * y - 2 * y.Pow(2) + x);
    Console.WriteLine($"[{p1}, {p2}]");
    var bs = Ring.ReducedGroebnerBasis(p1, p2);
    Console.WriteLine(bs.Select(Simplify).Glue("\n"));
    Console.WriteLine();
}

{
    var (x, y) = Ring.Polynomial("x", "y", Rational.KZero());
    x.Indeterminates.SetOrder(MonomOrder.GrLex);
    var (p1, p2) = (x.Pow(3) - 2 * x * y, x.Pow(2) * y - 2 * y.Pow(2) + x);
    Console.WriteLine($"[{p1}, {p2}]");
    var bs = Ring.GroebnerBasis(p1, p2);
    Console.WriteLine(bs.Select(Simplify).Glue("\n"));
    Console.WriteLine();
}

{
    var (x, y, z, t) = Ring.Polynomial("x", "y", "z", "t", Rational.KZero());
    x.Indeterminates.SetOrder(MonomOrder.Lex);
    var (p1, p2, p3) = (x + y - z, x * x - 2 * t * t, y * y - 5 * t * t);
    Console.WriteLine($"[{p1}, {p2}, {p3}]");
    var bs = Ring.GroebnerBasis(p1, p2, p3);
    Console.WriteLine(bs.Select(Simplify).Glue("\n"));
    Console.WriteLine();
}

{
    var (x, y, z, t) = Ring.Polynomial("x", "y", "z", "t", Rational.KZero());
    x.Indeterminates.SetOrder(MonomOrder.Lex);
    var (p1, p2, p3) = (x + y - z, x * x - 2 * t * t, y * y - 5 * t * t);
    Console.WriteLine($"[{p1}, {p2}, {p3}]");
    var bs = Ring.ReducedGroebnerBasis(p1, p2, p3);
    Console.WriteLine(bs.Select(Simplify).Glue("\n"));
    Console.WriteLine();
}

{
    var (x, y, z, t) = Ring.Polynomial("x", "y", "z", "t", Rational.KZero());
    x.Indeterminates.SetOrder(MonomOrder.GrLex);
    var (p1, p2, p3) = (x + y - z, x * x - 2 * t * t, y * y - 5 * t * t);
    Console.WriteLine($"[{p1}, {p2}, {p3}]");
    var bs = Ring.GroebnerBasis(p1, p2, p3);
    Console.WriteLine(bs.Select(Simplify).Glue("\n"));
    Console.WriteLine();
}

{
    var (x, y, z, t) = Ring.Polynomial("x", "y", "z", "t", Rational.KZero());
    x.Indeterminates.SetOrder(MonomOrder.GrLex);
    var (p1, p2, p3) = (x + y - z, x * x - 2 * t * t, y * y - 5 * t * t);
    Console.WriteLine($"[{p1}, {p2}, {p3}]");
    var bs = Ring.ReducedGroebnerBasis(p1, p2, p3);
    Console.WriteLine(bs.Select(Simplify).Glue("\n"));
    Console.WriteLine();
}

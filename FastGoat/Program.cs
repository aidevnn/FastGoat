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

//
// {
//     var cyclos = CyclotomicsPolys(12);
//     foreach(var p in Primes10000.Take(7))
//     {
//         var x = FG.QPoly('X');
//         var g = (x.Pow(p) - 1) / (x - 1);
//         var c = FG.EPoly(g, 'c');
//         var (Tr, Norm) = AlgebraicFactorization.TraceNorm(c - 1);
//         Console.WriteLine(new { p , g, Tr, Norm });
//     }
// }

{
    var n = 40;
    var q = 9;
    
    DisplayGroup.HeadElements(new NthRootQ(n));
    DisplayGroup.HeadElements(new NthRootFq(n, q));
}

//
// {
//     var a = FG.EQPoly('a', 1, 1, 1);
//     var x = FG.KPoly('x', a);
//     var f = x.Pow(6) - 5;
//     AlgebraicFactorization.AlgebraicFactors(f);
// }

// {
//     // f = X.Pow(36) + 1064340*X.Pow(30) + 352940951214*X.Pow(24) + 268377098275067940*X.Pow(18)
//     // + -12793029463642800406311*X.Pow(12) + 1708184999994920743701769920*X.Pow(6) + 30649678085347883679332069277696
//     var str = "1064340,352940951214,268377098275067940,-12793029463642800406311,1708184999994920743701769920,30649678085347883679332069277696";
//     var arr = str.Split(',').Select(BigInteger.Parse).ToArray();
//     var x = FG.QPoly('X');
//     var g = x.Pow(36) + arr.Reverse().Select((a, i) => x.Pow(i * 6) * new Rational(a))
//         .Aggregate(x.Zero, (acc, a) => a + acc);
//     Console.WriteLine(g);
//     var g1 = g.Coefs.Select((c, i) => c * x.Pow(i / 6)).Aggregate(x.Zero, (acc, a) => a + acc);
//     Console.WriteLine(g1);
//     var g2 = g1.Substitute(x * 243) / new Rational(BigInteger.Pow(243, 6));
//     Console.WriteLine(g2);
//     PolynomialFactorizationPart2.FirrQ(g2);
// }
//
// {
//     var x = Ring.Polynomial(Rational.KZero());
//     var ords = new[] { MonomOrder.Lex, MonomOrder.GrLex, MonomOrder.RevLex, MonomOrder.GrevLex };
//     
//     foreach (var o in ords)
//     {
//         x.Indeterminates.SetOrder(o);
//         Console.WriteLine(x.Indeterminates);
//         var X = x.Indeterminates[0];
//         var f = 2 * x.Pow(3) + 3 * x + 1;
//         var g = 7 * x.Pow(2) + x + 3;
//         var h = f * g;
//         Console.WriteLine(new { f, g, h });
//         Console.WriteLine(new { ldf = f.LeadingDetails, ldg = g.LeadingDetails, ldh = h.LeadingDetails });
//         Console.WriteLine(h.Div(g));
//         Console.WriteLine(h.Div(f));
//         Console.WriteLine(f.D(X));
//         Console.WriteLine(g.D(X));
//         Console.WriteLine(f.Div(x.Pow(2)));
//         Console.WriteLine();
//     }
// }

//
// {
//     Monom.Display = MonomDisplay.DotSuperscript;
//     var (x1, x2, x3, x4) = Ring.Polynomial("X1", "X2", "X3", "X4", ZnInt.KZero());
//     var lt = new[] { x1, x2 * x4, x3 };
//     
//     x1.Indeterminates.SetOrder(MonomOrder.Lex);
//     Console.WriteLine(lt.Order().Glue(" < ", "{0,5}"));
//     x1.Indeterminates.SetOrder(MonomOrder.RevLex);
//     Console.WriteLine(lt.Order().Glue(" < ", "{0,5}"));
//     Console.WriteLine();
// }
//
// {
//     Monom.Display = MonomDisplay.DotSuperscript;
//     var n = 3;
//     var xs = Ring.Polynomial(ZnInt.KZero(), n.Range(1).Select(i => $"X{i}").ToArray());
//     var lt = Partitions32[n].Select(l => l.Concat(Enumerable.Repeat(0, n - l.Count)).ToArray())
//         .SelectMany(l => GetPermutations(n).Select(p => p.Select(i => l[i - 1]).ToArray()))
//         .Distinct(new SequenceEquality<int>())
//         .Select(l => l.Select((c, i) => xs[i].Pow(c)).Aggregate((a, b) => a * b)).ToArray();
//     
//     xs[0].Indeterminates.SetOrder(MonomOrder.GrLex);
//     Console.WriteLine(lt.Order().Glue(" < ", "{0,8}"));
//     xs[0].Indeterminates.SetOrder(MonomOrder.GrevLex);
//     Console.WriteLine(lt.Order().Glue(" < ", "{0,8}"));
//     Console.WriteLine();
// }
//
// {
//     var (x, y) = Ring.Polynomial("X", "Y", Rational.KZero());
//     var f = y * x.Pow(2) - x;
//     var g1 = y * x - 1;
//     var g2 = x.Pow(2);
//     Console.WriteLine(new { f, g1, g2 });
//     Console.WriteLine(f.Div(g1));
//     Console.WriteLine(f.Div(g2));
//
//     var h = (x + 2 * y) * g1 + (y.Pow(2) - x) * g2;
//     Console.WriteLine(new { h });
//     Console.WriteLine(h.Div(g1));
//     Console.WriteLine(h.Div(g1).rem.Div(g2));
//     Console.WriteLine();
//     
//     Console.WriteLine(new { h });
//     Console.WriteLine(h.Div(g2));
//     Console.WriteLine(h.Div(g2).rem.Div(g1));
//
//     // Console.WriteLine(g1.Div(g2));
//     // Console.WriteLine(g2.Div(g1));
// }
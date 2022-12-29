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


{
    Monom.Display = MonomDisplay.DotSuperscript;
    var (x1, x2, x3, x4) = Ring.Polynomial("X1", "X2", "X3", "X4", ZnInt.KZero());
    var lt = new[] { x1, x2 * x4, x3 };
    
    x1.Indeterminates.SetOrder(MonomOrder.Lex);
    Console.WriteLine(lt.Order().Glue(" < ", "{0,5}"));
    x1.Indeterminates.SetOrder(MonomOrder.RevLex);
    Console.WriteLine(lt.Order().Glue(" < ", "{0,5}"));
    Console.WriteLine();
}

{
    Monom.Display = MonomDisplay.DotSuperscript;
    var n = 3;
    var xs = Ring.Polynomial(ZnInt.KZero(), n.Range(1).Select(i => $"X{i}").ToArray());
    var lt = Partitions32[n].Select(l => l.Concat(Enumerable.Repeat(0, n - l.Count)).ToArray())
        .SelectMany(l => GetPermutations(n).Select(p => p.Select(i => l[i - 1]).ToArray()))
        .Distinct(new SequenceEquality<int>())
        .Select(l => l.Select((c, i) => xs[i].Pow(c)).Aggregate((a, b) => a * b)).ToArray();
    
    xs[0].Indeterminates.SetOrder(MonomOrder.GrLex);
    Console.WriteLine(lt.Order().Glue(" < ", "{0,8}"));
    xs[0].Indeterminates.SetOrder(MonomOrder.GrevLex);
    Console.WriteLine(lt.Order().Glue(" < ", "{0,8}"));
    Console.WriteLine();
}

{
    var (x, y) = Ring.Polynomial("X", "Y", Rational.KZero());
    var f = y * x.Pow(2) - x;
    var g1 = y * x - 1;
    var g2 = x.Pow(2);
    Console.WriteLine(new { f, g1, g2 });
    Console.WriteLine(x.Indeterminates);
    Console.WriteLine(f.Div(g1));
    Console.WriteLine(f.Div(g2));
    Console.WriteLine();
}

{
    var x = Ring.Polynomial(Rational.KZero());
    var ords = new[] { MonomOrder.Lex, MonomOrder.GrLex, MonomOrder.RevLex, MonomOrder.GrevLex };
    
    foreach (var o in ords)
    {
        x.Indeterminates.SetOrder(o);
        Console.WriteLine(x.Indeterminates);
        var X = x.Indeterminates[0];
        var f = 2 * x.Pow(3) + 3 * x + 1;
        var g = 7 * x.Pow(2) + x + 3;
        var h = f * g;
        Console.WriteLine(new { f, g, h });
        Console.WriteLine(new { ldf = f.LeadingDetails, ldg = g.LeadingDetails, ldh = h.LeadingDetails });
        Console.WriteLine(h.Div(g));
        Console.WriteLine(h.Div(f));
        Console.WriteLine(f.D(X));
        Console.WriteLine(g.D(X));
        Console.WriteLine(f.Div(x.Pow(2)));
        Console.WriteLine();
    }
}

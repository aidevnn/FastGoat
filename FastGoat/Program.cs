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

Polynomial<Rational, Xi> Primitive(Polynomial<Rational, Xi> f)
{
    if (f.IsZero())
        return f;
    
    var arrGcd = f.Coefs.Values.Where(e => !e.IsZero()).Select(e => e.Absolute.Num).Distinct().Order().ToArray();
    var arrLcm = f.Coefs.Values.Select(e => e.Absolute.Denom).Distinct().Order().ToArray();
    return f * new Rational(f.LeadingDetails.lc.Sign * IntExt.LcmBigInt(arrLcm), IntExt.GcdBigInt(arrGcd));
}

(Polynomial<Rational, Xi> F, Polynomial<Rational, Xi>[] facts)
    GenerateRandomPolynomialFxy(int nbFactors, int maxDegree, bool scalarLT = true, bool multiplicity = false)
{
    // TODO
    var (X2, X1) = Ring.Polynomial(Rational.KZero(), MonomOrder.Lex, "X2", "X1").Deconstruct();
    var x1s = (1 + maxDegree).Range().Select(X1.Pow).ToArray();
    var x2s = maxDegree.Range().Select(X2.Pow).ToArray();
    var x1x2s = x1s.Grid2D(x2s).Select(e => e.t1 * e.t2).Where(f => f.Degree <= maxDegree).Order().ToArray();
    var nbMonoms = x1x2s.Length;
    var nbNonScalarLT = scalarLT ? 0 : Rng.Next(1, nbFactors);
    var (x2, x1) = X1.Indeterminates.Deconstruct();

    Polynomial<Rational, Xi> Choose()
    {
        var f = (1 + IntExt.Rng.Next(nbMonoms)).Range()
            .Select(_ => x1x2s[IntExt.Rng.Next(nbMonoms)] * IntExt.Rng.Next(-5, 5))
            .Where(e => !e.IsZero())
            .Aggregate(X1.Zero, (acc, e) => acc + e);
        
        f = Primitive(f);
        var degLT = int.Min(maxDegree, f.Degree + 1);
        var i = IntExt.Rng.Next(degLT, maxDegree + 1);
        var j = !scalarLT && nbNonScalarLT > 0 ? IntExt.Rng.Next(degLT, maxDegree + 1) : 0;
        --nbNonScalarLT;
        var g = f + X2.Pow(int.Max(i, j)) * X1.Pow(int.Min(i, j));
        if (!g.ConstTerm.IsZero())
            return Primitive(g);
        
        var mn1 = g.Coefs.Keys.Aggregate(Monom<Xi>.Gcd);
        if (mn1.Degree == 0)
            return Primitive(g);

        return Primitive(g) + IntExt.RngSign * IntExt.Rng.Next(1, 5);
    }

    var c = multiplicity ? (maxDegree == 2 ? 3 : 2) : 1;
    var facts = 1000.Range().Select(_ => Choose()).Distinct().Take(nbFactors)
        .Select(f => f.Pow(IntExt.Rng.Next(1, 1 + c))).Order().ToArray();
    var F = facts.Aggregate((a0, a1) => a0 * a1);
    var degX1 = F.DegreeOf(x1);
    var degX2 = F.DegreeOf(x2);
    if (degX1 > degX2)
        throw new UnreachableException(); 

    return (F, facts);
}

bool BatchTest(int nbTest, int nbFactors, int maxDegree, bool scalarLT)
{
    for (int i = 0; i < nbTest; i++)
    {
        var (F, facts) = GenerateRandomPolynomialFxy(nbFactors, maxDegree, scalarLT, multiplicity: false);
        var (y, x) = F.Indeterminates.Deconstruct();
        var testLT = scalarLT ? F.CoefMax(y).Degree > 0 : F.CoefMax(y).Degree == 0;
        if (testLT || facts.Length != nbFactors || facts.Any(fi => fi.DegreeOf(y) > maxDegree))
        {
            Console.WriteLine($"At i = {i} on {nbTest}");
            Console.WriteLine($"F = {F}");
            facts.Println();
            Console.WriteLine();
            return false;
        }
    }

    return true;
}

{
    Ring.DisplayPolynomial = MonomDisplay.StarCaret;
    RngSeed(678124);
    var nbTest = 500;
    foreach (var (n, m) in 4.Range(1).SelectMany(m => (5 - m).Range(2).Select(n => (n, m))))
    {
        var testName = $"NbFactors:{n} MaxDegree:{m} ScalarLT";
        var result = BatchTest(nbTest, nbFactors: n, maxDegree: m, scalarLT: true);
        Console.WriteLine($"Test {testName,-40} :{(result ? "Success" : "Failure")}");
    }

    Console.WriteLine();
    foreach (var (n, m) in 4.Range(2).SelectMany(m => (6 - m).Range(2).Select(n => (n, m))))
    {
        var testName = $"NbFactors:{n} MaxDegree:{m} Non-ScalarLT";
        var result = BatchTest(nbTest, nbFactors: n, maxDegree: m, scalarLT: false);
        Console.WriteLine($"Test {testName,-40} :{(result ? "Success" : "Failure")}");
    }
}
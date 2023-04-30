using System.Numerics;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Polynoms;

namespace FastGoat.Examples;

public static class AlgebraicIntegerRelationLLL
{
    static AlgebraicIntegerRelationLLL()
    {
        Ring.DisplayPolynomial = MonomDisplay.StarCaret;
    }
    public static Rational[] AlphaBetaPolynomial(BigCplx alpha, BigCplx beta, int d, int O)
    {
        Console.WriteLine("Start LLL algorithm");
        var pi = BigReal.Pi(alpha.O);
        var N = BigReal.FromBigInteger(BigInteger.Pow(10, O), pi.O);
        var mat = new KMatrix<Rational>(Rational.KZero(), d, d).Zero;
        var tmp = alpha.One;
        for (int i = 0; i < mat.M - 1; i++)
        {
            mat.Coefs[i, i] = Rational.KOne();
            var ai = tmp;
            var aipi = ai.RealPart + pi * ai.ImaginaryPart; // Re(𝛼^i) + π * Im(𝛼^i)
            mat.Coefs[i, mat.N - 1] = (aipi * N).ToRational.RoundEven;
            tmp *= alpha;
        }

        var bpi = beta.RealPart + pi * beta.ImaginaryPart;
        mat.Coefs[mat.N - 1, mat.N - 1] = (bpi * N).ToRational.RoundEven;
        Console.WriteLine(mat);
        var lll = IntFactorisation.LLL(mat.T);
        Console.WriteLine();
        Console.WriteLine(lll);

        var col = lll.Cols.OrderBy(l => l.Aggregate(Rational.KZero(), (acc, v) => acc + v.Pow(2))).First();
        Console.WriteLine("End LLL algorithm");
        Console.WriteLine("Possible Solution");
        Console.WriteLine(col.T);
        Console.WriteLine();
        return col.SkipLast(1).Select(c => -c).ToArray();
    }

    public static Rational[] AlphaBetaPolynomial(BigReal alpha, BigReal beta, int d, int O)
    {
        return AlphaBetaPolynomial(BigCplx.FromBigReal(alpha), BigCplx.FromBigReal(beta), d, O);
    }

    public static void Example1()
    {
        var d = 8;
        var O = 40;
        var pi = BigReal.Pi(O);
        var beta = (5 * pi.Pow(2) / 24) * (3 * pi.Pow(4) - 28 * pi.Pow(2) - 24);
        Console.WriteLine(new { alpha = pi, beta });
        var coefs = AlphaBetaPolynomial(pi, beta, d, O);
        Console.WriteLine(coefs.Glue("; "));
        
        var poly = new KPoly<Rational>('x', Rational.KZero(), coefs.ToArray());
        Console.WriteLine(poly.Substitute(pi));
        Console.WriteLine(beta.ToBigReal(O));
        var fact = (poly.Substitute(pi) / beta);
        Console.WriteLine(fact);
        var P = poly / fact.ToRational.RoundEven;
        Console.WriteLine($"P = {P.SubstituteChar('π')}");
    }

    public static void Example2()
    {
        var d = 17;
        var O = 100;
        var alpha = BigReal.NthRoot(3, 4, 3 * O / 2) - BigReal.NthRoot(2, 4, 3 * O / 2);
        var beta = alpha.Pow(d - 1);
        Console.WriteLine(new { alpha, beta });
        var coefs = AlphaBetaPolynomial(alpha, beta, d, O);
        Console.WriteLine(coefs.Glue("; "));
        
        var poly = new KPoly<Rational>('x', Rational.KZero(), coefs.ToArray());
        Console.WriteLine(poly.Substitute(alpha));
        Console.WriteLine(beta.ToBigReal(O));
        var fact = (poly.Substitute(alpha) / beta);
        Console.WriteLine(fact);
        var P = poly.X.Pow(16) - poly / fact.ToRational.RoundEven;
        Console.WriteLine($"P = {P.SubstituteChar('X')}");

        var x = FG.QPoly();
        IntFactorisation.PrimitiveElt(x.Pow(4) - 2, x.Pow(4) - 3).Println(); // more faster
    }

    static EPoly<Rational>[] ConjugatesOfBeta(EPoly<Rational> scalar, BigCplx alpha, BigCplx beta, int d, int O)
    {
        GlobalStopWatch.AddLap();
        var coefs = AlphaBetaPolynomial(alpha, beta, d, O);
        var P = new KPoly<Rational>('x', Rational.KZero(), coefs.ToArray());
        var fact = (P.Substitute(alpha) / beta).ToBigCplx(O);
        P /= fact.RealPart.ToRational.RoundEven;
        Console.WriteLine(P);
        Console.WriteLine(P.Substitute(alpha));
        Console.WriteLine((alpha, "=> ", beta));
        Console.WriteLine();
        GlobalStopWatch.Show("Conjugates");
        
        if (!beta.ToBigCplx(O).Equals(P.Substitute(alpha).ToBigCplx(O)))
            throw new();
        
        var stack = new Stack<EPoly<Rational>>(new[] { P.Substitute(scalar.X) });
        for (int i = 0; i < d; i++)
        {
            var p0 = P.Substitute(stack.Peek());
            if (stack.Contains(p0))
                break;
        
            stack.Push(p0);
        }
        
        stack.Println("Conjugates");

        return stack.ToArray();
    }

    public static void Example3()
    {
        var x = FG.QPoly();
        var P = x.Pow(6) + 108; // S3
        
        var O = (P.Degree + 1).Pow(2);
        GlobalStopWatch.Restart();
        GlobalStopWatch.AddLap();
        var roots = FG.NRoots(P.ToBcPoly(3 * O / 2));
        GlobalStopWatch.Show("Roots");
        Console.WriteLine(P);
        roots.Println();

        var (X, y) = FG.EPolyXc(P, 'y');
        var alpha = roots[0];
        var beta = roots.First(r => !alpha.Equals(r) && !alpha.Equals(-r));
        var stack1 = ConjugatesOfBeta(y, alpha, beta, P.Degree + 1, O);

        var res = stack1.Select(fy => fy.Poly.Substitute(alpha)).ToArray();
        roots.Println();
        res.Println();
        roots.Where(r1 => res.All(r2 => (r1 - r2).Magnitude > 1e-14)).Println(); // roots.Except(res) wont work due to precision
        var beta2 = roots.Where(r1 => !alpha.Equals(r1) && res.All(r2 => (r1 - r2).Magnitude > 1e-14)).First();
        var stack2 = ConjugatesOfBeta(y, alpha, beta2, P.Degree + 1, O);

        var roots3 = stack1.Grid2D(stack2).Select(e => e.t1.Substitute(e.t2)).Order().ToHashSet();
        roots3.Println($"All Roots of P = {P}");
        Console.WriteLine("Prod[X - ri] = {0}", roots3.Aggregate(X.One, (acc, r) => acc * (X - r)));
        GlobalStopWatch.Show("END"); // Time:1478 ms
    }
    
    public static void Example4()
    {
        var x = FG.QPoly();
        var P = x.Pow(8) + 4324 * x.Pow(4) + 7496644; // D8

        var O = (P.Degree + 1).Pow(2);
        GlobalStopWatch.Restart();
        GlobalStopWatch.AddLap();
        var roots = FG.NRoots(P.ToBcPoly(3 * O / 2));
        GlobalStopWatch.Show("Roots");
        Console.WriteLine(P);
        roots.Println();

        var (X, y) = FG.EPolyXc(P, 'y');
        var alpha = roots[0];
        var beta1 = roots.First(r => !alpha.Equals(r));
        var stack1 = ConjugatesOfBeta(y, alpha, beta1, P.Degree + 1, O);

        var res1 = stack1.Select(fy => fy.Poly.Substitute(alpha)).ToArray();
        roots.Println();
        res1.Println();
        var remains1 = roots.Where(r1 => res1.All(r2 => (r1 - r2).Magnitude > 1e-14)).ToArray();
        remains1.Println(); // roots.Except(res) wont work due to precision
        var beta2 = remains1.First();
        var stack2 = ConjugatesOfBeta(y, alpha, beta2, P.Degree + 1, O);
        var roots3 = stack1.Grid2D(stack2).Select(e => e.t1.Substitute(e.t2)).Order().ToHashSet();
        roots3.Println("Roots");

        var res2 = roots3.Select(fy => fy.Poly.Substitute(alpha)).ToArray();
        roots.Println();
        res2.Println();
        var remains2 = remains1.Where(r1 => res2.All(r2 => (r1 - r2).Magnitude > 1e-14)).ToArray();
        remains2.Println(); // roots.Except(res) wont work due to precision
        var beta3 = remains2.First();
        var stack3 = ConjugatesOfBeta(y, alpha, beta3, P.Degree + 1, O);

        var roots4 = roots3.Grid2D(stack3).Select(e => e.t1.Substitute(e.t2)).Order().ToHashSet();
        var roots5 = roots4.Grid2D(roots4).Select(e => e.t1.Substitute(e.t2)).Order().ToHashSet();
        
        roots5.Println($"All Roots of P = {P}");
        Console.WriteLine("Prod[X - ri] = {0}", roots5.Aggregate(X.One, (acc, r) => acc * (X - r)));
        GlobalStopWatch.Show("END"); // Time:6028 ms
    }
    
    public static void Example5()
    {
        var x = FG.QPoly();
        var P = x.Pow(10) + -10 * x.Pow(8) + -75 * x.Pow(6) + 1500 * x.Pow(4) + -5500 * x.Pow(2) + 16000; // D10

        var O = (P.Degree + 1).Pow(2);
        GlobalStopWatch.Restart();
        GlobalStopWatch.AddLap();
        var roots = FG.NRoots(P.ToBcPoly(3 * O / 2));
        GlobalStopWatch.Show("Roots");
        Console.WriteLine(P);
        roots.Println();

        var (X, y) = FG.EPolyXc(P, 'y');
        var alpha = roots[0];
        var beta1 = roots.First(r => !alpha.Equals(r));
        var stack1 = ConjugatesOfBeta(y, alpha, beta1, P.Degree + 1, O);

        var res1 = stack1.Select(fy => fy.Poly.Substitute(alpha)).ToArray();
        roots.Println();
        res1.Println();
        var remains1 = roots.Where(r1 => res1.All(r2 => (r1 - r2).Magnitude > 1e-14)).ToArray();
        remains1.Println(); // roots.Except(res) wont work due to precision
        var beta2 = remains1.First();
        var stack2 = ConjugatesOfBeta(y, alpha, beta2, P.Degree + 1, O);
        var roots3 = stack1.Grid2D(stack2).Select(e => e.t1.Substitute(e.t2)).Order().ToHashSet();
        roots3.Println($"All Roots of P = {P}");
        Console.WriteLine("Prod[X - ri] = {0}", roots3.Aggregate(X.One, (acc, r) => acc * (X - r)));
        GlobalStopWatch.Show("END"); // Time:13124 ms
    }

    public static void Example6()
    {
        var x = FG.QPoly();
        var P = x.Pow(12) + 96 * x.Pow(8) + 1664 * x.Pow(6) + -16128 * x.Pow(4) + 165888 * x.Pow(2) + 331776; // A4

        var O = (P.Degree + 1).Pow(2);
        GlobalStopWatch.Restart();
        GlobalStopWatch.AddLap();
        var roots = FG.NRoots(P.ToBcPoly(3 * O / 2));
        GlobalStopWatch.Show("Roots");
        Console.WriteLine(P);
        roots.Println();

        var (X, y) = FG.EPolyXc(P, 'y');
        var alpha = roots[0];
        var beta1 = roots.First(r => !alpha.Equals(r));
        var stack1 = ConjugatesOfBeta(y, alpha, beta1, P.Degree + 1, O);

        var res1 = stack1.Select(fy => fy.Poly.Substitute(alpha)).ToArray();
        roots.Println();
        res1.Println();
        var remains1 = roots.Where(r1 => res1.All(r2 => (r1 - r2).Magnitude > 1e-14)).ToArray();
        remains1.Println(); // roots.Except(res) wont work due to precision
        var beta2 = remains1.First();
        var stack2 = ConjugatesOfBeta(y, alpha, beta2, P.Degree + 1, O);
        var roots3 = stack1.Grid2D(stack2).Select(e => e.t1.Substitute(e.t2)).Order().ToHashSet();
        roots3.Println("Roots");

        var res2 = roots3.Select(fy => fy.Poly.Substitute(alpha)).ToArray();
        roots.Println();
        res2.Println();
        var remains2 = remains1.Where(r1 => res2.All(r2 => (r1 - r2).Magnitude > 1e-14)).ToArray();
        remains2.Println(); // roots.Except(res) wont work due to precision
        var beta3 = remains2.First();
        var stack3 = ConjugatesOfBeta(y, alpha, beta3, P.Degree + 1, O);

        var roots4 = roots3.Grid2D(stack3).Select(e => e.t1.Substitute(e.t2)).Order().ToHashSet();
        var roots5 = roots4.Grid2D(roots4).Select(e => e.t1.Substitute(e.t2)).Order().ToHashSet();
        
        roots5.Println($"All Roots of P = {P}");
        Console.WriteLine("Prod[X - ri] = {0}", roots5.Aggregate(X.One, (acc, r) => acc * (X - r)));
        GlobalStopWatch.Show("END"); // Time:65353 ms
    }

    public static void Example7()
    {
        var x = FG.QPoly();
        var P = x.Pow(20) + 2500 * x.Pow(10) + 50000; // C5x:C4

        var O = 160; // enough
        GlobalStopWatch.Restart();
        GlobalStopWatch.AddLap();
        var roots = FG.NRoots(P.ToBcPoly(3 * O / 2));
        GlobalStopWatch.Show("Roots");
        Console.WriteLine(P);
        roots.Println();

        var (X, y) = FG.EPolyXc(P, 'y');
        var alpha = roots[0];
        var beta1 = roots.First(r => !alpha.Equals(r));
        var stack1 = ConjugatesOfBeta(y, alpha, beta1, P.Degree + 1, O);

        var res1 = stack1.Select(fy => fy.Poly.Substitute(alpha)).ToArray();
        roots.Println();
        res1.Println();
        var remains1 = roots.Where(r1 => res1.All(r2 => (r1 - r2).Magnitude > 1e-14)).ToArray();
        remains1.Println(); // roots.Except(res) wont work due to precision
        var beta2 = remains1.First();
        var stack2 = ConjugatesOfBeta(y, alpha, beta2, P.Degree + 1, O);
        var roots3 = stack1.Grid2D(stack2).Select(e => e.t1.Substitute(e.t2)).Order().ToHashSet();
        roots3.Println("Roots");
       
        var res2 = roots3.Select(fy => fy.Poly.Substitute(alpha)).ToArray();
        roots.Println();
        res2.Println();
        var remains2 = remains1.Where(r1 => res2.All(r2 => (r1 - r2).Magnitude > 1e-14)).ToArray();
        remains2.Println(); // roots.Except(res) wont work due to precision
        var beta3 = remains2.First();
        var stack3 = ConjugatesOfBeta(y, alpha, beta3, P.Degree + 1, O);

        var roots4 = roots3.Grid2D(stack3).Select(e => e.t1.Substitute(e.t2)).Order().ToHashSet();
        var roots5 = roots4.Grid2D(roots4).Select(e => e.t1.Substitute(e.t2)).Order().ToHashSet();
        
        roots5.Println($"All Roots of P = {P}");
        Console.WriteLine("Prod[X - ri] = {0}", roots5.Aggregate(X.One, (acc, r) => acc * (X - r)));
        GlobalStopWatch.Show("END"); // Time:160207 ms
    }

}
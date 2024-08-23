using System.Numerics;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Floats;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Polynoms;

namespace FastGoat.Examples;

public static class AlgebraicIntegerRelationPSLQ
{
    static AlgebraicIntegerRelationPSLQ()
    {
        Ring.DisplayPolynomial = MonomDisplay.StarCaret;
        GlobalStopWatch.Restart();
    }

    public static Rational[] AlphaBetaPolynomial(BigCplx alpha, BigCplx beta, int d, int O)
    {
        if (Logger.Level != LogLevel.Off)
            Console.WriteLine("Start PSLQM2 algorithm");

        var pi = BigReal.Pi(alpha.O);
        var mat = new KMatrix<BigReal>(pi.Zero.ToBigReal(O), 1, d).Zero;
        var ai = alpha.One;
        for (int i = 0; i < mat.N - 1; i++)
        {
            var aipi = ai.RealPart + pi * ai.ImaginaryPart; // Re(𝛼^i) + π * Im(𝛼^i)
            mat.Coefs[0, i] = aipi.ToBigReal(O);
            ai *= alpha;
        }

        mat.Coefs[0, mat.N - 1] = (beta.RealPart + pi * beta.ImaginaryPart).ToBigReal(O); // Re(β) + π * Im(β)

        if (Logger.Level != LogLevel.Off)
        {
            Console.WriteLine(mat);
            Console.WriteLine();
        }

        var gamma = 2 / BigReal.Sqrt(BigReal.FromBigInteger(3, O));
        var rel = PSLQM2<Dcml>.TwoLevelMultipair(mat, gamma);

        if (Logger.Level != LogLevel.Off)
        {
            Console.WriteLine("End PSLQM2 algorithm");
            Console.WriteLine("Possible Solution");
            Console.WriteLine(rel.ToKMatrix());
            Console.WriteLine();
        }

        return rel;
    }

    public static Rational[] AlphaBetaPolynomial(BigReal alpha, BigReal beta, int d, int O)
    {
        return AlphaBetaPolynomial(BigCplx.FromBigReal(alpha), BigCplx.FromBigReal(beta), d, O);
    }

    static KPoly<Rational> PslqMinPoly<T>(int r, int s, int O)
        where T : struct, IElt<T>, IRingElt<T>, IFieldElt<T>, IFloatElt<T>, IFixedPrecisionElt<T>
    {
        var n = r * s + 1; // Expected polynomial degree plus one
        var alpha = BigReal.NthRoot(3, r, O) - BigReal.NthRoot(2, s, O); // a = 3^(1/r) - 2^(1/s)
        var ai = n.Range().Select(k => alpha.Pow(k)).ToKMatrix();

        var gamma = 2 / BigReal.Sqrt(BigReal.FromBigInteger(3, O));

        GlobalStopWatch.AddLap();
        var coefs = PSLQM2<T>.TwoLevelMultipair(ai, gamma);
        Console.WriteLine(coefs.ToKMatrix());
        var P = FG.KPoly('X', coefs).Monic;
        GlobalStopWatch.Show($"Two level Multipair PSLQ<{typeof(T).Name}> min poly a = 3^(1/{r}) - 2^(1/{s})");
        Console.WriteLine($"P = {P} and P(a) = {P.Substitute(alpha).ToBigReal(3 * O / 4)}");
        Console.WriteLine();

        return P;
    }


    public static void Example1()
    {
        Logger.Level = LogLevel.Level2;
        var d = 8; // Expected polynomial degree plus one
        var O = 30; // maximum precision digits
        var pi = BigReal.Pi(O + 2 * d);
        var beta = BigReal.FromBigIntegerAndExponent(BigInteger.Parse("-1669947371922907049618724340073146784"), 1, O);
        var coefs = AlphaBetaPolynomial(pi, beta, d, O);

        var poly = new KPoly<Rational>('x', Rational.KZero(), coefs.SkipLast(1).ToArray());
        var fact = -coefs.Last();
        Console.WriteLine("factor : {0}", fact.ToDouble);
        var P = poly / fact;
        Console.WriteLine($"beta = {P.SubstituteChar('π')}");

        var betaExpected = (5 * pi.Pow(2) / 24) * (3 * pi.Pow(4) - 28 * pi.Pow(2) - 24);
        Console.WriteLine("Expected beta : {0:F24}", betaExpected.ToDcml.K);
        Console.WriteLine("Actual   beta : {0:F24}", P.Substitute(pi).ToDcml.K);
        Console.WriteLine("Are Equals {0}", betaExpected.ToBigReal(O).Equals(P.Substitute(pi).ToBigReal(O)));
    }

    public static void Example2()
    {
        Logger.Level = LogLevel.Level1;
        var d = 17; // Expected polynomial degree plus one
        var O1 = 80; // rounding digits
        var O2 = 120; // maximum precision digits
        var alpha = BigReal.NthRoot(3, 4, O2) - BigReal.NthRoot(2, 4, O2); // alpha = Qtrt(3) - Qtrt(2)
        var beta = alpha.Pow(d - 1);
        Console.WriteLine(new { alpha, beta });
        var coefs = AlphaBetaPolynomial(alpha, beta, d, O1);
        Console.WriteLine(coefs.Glue("; "));

        var poly = new KPoly<Rational>('x', Rational.KZero(), coefs.SkipLast(1).ToArray());
        var sum = poly.Substitute(alpha);
        var fact = -coefs.Last();
        Console.WriteLine("factor : {0}", fact.ToDouble);
        var P = poly.X.Pow(16) - poly / fact;
        Console.WriteLine($"P = {P.SubstituteChar('X')}");

        var x = FG.QPoly();
        IntFactorisation.PrimitiveElt(x.Pow(4) - 2, x.Pow(4) - 3).Println(); // more faster
    }

    static ConcreteGroup<KAut<Rational>> ConjugatesOfBeta(KAutGroup<Rational> bsKAutGroup, BigCplx alpha, BigCplx beta,
        int d, int O)
    {
        if (Logger.Level != LogLevel.Off)
            GlobalStopWatch.AddLap();

        var coefs = AlphaBetaPolynomial(alpha, beta, d, O);
        var P = new KPoly<Rational>('x', Rational.KZero(), coefs.SkipLast(1).ToArray());
        P /= -coefs.Last();
        if (Logger.Level != LogLevel.Off)
        {
            Console.WriteLine(P);
            Console.WriteLine(P.Substitute(alpha));
            Console.WriteLine((alpha, "=> ", beta));
            Console.WriteLine();
        }

        if (!beta.ToBigCplx(O).Equals(P.Substitute(alpha).ToBigCplx(O)))
            throw new();

        var fy = P.Substitute(bsKAutGroup.Neutral().E);
        var subGr = Group.Generate("Conjugates", bsKAutGroup, fy);
        if (Logger.Level != LogLevel.Off)
        {
            DisplayGroup.HeadElements(subGr);
            GlobalStopWatch.Show("Conjugates");
            Console.WriteLine();
        }

        return subGr;
    }

    public static ConcreteGroup<KAut<Rational>> GaloisGroupNumericRoots(BigCplx alpha, BigCplx[] cplxRoots,
        KPoly<Rational> P, int O)
    {
        P = P.SubstituteChar('y');
        var y = FG.EPoly(P, 'y');
        var kAut = new KAutGroup<Rational>(P);
        var subGrGal = Group.Generate(kAut, kAut.Neutral());

        if (P.Substitute(-y).IsZero())
            subGrGal = Group.Generate("Conjugates", kAut, -y);

        if (Logger.Level != LogLevel.Off)
        {
            Console.WriteLine(new { alpha, O });
            DisplayGroup.HeadElements(subGrGal);
            cplxRoots.Println("All Complex Roots");
        }

        while (subGrGal.Count() < cplxRoots.Length)
        {
            var remains = cplxRoots.ToHashSet();
            var sz = 0;
            while (sz != remains.Count)
            {
                sz = remains.Count;
                var tmp0 = subGrGal.Select(fy => fy.E.Poly.Substitute(alpha).ToBigCplx(O)).ToHashSet();
                var tmp1 = new HashSet<BigCplx>();
                foreach (var e in remains)
                {
                    var res0 = subGrGal.Select(fy => fy.E.Poly.Substitute(e).ToBigCplx(O)).ToHashSet();
                    if (res0.All(e1 => !tmp0.Contains(e1)))
                    {
                        tmp0.UnionWith(res0);
                        tmp1.Add(e);
                    }
                }

                remains = tmp1.ToHashSet();
            }

            if (Logger.Level != LogLevel.Off)
                remains.Println($"Remaining roots {remains.Count}");

            if (remains.Count == 0)
                break;

            var beta2 = remains.First(b => (alpha - b).Magnitude > 1e-12);
            var gr = ConjugatesOfBeta(kAut, alpha, beta2, P.Degree + 1, O);
            subGrGal = Group.DirectProduct("SubGr(Gal(P))", gr, subGrGal);
        }

        subGrGal.Name = "Gal( Q(y)/Q )";
        return subGrGal;
    }

    public static ConcreteGroup<KAut<Rational>> GaloisGroupPSLQ(KPoly<Rational> P, int O1, int O2)
    {
        if (Logger.Level != LogLevel.Off)
        {
            Console.WriteLine(P);
            GlobalStopWatch.AddLap();
        }

        var nRoots = FG.NRoots(P.ToBcPoly(O2));
        if (Logger.Level != LogLevel.Off)
            GlobalStopWatch.Show("Roots");

        var alpha = nRoots[0];
        return GaloisGroupNumericRoots(alpha, nRoots, P, O1);
    }

    public static void Example3()
    {
        var x = FG.QPoly();
        var P = x.Pow(3) + 2;
        var roots = IntFactorisation.SplittingField(P); // S3
        var minPoly = roots[0].F.SubstituteChar('X');

        var O1 = 20; // rounding digits
        var O2 = 30; // maximum precision digits

        GlobalStopWatch.AddLap();
        var galGr = GaloisGroupPSLQ(minPoly, O1, O2);
        DisplayGroup.HeadElements(galGr);
        var X = FG.KPoly('X', galGr.Neutral().E);
        Console.WriteLine("Prod[X - ri] = {0}", galGr.Aggregate(X.One, (acc, r) => acc * (X - r)));
        GlobalStopWatch.Show("END"); // Time:63ms
    }

    public static void Example4()
    {
        var x = FG.QPoly();
        var (minPoly, _, _, _) =
            IntFactorisation.PrimitiveElt(x.Pow(4) - 2, x.Pow(2) + 1).First(); // Gal(Q(i, √2)/Q) = D8

        var O1 = 30; // rounding digits
        var O2 = 40; // maximum precision digits

        GlobalStopWatch.AddLap();
        var galGr = GaloisGroupPSLQ(minPoly, O1, O2);
        DisplayGroup.HeadElements(galGr);
        var X = FG.KPoly('X', galGr.Neutral().E);
        Console.WriteLine("Prod[X - ri] = {0}", galGr.Aggregate(X.One, (acc, r) => acc * (X - r)));
        GlobalStopWatch.Show("END"); // Time:165ms
    }

    public static void Example5()
    {
        var x = FG.QPoly();
        var P = x.Pow(5) - 5 * x + 12;
        var roots = IntFactorisation.SplittingField(P); // D10
        var minPoly = roots[0].F.SubstituteChar('X');

        var O1 = 60; // rounding digits
        var O2 = 80; // maximum precision digits

        Logger.SetOff();
        GlobalStopWatch.AddLap();
        var galGr = GaloisGroupPSLQ(minPoly, O1, O2);
        DisplayGroup.HeadElements(galGr);
        var X = FG.KPoly('X', galGr.Neutral().E);
        Console.WriteLine("Prod[X - ri] = {0}", galGr.Aggregate(X.One, (acc, r) => acc * (X - r)));
        GlobalStopWatch.Show("END"); // Time:436ms
    }

    public static void Example6()
    {
        var x = FG.QPoly();
        var P = x.Pow(4) + 8 * x + 12;
        var roots = IntFactorisation.SplittingField(P); // A4
        var minPoly = roots[0].F.SubstituteChar('X');

        var O1 = 110; // rounding digits
        var O2 = 130; // maximum precision digits

        Logger.Level = LogLevel.Level1;
        GlobalStopWatch.AddLap();
        var galGr = GaloisGroupPSLQ(minPoly, O1, O2);
        DisplayGroup.HeadElements(galGr);
        var X = FG.KPoly('X', galGr.Neutral().E);
        Console.WriteLine("Prod[X - ri] = {0}", galGr.Aggregate(X.One, (acc, r) => acc * (X - r)));
        GlobalStopWatch.Show("END Roots"); // Time:1.417s
    }

    public static void Example7()
    {
        var x = FG.QPoly();
        var P = x.Pow(5) + 2;
        var roots = IntFactorisation.SplittingField(P); // C5x:C4
        var minPoly = roots[0].F.SubstituteChar('X');

        var O1 = 120; // rounding digits
        var O2 = 140; // maximum precision digits
        Logger.Level = LogLevel.Level1;
        GlobalStopWatch.AddLap();
        var galGr = GaloisGroupPSLQ(minPoly, O1, O2);
        DisplayGroup.HeadElements(galGr);
        var X = FG.KPoly('X', galGr.Neutral().E);
        Console.WriteLine("Prod[X - ri] = {0}", galGr.Aggregate(X.One, (acc, r) => acc * (X - r)));
        GlobalStopWatch.Show("END Roots"); // Time:8.123s
    }

    public static void Example8()
    {
        var x = FG.QPoly();
        // pari/gp
        // ? nfsplitting(x^6+3*x^3+3)
        // time = 11 ms. !!! pari is unbeatable
        // %1 = x^18 + 171*x^12 + 5130*x^6 + 27
        var P = x.Pow(18) + 171 * x.Pow(12) + 5130 * x.Pow(6) + 27; // C3x:C6

        var O1 = 140; // rounding digits
        var O2 = 160; // maximum precision digits

        Logger.Level = LogLevel.Level1;
        GlobalStopWatch.AddLap();
        var galGr = GaloisGroupPSLQ(P, O1, O2);
        DisplayGroup.HeadElements(galGr);
        var X = FG.KPoly('X', galGr.Neutral().E);
        Console.WriteLine("Prod[X - ri] = {0}", galGr.Aggregate(X.One, (acc, r) => acc * (X - r)));
        GlobalStopWatch.Show("END Roots"); // Time:10.873s
    }

    public static void Example9()
    {
        var x = FG.QPoly();
        /*
            pari/gp
            ? galoisgetname(21,1)
            %5 = "C7 : C3"
            ? galoisgetpol(21,1)
            %6 = [x^21 - 84*x^19 + 2436*x^17 - 31136*x^15 + 2312*x^14 + 203840*x^13 - 30688*x^12 - 733824*x^11 + 152992*x^10 + 1480192*x^9 - 359296*x^8 - 1628096*x^7 + 413952*x^6 + 892416*x^5 - 225792*x^4 - 189952*x^3 + 50176*x^2 + 3584*x - 512, 9219840]
         */
        var P = x.Pow(21) - 84 * x.Pow(19) + 2436 * x.Pow(17) - 31136 * x.Pow(15) + 2312 * x.Pow(14) +
            203840 * x.Pow(13) -
            30688 * x.Pow(12) - 733824 * x.Pow(11) + 152992 * x.Pow(10) + 1480192 * x.Pow(9) - 359296 * x.Pow(8) -
            1628096 * x.Pow(7) + 413952 * x.Pow(6) + 892416 * x.Pow(5) - 225792 * x.Pow(4) - 189952 * x.Pow(3) +
            50176 * x.Pow(2) +
            3584 * x - 512; // C7x:C3

        var O1 = 300; // rounding digits
        var O2 = 350; // maximum precision digits

        var galGr = GaloisGroupPSLQ(P, O1, O2);
        DisplayGroup.HeadElements(galGr);
        var X = FG.KPoly('X', galGr.Neutral().E);
        Console.WriteLine("Prod[X - ri] = {0}", galGr.Aggregate(X.One, (acc, r) => acc * (X - r)));
        GlobalStopWatch.Show("END Roots"); // Time:1m7s
    }

    public static void Example_PSLQM2_Dble_and_Dcml()
    {
        Logger.Level = LogLevel.Level1;
        Ring.DisplayPolynomial = MonomDisplay.StarCaret;

        PslqMinPoly<Dcml>(2, 2, 30);
        PslqMinPoly<Dble>(2, 2, 30);

        PslqMinPoly<Dcml>(3, 3, 50);
        PslqMinPoly<Dble>(3, 3, 50);

        PslqMinPoly<Dcml>(4, 4, 90);
        PslqMinPoly<Dble>(4, 4, 90);

        PslqMinPoly<Dcml>(5, 5, 180);
        PslqMinPoly<Dble>(5, 5, 180);

        PslqMinPoly<Dcml>(5, 6, 250);
        PslqMinPoly<Dble>(5, 6, 250);

        PslqMinPoly<Dcml>(6, 6, 320);
        PslqMinPoly<Dble>(6, 6, 320);
        Console.Beep();
    }
    // Possible Solution step:840
    // [697, -1440, -20520, -98280, -102060, -1458, 80, -43920, 538380, -336420, 1215, 0, -80, -56160, -135540, -540, 0, 0, 40, -7380, 135, 0, 0, 0, -10, -18, 0, 0, 0, 0, 1]
    // # Two level Multipair PSLQ<Dcml> min poly a = 3^(1/5) - 2^(1/6) Time:18.452s
    // P = X^30 - 18*X^25 - 10*X^24 + 135*X^20 - 7380*X^19 + 40*X^18 - 540*X^15 - 135540*X^14 - 56160*X^13 - 80*X^12 + 1215*X^10 - 336420*X^9 + 538380*X^8 - 43920*X^7 + 80*X^6 - 1458*X^5 - 102060*X^4 - 98280*X^3 - 20520*X^2 - 1440*X + 697 and P(a) = 0
    // 
    // Possible Solution step:840
    // [697, -1440, -20520, -98280, -102060, -1458, 80, -43920, 538380, -336420, 1215, 0, -80, -56160, -135540, -540, 0, 0, 40, -7380, 135, 0, 0, 0, -10, -18, 0, 0, 0, 0, 1]
    // # Two level Multipair PSLQ<Dble> min poly a = 3^(1/5) - 2^(1/6) Time:18.912s
    // P = X^30 - 18*X^25 - 10*X^24 + 135*X^20 - 7380*X^19 + 40*X^18 - 540*X^15 - 135540*X^14 - 56160*X^13 - 80*X^12 + 1215*X^10 - 336420*X^9 + 538380*X^8 - 43920*X^7 + 80*X^6 - 1458*X^5 - 102060*X^4 - 98280*X^3 - 20520*X^2 - 1440*X + 697 and P(a) = 0
    // 
    // Possible Solution step:1180
    // [1, 0, 0, 0, 0, 0, -4281690, 0, 0, 0, 0, 0, 3137919, 0, 0, 0, 0, 0, -618280, 0, 0, 0, 0, 0, -16221, 0, 0, 0, 0, 0, -30, 0, 0, 0, 0, 0, 1]
    // # Two level Multipair PSLQ<Dcml> min poly a = 3^(1/6) - 2^(1/6) Time:54.903s
    // P = X^36 - 30*X^30 - 16221*X^24 - 618280*X^18 + 3137919*X^12 - 4281690*X^6 + 1 and P(a) = 0
    // 
    // Possible Solution step:1180
    // [1, 0, 0, 0, 0, 0, -4281690, 0, 0, 0, 0, 0, 3137919, 0, 0, 0, 0, 0, -618280, 0, 0, 0, 0, 0, -16221, 0, 0, 0, 0, 0, -30, 0, 0, 0, 0, 0, 1]
    // # Two level Multipair PSLQ<Dble> min poly a = 3^(1/6) - 2^(1/6) Time:58.721s
    // P = X^36 - 30*X^30 - 16221*X^24 - 618280*X^18 + 3137919*X^12 - 4281690*X^6 + 1 and P(a) = 0
    // 
}
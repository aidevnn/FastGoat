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

public static class AlgebraicIntegerRelationLLL
{
    static AlgebraicIntegerRelationLLL()
    {
        Ring.DisplayPolynomial = MonomDisplay.StarCaret;
    }

    public static Rational[] AlphaBetaPolynomial(BigCplx alpha, BigCplx beta, int d, int O, bool details = true)
    {
        if (details)
            Console.WriteLine("Start LLL algorithm");

        var pi = BigReal.Pi(alpha.O);
        var N = BigReal.FromBigInteger(BigInteger.Pow(10, O), alpha.O);
        var mat = new KMatrix<BigReal>(pi.Zero.ToBigReal(O), d, d).Zero;
        var ai = alpha.One;
        var aipow = new List<BigCplx>();
        for (int i = 0; i < mat.M - 1; i++)
        {
            aipow.Add(ai);
            mat.Coefs[i, i] = pi.One;
            var aipi = ai.RealPart + pi * ai.ImaginaryPart; // Re(𝛼^i) + π * Im(𝛼^i)
            mat.Coefs[i, mat.N - 1] = (aipi * N).RoundEven;
            ai *= alpha;
        }

        var bpi = beta.RealPart + pi * beta.ImaginaryPart; // Re(β) + π * Im(β)
        mat.Coefs[mat.N - 1, mat.N - 1] = (bpi * N).RoundEven;
        if (details)
            Console.WriteLine(mat);
        var lll = IntFactorisation.LLL(mat.T);

        if (details)
        {
            Console.WriteLine();
            Console.WriteLine(lll);
        }

        var col = lll.Cols
            .OrderBy(l => l.SkipLast(1).Zip(aipow).Aggregate(-beta, (acc, v) => acc + BigCplx.FromBigReal(v.First) * v.Second)
                .Magnitude2).First();
        if (details)
        {
            Console.WriteLine("End LLL algorithm");
            Console.WriteLine("Possible Solution");
            Console.WriteLine(col.T);
            Console.WriteLine();
        }

        return col.SkipLast(1).Select(c => -c.RoundEven.ToRational).ToArray();
    }

    public static Rational[] AlphaBetaPolynomial(BigReal alpha, BigReal beta, int d, int O, bool details = true)
    {
        return AlphaBetaPolynomial(BigCplx.FromBigReal(alpha), BigCplx.FromBigReal(beta), d, O, details);
    }

    public static void Example1()
    {
        var d = 8; // Expected polynomial degree plus one
        var O1 = 20; // rounding digits
        var O2 = 30; // maximum precision digits
        var pi = BigReal.Pi(O2);
        var beta = BigReal.FromBigIntegerAndExponent(BigInteger.Parse("-1669947371922907049619"), 1, O2);
        var coefs = AlphaBetaPolynomial(pi, beta, d, O1);

        var poly = new KPoly<Rational>('x', Rational.KZero(), coefs.ToArray());
        var sum = poly.Substitute(pi);
        Console.WriteLine("Sum[ci*ai] : {0}", sum.ToDouble);
        Console.WriteLine("Actual beta : {0}", beta.ToDouble);
        var fact = (sum / beta);
        Console.WriteLine("factor : {0}", fact.ToDouble);
        var P = poly / fact.ToRational.RoundEven;
        Console.WriteLine($"beta = {P.SubstituteChar('π')}");

        var betaExpected = (5 * pi.Pow(2) / 24) * (3 * pi.Pow(4) - 28 * pi.Pow(2) - 24);
        Console.WriteLine("Expected beta : {0}", betaExpected.ToDouble);
        Console.WriteLine("Are Equals {0}", betaExpected.ToBigReal(O1).Equals(P.Substitute(pi).ToBigReal(O1)));
    }

    public static void Example2()
    {
        var d = 17; // Expected polynomial degree plus one
        var O1 = 80; // rounding digits
        var O2 = 120; // maximum precision digits
        var alpha = BigReal.NthRoot(3, 4, O2) - BigReal.NthRoot(2, 4, O2); // alpha = Qtrt(3) - Qtrt(2)
        var beta = alpha.Pow(d - 1);
        Console.WriteLine(new { alpha, beta });
        var coefs = AlphaBetaPolynomial(alpha, beta, d, O1);
        Console.WriteLine(coefs.Glue("; "));

        var poly = new KPoly<Rational>('x', Rational.KZero(), coefs.ToArray());
        var sum = poly.Substitute(alpha);
        Console.WriteLine("Sum[ci*ai] : {0}", sum.ToDouble);
        Console.WriteLine("Actual beta : {0}", beta.ToDouble);
        var fact = (sum / beta);
        Console.WriteLine("factor : {0}", fact.ToDouble);
        var P = poly.X.Pow(16) - poly / fact.ToRational.RoundEven;
        Console.WriteLine($"P = {P.SubstituteChar('X')}");

        var x = FG.QPoly();
        IntFactorisation.PrimitiveElt(x.Pow(4) - 2, x.Pow(4) - 3).Println(); // more faster
    }

    static ConcreteGroup<KAut<Rational>> ConjugatesOfBeta(KAutGroup<Rational> bsKAutGroup, BigCplx alpha, BigCplx beta, int d, int O,
        bool details = true)
    {
        if (details)
            GlobalStopWatch.AddLap();

        var coefs = AlphaBetaPolynomial(alpha, beta, d, O, details);
        var P = new KPoly<Rational>('x', Rational.KZero(), coefs.ToArray());
        var fact = (P.Substitute(alpha) / beta).ToBigCplx(O);
        P /= fact.RealPart.ToRational.RoundEven;
        if (details)
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
        if (details)
        {
            DisplayGroup.HeadElements(subGr);
            GlobalStopWatch.Show("Conjugates");
            Console.WriteLine();
        }

        return subGr;
    }

    public static ConcreteGroup<KAut<Rational>> GaloisGroupNumericRoots(BigCplx alpha, BigCplx[] cplxRoots, KPoly<Rational> P, int O,
        bool details = true)
    {
        P = P.SubstituteChar('y');
        var y = FG.EPoly(P, 'y');
        var kAut = new KAutGroup<Rational>(P);
        var subGrGal = Group.Generate(kAut, kAut.Neutral());

        if (P.Substitute(-y).IsZero())
            subGrGal = Group.Generate("Conjugates", kAut, -y);

        if (details)
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

            if (details)
                remains.Println($"Remaining roots {remains.Count}");

            if (remains.Count == 0)
                break;

            var beta2 = remains.First(b => (alpha - b).Magnitude > 1e-12);
            var gr = ConjugatesOfBeta(kAut, alpha, beta2, P.Degree + 1, O, details);
            subGrGal = Group.DirectProduct("SubGr(Gal(P))", gr, subGrGal);
        }

        subGrGal.SetName("Gal( Q(y)/Q )");
        return subGrGal;
    }

    public static ConcreteGroup<KAut<Rational>> GaloisGroupLLL(KPoly<Rational> P, int O1, int O2, bool details = true)
    {
        if (details)
        {
            Console.WriteLine(P);
            GlobalStopWatch.AddLap();
        }

        var nRoots = FG.NRoots(P.ToBcPoly(O2));
        if (details)
            GlobalStopWatch.Show("Roots");

        var alpha = nRoots[0];
        return GaloisGroupNumericRoots(alpha, nRoots, P, O1, details);
    }

    public static void Example3()
    {
        var x = FG.QPoly();
        var P = x.Pow(3) + 2;
        var roots = IntFactorisation.SplittingField(P, details: true); // S3
        var minPoly = roots[0].F.SubstituteChar('X');

        var O1 = 20; // rounding digits
        var O2 = 30; // maximum precision digits
        GlobalStopWatch.Restart();
        var galGr = GaloisGroupLLL(minPoly, O1, O2);
        DisplayGroup.HeadElements(galGr);
        var X = FG.KPoly('X', galGr.Neutral().E);
        Console.WriteLine("Prod[X - ri] = {0}", galGr.Aggregate(X.One, (acc, r) => acc * (X - r)));
        GlobalStopWatch.Show("END"); // Time:136 ms
        GlobalStopWatch.AddLap();
        GaloisApplications.GaloisCorrespondence(galGr.ToList());
        GlobalStopWatch.Show("END GaloisCorrespondence"); // Time:39 ms
        GlobalStopWatch.Show("END S3"); // Time:175 ms
    }

    public static void Example4()
    {
        var x = FG.QPoly();
        var (minPoly, _, _, _) = IntFactorisation.PrimitiveElt(x.Pow(4) - 2, x.Pow(2) + 1).First(); // Gal(Q(i, √2)/Q) = D8

        var O1 = 20; // rounding digits
        var O2 = 30; // maximum precision digits
        GlobalStopWatch.Restart();
        var galGr = GaloisGroupLLL(minPoly, O1, O2);
        DisplayGroup.HeadElements(galGr);
        var X = FG.KPoly('X', galGr.Neutral().E);
        Console.WriteLine("Prod[X - ri] = {0}", galGr.Aggregate(X.One, (acc, r) => acc * (X - r)));
        GlobalStopWatch.Show("END"); // Time:266 ms
        GlobalStopWatch.AddLap();
        GaloisApplications.GaloisCorrespondence(galGr.ToList());
        GlobalStopWatch.Show("END GaloisCorrespondence"); // Time:166 ms
        GlobalStopWatch.Show("END D8"); // Time:432 ms
    }

    public static void Example5()
    {
        var x = FG.QPoly();
        var P = x.Pow(5) - 5 * x + 12;
        var roots = IntFactorisation.SplittingField(P, details: true); // D10
        var minPoly = roots[0].F.SubstituteChar('X');

        var O1 = 40; // rounding digits
        var O2 = 50; // maximum precision digits
        GlobalStopWatch.Restart();
        var galGr = GaloisGroupLLL(minPoly, O1, O2);
        DisplayGroup.HeadElements(galGr);
        var X = FG.KPoly('X', galGr.Neutral().E);
        Console.WriteLine("Prod[X - ri] = {0}", galGr.Aggregate(X.One, (acc, r) => acc * (X - r)));
        GlobalStopWatch.Show("END"); // Time:814 ms
        GlobalStopWatch.AddLap();
        GaloisApplications.GaloisCorrespondence(galGr.ToList());
        GlobalStopWatch.Show("END GaloisCorrespondence"); // Time:449 ms
        GlobalStopWatch.Show("END D10"); // Time:1263 ms
    }

    public static void Example6()
    {
        var x = FG.QPoly();
        var P = x.Pow(4) + 8 * x + 12;
        var roots = IntFactorisation.SplittingField(P, details: true); // A4
        var minPoly = roots[0].F.SubstituteChar('X');

        var O1 = 110; // rounding digits
        var O2 = 130; // maximum precision digits
        GlobalStopWatch.Restart();
        var galGr = GaloisGroupLLL(minPoly, O1, O2);
        DisplayGroup.HeadElements(galGr);
        var X = FG.KPoly('X', galGr.Neutral().E);
        Console.WriteLine("Prod[X - ri] = {0}", galGr.Aggregate(X.One, (acc, r) => acc * (X - r)));
        GlobalStopWatch.Show("END Roots"); // Time:4956 ms
        GlobalStopWatch.AddLap();
        GaloisApplications.GaloisCorrespondence(galGr.ToList());
        GlobalStopWatch.Show("END GaloisCorrespondence"); // Time:2302 ms
        GlobalStopWatch.Show("END A4"); // Time:7258 ms
    }

    public static void Example7()
    {
        var x = FG.QPoly();
        var P = x.Pow(5) + 2;
        var roots = IntFactorisation.SplittingField(P, details: true); // C5x:C4
        var minPoly = roots[0].F.SubstituteChar('X');

        var O1 = 80; // rounding digits
        var O2 = 100; // maximum precision digits
        GlobalStopWatch.Restart();
        var galGr = GaloisGroupLLL(minPoly, O1, O2);
        DisplayGroup.HeadElements(galGr);
        var X = FG.KPoly('X', galGr.Neutral().E);
        Console.WriteLine("Prod[X - ri] = {0}", galGr.Aggregate(X.One, (acc, r) => acc * (X - r)));
        GlobalStopWatch.Show("END Roots"); // Time:18342 ms
        GlobalStopWatch.AddLap();
        GaloisApplications.GaloisCorrespondence(galGr.ToList());
        GlobalStopWatch.Show("END GaloisCorrespondence"); // Time:5462 ms
        GlobalStopWatch.Show("END C5x:C4"); // Time:23804 ms
    }

    public static void Example8()
    {
        var x = FG.QPoly();
        // pari/gp
        // ? nfsplitting(x^6+3*x^3+3)
        // time = 11 ms. !!! pari is unbeatable
        // %1 = x^18 + 171*x^12 + 5130*x^6 + 27
        var P = x.Pow(18) + 171 * x.Pow(12) + 5130 * x.Pow(6) + 27; // C3x:C6

        var O1 = 120; // rounding digits
        var O2 = 140; // maximum precision digits
        GlobalStopWatch.Restart();
        var galGr = GaloisGroupLLL(P, O1, O2);
        DisplayGroup.HeadElements(galGr);
        var X = FG.KPoly('X', galGr.Neutral().E);
        Console.WriteLine("Prod[X - ri] = {0}", galGr.Aggregate(X.One, (acc, r) => acc * (X - r)));
        GlobalStopWatch.Show("END Roots"); // Time:29501 ms
        GlobalStopWatch.AddLap();
        GaloisApplications.GaloisCorrespondence(galGr.ToList());
        GlobalStopWatch.Show("END GaloisCorrespondence"); // Time:5758 ms
        GlobalStopWatch.Show("END C3x:C6"); // Time:35259 ms
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
        var P = x.Pow(21) - 84 * x.Pow(19) + 2436 * x.Pow(17) - 31136 * x.Pow(15) + 2312 * x.Pow(14) + 203840 * x.Pow(13) -
            30688 * x.Pow(12) - 733824 * x.Pow(11) + 152992 * x.Pow(10) + 1480192 * x.Pow(9) - 359296 * x.Pow(8) -
            1628096 * x.Pow(7) + 413952 * x.Pow(6) + 892416 * x.Pow(5) - 225792 * x.Pow(4) - 189952 * x.Pow(3) + 50176 * x.Pow(2) +
            3584 * x - 512; // C7x:C3

        var O1 = 250; // rounding digits
        var O2 = 300; // maximum precision digits
        GlobalStopWatch.Restart();
        var galGr = GaloisGroupLLL(P, O1, O2);
        DisplayGroup.HeadElements(galGr);
        var X = FG.KPoly('X', galGr.Neutral().E);
        Console.WriteLine("Prod[X - ri] = {0}", galGr.Aggregate(X.One, (acc, r) => acc * (X - r)));
        GlobalStopWatch.Show("END Roots"); // Time:245963 ms ~ 4 min
        GlobalStopWatch.AddLap();
        GaloisApplications.GaloisCorrespondence(galGr.ToList());
        GlobalStopWatch.Show("END GaloisCorrespondence"); // Time:119329 ms ~ 2 min
        GlobalStopWatch.Show("END C7x:C3"); // Time:365292 ms ~ 6 min
    }
/*
With P = X^21 + -84*X^19 + 2436*X^17 + -31136*X^15 + 2312*X^14 + 203840*X^13 + -30688*X^12 + -733824*X^11 + 152992*X^10 + 1480192*X^9 + -359296*X^8 + -1628096*X^7 + 413952*X^6 + 892416*X^5 + -225792*X^4 + -189952*X^3 + 50176*X^2 + 3584*X + -512
Factorization in Q(y)/Q
    X + -y
    X + 1361/1843968*y^20 + 3533/9219840*y^19 + -94909/1536640*y^18 + -3513/109760*y^17 + 2049643/1152480*y^16 + 2116697/2304960*y^15 + -25847291/1152480*y^14 + -94067/9604*y^13 + 16659499/115248*y^12 + 14613649/288120*y^11 + -2005019/3920*y^10 + -8219497/57624*y^9 + 19433377/19208*y^8 + 32806061/144060*y^7 + -31322509/28812*y^6 + -29227601/144060*y^5 + 20857639/36015*y^4 + 299659/3430*y^3 + -4321439/36015*y^2 + -372361/36015*y + 5540/2401
    X + -1257/3073280*y^20 + -181/921984*y^19 + 6623/192080*y^18 + 10819/658560*y^17 + -2320181/2304960*y^16 + -5678/12005*y^15 + 312741/24010*y^14 + 962943/192080*y^13 + -50350081/576240*y^12 + -7259411/288120*y^11 + 26970427/82320*y^10 + 2340124/36015*y^9 + -33746919/48020*y^8 + -11894461/144060*y^7 + 60642787/72030*y^6 + 1316096/36015*y^5 + -18449351/36015*y^4 + 77582/5145*y^3 + 857587/7203*y^2 + -164509/12005*y + -33622/36015
    X + -883/2304960*y^20 + -341/2304960*y^19 + 12379/384160*y^18 + 137/10976*y^17 + -2159081/2304960*y^16 + -420403/1152480*y^15 + 432853/36015*y^14 + 1474283/384160*y^13 + -45718607/576240*y^12 + -5634619/288120*y^11 + 7938703/27440*y^10 + 1984102/36015*y^9 + -7109447/12005*y^8 + -2643583/28812*y^7 + 9452249/14406*y^6 + 6495973/72030*y^5 + -25148009/72030*y^4 + -30455/686*y^3 + 2242091/36015*y^2 + 210376/36015*y + 8438/12005
    X + 2123/4609920*y^20 + 643/9219840*y^19 + -178529/4609920*y^18 + -3937/658560*y^17 + 2593627/2304960*y^16 + 414373/2304960*y^15 + -16626439/1152480*y^14 + -797029/576240*y^13 + 9114871/96040*y^12 + 94387/28812*y^11 + -14144113/41160*y^10 + 292613/288120*y^9 + 25049828/36015*y^8 + -84219/12005*y^7 + -109302059/144060*y^6 + -94943/9604*y^5 + 2852819/7203*y^4 + 224773/10290*y^3 + -2480644/36015*y^2 + -54828/12005*y + -8606/12005
    X + 991/3073280*y^20 + 361/9219840*y^19 + -125513/4609920*y^18 + -727/219520*y^17 + 612487/768320*y^16 + 44851/460992*y^15 + -1494211/144060*y^14 + -609467/1152480*y^13 + 13428531/192080*y^12 + -111101/72030*y^11 + -21652387/82320*y^10 + 1947119/96040*y^9 + 20362934/36015*y^8 + -297529/4802*y^7 + -97744189/144060*y^6 + 11987431/144060*y^5 + 29797223/72030*y^4 + -274373/5145*y^3 + -1150434/12005*y^2 + 526567/36015*y + 22688/36015
    X + -6721/9219840*y^20 + -1363/9219840*y^19 + 281269/4609920*y^18 + 2719/219520*y^17 + -506389/288120*y^16 + -274781/768320*y^15 + 4260759/192080*y^14 + 1645529/576240*y^13 + -41101813/288120*y^12 + -147939/19208*y^11 + 10314869/20580*y^10 + 369707/288120*y^9 + -281694509/288120*y^8 + 2242813/144060*y^7 + 147850729/144060*y^6 + 27199/9604*y^5 + -1266466/2401*y^4 + -39049/1470*y^3 + 3723359/36015*y^2 + 329444/36015*y + -71662/36015
    X + 5771/9219840*y^20 + 3697/9219840*y^19 + -16141/307328*y^18 + -7307/219520*y^17 + 175211/115248*y^16 + 2175949/2304960*y^15 + -5583703/288120*y^14 + -1943043/192080*y^13 + 9190667/72030*y^12 + 14850167/288120*y^11 + -6437313/13720*y^10 + -38396917/288120*y^9 + 95286763/96040*y^8 + 4713473/28812*y^7 + -171817759/144060*y^6 + -8088463/144060*y^5 + 27367502/36015*y^4 + -162177/3430*y^3 + -7337941/36015*y^2 + 1022593/36015*y + 63514/12005
    X + -3/31360*y^20 + -517/1317120*y^19 + 251/31360*y^18 + 21559/658560*y^17 + -37921/164640*y^16 + -20555/21952*y^15 + 159077/54880*y^14 + 89559/7840*y^13 + -1594051/82320*y^12 + -410171/5880*y^11 + 6163963/82320*y^10 + 9327161/41160*y^9 + -587907/3430*y^8 + -1596263/4116*y^7 + 659249/2940*y^6 + 6581633/20580*y^5 + -220873/1470*y^4 + -509608/5145*y^3 + 204959/5145*y^2 + 1902/1715*y + -3196/5145
    X + -5687/9219840*y^20 + -509/3073280*y^19 + 238657/4609920*y^18 + 75/5488*y^17 + -1728023/1152480*y^16 + -176719/460992*y^15 + 4410415/230496*y^14 + 915319/288120*y^13 + -72431467/576240*y^12 + -2455603/288120*y^11 + 7555591/16464*y^10 + -2330861/288120*y^9 + -276054143/288120*y^8 + 11481509/144060*y^7 + 54351981/48020*y^6 + -20731153/144060*y^5 + -25506853/36015*y^4 + 222433/2058*y^3 + 6739861/36015*y^2 + -71034/2401*y + -179104/36015
    X + 323/9219840*y^20 + 109/384160*y^19 + -6931/2304960*y^18 + -97/4116*y^17 + 34879/384160*y^16 + 128477/192080*y^15 + -481003/384160*y^14 + -9351791/1152480*y^13 + 5725343/576240*y^12 + 3504839/72030*y^11 + -571463/11760*y^10 + -7262659/48020*y^9 + 10415233/72030*y^8 + 1686271/7203*y^7 + -1776245/7203*y^6 + -10546027/72030*y^5 + 15471991/72030*y^4 + 2729/2058*y^3 + -881308/12005*y^2 + 697001/36015*y + 23458/12005
    X + 479/2304960*y^20 + -283/9219840*y^19 + -15941/921984*y^18 + 1819/658560*y^17 + 227173/460992*y^16 + -206581/2304960*y^15 + -7040507/1152480*y^14 + 1058191/576240*y^13 + 3680231/96040*y^12 + -600476/36015*y^11 + -1350941/10290*y^10 + 21051853/288120*y^9 + 73572629/288120*y^8 + -779199/4802*y^7 + -41770369/144060*y^6 + 8665179/48020*y^5 + 7057502/36015*y^4 + -957001/10290*y^3 + -2451926/36015*y^2 + 178656/12005*y + 48494/12005
    X + -1399/4609920*y^20 + -139/384160*y^19 + 58301/2304960*y^18 + 47/1568*y^17 + -1666759/2304960*y^16 + -975847/1152480*y^15 + 2595281/288120*y^14 + 11048791/1152480*y^13 + -33239743/576240*y^12 + -15294571/288120*y^11 + 16971491/82320*y^10 + 22119137/144060*y^9 + -122956817/288120*y^8 + -6498647/28812*y^7 + 2469119/4802*y^6 + 10037947/72030*y^5 + -24937781/72030*y^4 + -12443/2058*y^3 + 3886639/36015*y^2 + -177347/12005*y + -167854/36015
    X + 1357/9219840*y^20 + 613/2304960*y^19 + -28237/2304960*y^18 + -14677/658560*y^17 + 40217/115248*y^16 + 185309/288120*y^15 + -154859/36015*y^14 + -8981573/1152480*y^13 + 2582917/96040*y^12 + 6891419/144060*y^11 + -3740881/41160*y^10 + -23138261/144060*y^9 + 23650649/144060*y^8 + 3580343/12005*y^7 + -10159843/72030*y^6 + -3519266/12005*y^5 + 490453/14406*y^4 + 199879/1470*y^3 + 372578/36015*y^2 + -232651/12005*y + -12356/12005
    X + 69/1536640*y^20 + -843/3073280*y^19 + -3863/921984*y^18 + 187/8232*y^17 + 165901/1152480*y^16 + -1478069/2304960*y^15 + -683209/288120*y^14 + 1816903/230496*y^13 + 2919331/144060*y^12 + -14114183/288120*y^11 + -1305869/13720*y^10 + 9446291/57624*y^9 + 71801441/288120*y^8 + -42742559/144060*y^7 + -50433703/144060*y^6 + 13263629/48020*y^5 + 3322027/14406*y^4 + -40765/343*y^3 + -1610507/36015*y^2 + 734444/36015*y + -17246/7203
    X + -1063/3073280*y^20 + 1/9219840*y^19 + 132721/4609920*y^18 + -3/31360*y^17 + -630677/768320*y^16 + 16463/2304960*y^15 + 5856629/576240*y^14 + -1131791/1152480*y^13 + -12161673/192080*y^12 + 1707767/144060*y^11 + 17443367/82320*y^10 + -5167703/96040*y^9 + -111067781/288120*y^8 + 1288394/12005*y^7 + 53328239/144060*y^6 + -12583961/144060*y^5 + -12845971/72030*y^4 + 91741/5145*y^3 + 493756/12005*y^2 + 83053/36015*y + -96976/36015
    X + 61/9219840*y^20 + 227/1536640*y^19 + -131/2304960*y^18 + -577/47040*y^17 + -4863/192080*y^16 + 133711/384160*y^15 + 74531/76832*y^14 + -980183/230496*y^13 + -7338179/576240*y^12 + 3899321/144060*y^11 + 1303483/16464*y^10 + -940595/9604*y^9 + -73884331/288120*y^8 + 2983747/14406*y^7 + 31208609/72030*y^6 + -8987662/36015*y^5 + -24980647/72030*y^4 + 1602887/10290*y^3 + 1139266/12005*y^2 + -278491/7203*y + 1786/2401
    X + 17/219520*y^20 + -103/1317120*y^19 + -4283/658560*y^18 + 4283/658560*y^17 + 31039/164640*y^16 + -60919/329280*y^15 + -396449/164640*y^14 + 404113/164640*y^13 + 1281517/82320*y^12 + -223369/13720*y^11 + -4472117/82320*y^10 + 2309347/41160*y^9 + 531641/5145*y^8 + -2032649/20580*y^7 + -703789/6860*y^6 + 1652543/20580*y^5 + 160961/3430*y^4 + -116026/5145*y^3 + -34079/5145*y^2 + 11701/5145*y + -24/1715
    X + 23/9219840*y^20 + 101/1152480*y^19 + -799/1152480*y^18 + -4843/658560*y^17 + 3551/76832*y^16 + 61261/288120*y^15 + -700351/576240*y^14 + -3096067/1152480*y^13 + 527003/36015*y^12 + 765987/48020*y^11 + -1227223/13720*y^10 + -521247/12005*y^9 + 83690477/288120*y^8 + 1302851/36015*y^7 + -35589517/72030*y^6 + 1782758/36015*y^5 + 1913463/4802*y^4 + -325751/3430*y^3 + -1338626/12005*y^2 + 449846/12005*y + -15352/36015
    X + 853/3073280*y^20 + 289/921984*y^19 + -53099/2304960*y^18 + -17141/658560*y^17 + 1505659/2304960*y^16 + 425581/576240*y^15 + -2305213/288120*y^14 + -4846991/576240*y^13 + 28608269/576240*y^12 + 4556593/96040*y^11 + -13817173/82320*y^10 + -20695049/144060*y^9 + 91102031/288120*y^8 + 33816689/144060*y^7 + -7885081/24010*y^6 + -6950864/36015*y^5 + 2197848/12005*y^4 + 337102/5145*y^3 + -368159/7203*y^2 + -243892/36015*y + 36306/12005
    X + -97/1536640*y^20 + -1811/9219840*y^19 + 26231/4609920*y^18 + 5441/329280*y^17 + -42815/230496*y^16 + -1106639/2304960*y^15 + 329831/115248*y^14 + 6909449/1152480*y^13 + -6932531/288120*y^12 + -711663/19208*y^11 + 476353/4116*y^10 + 34224101/288120*y^9 + -91413733/288120*y^8 + -27355189/144060*y^7 + 4530489/9604*y^6 + 3569669/28812*y^5 + -8017577/24010*y^4 + -14159/5145*y^3 + 2806667/36015*y^2 + -122519/7203*y + 21118/12005

Galois Group
|G1| = 21
Type        NonAbelianGroup
BaseGroup   S21
SuperGroup  |Gal( Q(y)/Q )| = 21

Elements
( 1)[1] = []
( 2)[3] = [(1 3 10)(2 20 13)(4 15 12)(5 11 14)(6 18 8)(7 17 19)(9 21 16)]
( 3)[3] = [(1 6 7)(2 12 8)(3 5 13)(4 10 16)(9 11 18)(14 15 19)(17 20 21)]
( 4)[3] = [(1 7 6)(2 8 12)(3 13 5)(4 16 10)(9 18 11)(14 19 15)(17 21 20)]
( 5)[3] = [(1 8 13)(2 15 17)(3 14 4)(5 20 9)(6 19 11)(7 10 21)(12 18 16)]
( 6)[3] = [(1 9 15)(2 10 11)(3 18 17)(4 6 20)(5 12 7)(8 21 14)(13 16 19)]
( 7)[3] = [(1 10 3)(2 13 20)(4 12 15)(5 14 11)(6 8 18)(7 19 17)(9 16 21)]
( 8)[3] = [(1 11 12)(2 16 14)(3 9 19)(4 7 13)(5 8 17)(6 21 15)(10 18 20)]
( 9)[3] = [(1 12 11)(2 14 16)(3 19 9)(4 13 7)(5 17 8)(6 15 21)(10 20 18)]
(10)[3] = [(1 13 8)(2 17 15)(3 4 14)(5 9 20)(6 11 19)(7 21 10)(12 16 18)]
(11)[3] = [(1 14 20)(2 9 7)(3 21 12)(4 17 11)(5 6 16)(8 19 10)(13 15 18)]
(12)[3] = [(1 15 9)(2 11 10)(3 17 18)(4 20 6)(5 7 12)(8 14 21)(13 19 16)]
(13)[3] = [(1 16 17)(2 3 6)(4 8 9)(5 15 10)(7 14 18)(11 13 21)(12 19 20)]
(14)[3] = [(1 17 16)(2 6 3)(4 9 8)(5 10 15)(7 18 14)(11 21 13)(12 20 19)]
(15)[3] = [(1 20 14)(2 7 9)(3 12 21)(4 11 17)(5 16 6)(8 10 19)(13 18 15)]
(16)[7] = [(1 2 19 18 4 5 21)(3 15 20 11 7 16 8)(6 14 10 17 12 9 13)]
(17)[7] = [(1 4 2 5 19 21 18)(3 7 15 16 20 8 11)(6 12 14 9 10 13 17)]
(18)[7] = [(1 5 18 2 21 4 19)(3 16 11 15 8 7 20)(6 9 17 14 13 12 10)]
(19)[7] = [(1 18 21 19 5 2 4)(3 11 8 20 16 15 7)(6 17 13 10 9 14 12)]
(20)[7] = [(1 19 4 21 2 18 5)(3 20 7 8 15 11 16)(6 10 12 13 14 17 9)]
(21)[7] = [(1 21 5 4 18 19 2)(3 8 16 7 11 20 15)(6 13 9 12 17 10 14)]

Tower 1
  |G10| = 1  => [Q(y):Q] = 21
  |G3| = 3   => [Q(b):Q] = 7 with b=y^15 + -385919586248/272963387143*y^14 + -22423754495016/272963387143*y^13 + 30436879698512/272963387143*y^12 + 621842228060172/272963387143*y^11 + -25286444780464/8805270553*y^10 + -7224879886459704/272963387143*y^9 + 278006675821712/8805270553*y^8 + 38334910679506448/272963387143*y^7 + -350747216514240/2149318009*y^6 + -92019605722181760/272963387143*y^5 + 109501814104171328/272963387143*y^4 + 79462071607347520/272963387143*y^3 + -114624157244216704/272963387143*y^2 + 22778674912316608/272963387143*y
  |G1| = 21  => [Q:Q]    = 1
Tower 2
  |G10| = 1  => [Q(y):Q] = 21
  |G4| = 3   => [Q(c):Q] = 7 with c=y^15 + 2178925/39979708*y^14 + -794116848/9994927*y^13 + -605460268/129934051*y^12 + 20665642956/9994927*y^11 + 1350830236/9994927*y^10 + -211638240324/9994927*y^9 + 7474347344/9994927*y^8 + 931498385636/9994927*y^7 + -140208336912/9994927*y^6 + -1646403425296/9994927*y^5 + 444967592832/9994927*y^4 + 816724473600/9994927*y^3 + -231278938912/9994927*y^2 + -16607139840/9994927*y
  |G1| = 21  => [Q:Q]    = 1
Tower 3
  |G10| = 1  => [Q(y):Q] = 21
  |G5| = 3   => [Q(d):Q] = 7 with d=y^15 + -14027780265/2034773137*y^14 + -209195730541/2034773137*y^13 + 1073631375012/2034773137*y^12 + 7937243721288/2034773137*y^11 + -25773444164230/2034773137*y^10 + -137766297961824/2034773137*y^9 + 218813116574544/2034773137*y^8 + 1073678700615264/2034773137*y^7 + -637848173900504/2034773137*y^6 + -3724615546219456/2034773137*y^5 + -44206436042384/2034773137*y^4 + 47237654832816/19755079*y^3 + 106637631252992/107093323*y^2 + -461208255634368/2034773137*y
  |G1| = 21  => [Q:Q]    = 1
Tower 4
  |G10| = 1  => [Q(y):Q] = 21
  |G6| = 3   => [Q(e):Q] = 7 with e=y^15 + 270755177/662021746*y^14 + -25888776599/331010873*y^13 + -9833761441/331010873*y^12 + 652879562034/331010873*y^11 + 208167820434/331010873*y^10 + -6240734903060/331010873*y^9 + -310605071410/331010873*y^8 + 24183941335064/331010873*y^7 + -9909306335152/331010873*y^6 + -27654204706968/331010873*y^5 + 50279379125024/331010873*y^4 + -36104543030416/331010873*y^3 + -66910847658144/331010873*y^2 + 71748660854928/331010873*y
  |G1| = 21  => [Q:Q]    = 1
Tower 5
  |G10| = 1  => [Q(y):Q] = 21
  |G7| = 3   => [Q(f):Q] = 7 with f=y^15 + -18613119777/536743812785*y^14 + -43213394536676/536743812785*y^13 + 1791623213478/536743812785*y^12 + 1156023637377644/536743812785*y^11 + -12479589658788/107348762557*y^10 + -12622016704072216/536743812785*y^9 + 2212137578021912/536743812785*y^8 + 63778266805028496/536743812785*y^7 + -3782704466332736/107348762557*y^6 + -152282863133656976/536743812785*y^5 + 58331789241131168/536743812785*y^4 + 158482304637622784/536743812785*y^3 + -62523285303412864/536743812785*y^2 + -58211461304048896/536743812785*y
  |G1| = 21  => [Q:Q]    = 1
Tower 6
  |G10| = 1  => [Q(y):Q] = 21
  |G8| = 3   => [Q(g):Q] = 7 with g=y^15 + 93228046639/7778131624*y^14 + -433585876969/5833598718*y^13 + -22326743691343/23334394872*y^12 + 6412272774103/3889065812*y^11 + 146724427249813/5833598718*y^10 + -56502013826951/5833598718*y^9 + -254324375545048/972266453*y^8 + -18761338036175/1944532906*y^7 + 3552424377316738/2916799359*y^6 + 252895110426000/972266453*y^5 + -2454417480518682/972266453*y^4 + -2281166794612426/2916799359*y^3 + 5753351366491568/2916799359*y^2 + 1700355995323100/2916799359*y
  |G1| = 21  => [Q:Q]    = 1
Tower 7
  |G10| = 1  => [Q(y):Q] = 21
  |G9| = 3   => [Q(h):Q] = 7 with h=y^15 + 318834727/342275249*y^14 + -27336096046/342275249*y^13 + -25481561980/342275249*y^12 + 719247454576/342275249*y^11 + 671278748364/342275249*y^10 + -7561532813128/342275249*y^9 + -6280968784448/342275249*y^8 + 36010935682896/342275249*y^7 + 25834386777024/342275249*y^6 + -77534174991504/342275249*y^5 + -47586440684720/342275249*y^4 + 64509006038848/342275249*y^3 + 33834537544640/342275249*y^2 + -9520157954816/342275249*y
  |G1| = 21  => [Q:Q]    = 1
Tower 8
  |G10| = 1  => [Q(y):Q] = 21
  |G2| = 7   => [Q(a):Q] = 3 with a=y^19 + -72991/83766*y^18 + -3516995/41883*y^17 + 3034055/41883*y^16 + 101930956/41883*y^15 + -86248862/41883*y^14 + -1301411636/41883*y^13 + 1156131940/41883*y^12 + 8421000008/41883*y^11 + -2579914384/13961*y^10 + -29473382536/41883*y^9 + 26836469680/41883*y^8 + 56403794944/41883*y^7 + -47373532192/41883*y^6 + -18880507808/13961*y^5 + 38461198816/41883*y^4 + 8824608320/13961*y^3 + -10692074240/41883*y^2 + -3957553024/41883*y
  |G1| = 21  => [Q:Q]    = 1

 */
}
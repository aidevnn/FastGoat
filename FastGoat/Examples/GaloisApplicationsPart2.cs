using System.Numerics;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Perms;
using FastGoat.UserGroup.Polynoms;

namespace FastGoat.Examples;

public static class GaloisApplicationsPart2
{
    static void CheckChebotarev(KPoly<Rational> P, int nb, ConcreteGroup<Perm> gal, bool detail = false)
    {
        var lt = IntFactorisation.ChebotarevTypes(P, nb, detail)
            .GroupBy(e => e.Deconstruct())
            .ToDictionary(e => e.Key, e => e.Count());

        var types = gal.Select(perm => IntExt.PermutationToCycles(perm.Sn.N, perm.Table).Select(l => l.Length).Order().Deconstruct())
            .GroupBy(e => e).ToDictionary(e => e.Key, e => e.Count());
        types.Keys.OrderBy(e => e.ToString()).ToDictionary(e => e, e => types[e]).Println($"Expected types");

        var d = gal.Count();
        lt.Keys.OrderBy(e => e.ToString()).ToDictionary(e => e, e => (1.0 * d * lt[e]) / nb).Println("Actual types");

        Console.WriteLine();
    }

    static (KPoly<Rational> newP, KPoly<Rational> c) RewriteQuadraticPolynomial(KPoly<Rational> p, char c = 'a')
    {
        if (p.Degree > 2)
            throw new();

        var a = p[1] / 2;
        var p0 = p.Substitute(p.X - a).SubstituteChar('Y');
        var b = p0[0].Denom;
        var decomp = IntExt.PrimesDecompositionBigInt(b).GroupBy(i => i).ToDictionary(e => e.Key, e => e.Count());
        var b0 = decomp.Where(e => e.Value % 2 == 0).Aggregate(1, (prod, e) => prod * e.Key.Pow(e.Value / 2));
        var p1 = p.Substitute(p.X / b0 - a).SubstituteChar(c);
        Console.WriteLine($"{p} = 0 <=> {p1.Monic} = 0 with {c} = {(p.X + a) * b0}");
        return (p1.Monic, p.X / b0 - a);
    }

    static (KPoly<EPoly<Rational>> newP, KPoly<EPoly<Rational>> Y) RewriteQuadraticPolynomial2(KPoly<EPoly<Rational>> p, char c = 'a')
    {
        if (p.Degree > 2)
            throw new();

        var a = p[1] / 2;
        var p0 = p.Substitute(p.X - a).SubstituteChar(c);
        var b = IntExt.GcdBigInt(p0[0].Poly.Coefs.Select(e => e.Denom).ToArray());
        var decomp = IntExt.PrimesDecompositionBigInt(b).GroupBy(i => i).ToDictionary(e => e.Key, e => e.Count());
        var b0 = decomp.Where(e => e.Value >= 2)
            .Aggregate(1, (prod, e) => prod * e.Key.Pow(e.Value % 2 == 0 ? e.Value / 2 : (e.Value - 1) / 2));
        var p1 = p.Substitute(p.X / b0 - a).SubstituteChar(c);
        Console.WriteLine($"{p} = 0 <=> {p1.Monic} = 0 with {c} = {(p.X + a) * b0}");
        return (p1.Monic, p.X / b0 - a);
    }

    public static void ChebotarevExamples()
    {
        Ring.DisplayPolynomial = MonomDisplay.Caret;
        var x = FG.QPoly('X');

        {
            var P = x.Pow(3) - 3 * x - 1;
            var rootsK = IntFactorisation.AlgebraicRoots(P);
            var gal = GaloisTheory.GaloisGroup(rootsK, details: true);

            CheckChebotarev(P, 20, gal);
        }

        {
            var P = x.Pow(4) + x.Pow(3) + x.Pow(2) + x + 1;
            var rootsK = IntFactorisation.AlgebraicRoots(P);
            var gal = GaloisTheory.GaloisGroup(rootsK, details: true);

            CheckChebotarev(P, 30, gal);
        }

        {
            var P = x.Pow(5) + x.Pow(4) - 4 * x.Pow(3) - 3 * x.Pow(2) + 3 * x + 1;
            var rootsK = IntFactorisation.AlgebraicRoots(P);
            var gal = GaloisTheory.GaloisGroup(rootsK, details: true);

            CheckChebotarev(P, 20, gal);
        }

        {
            var P = x.Pow(6) + 243;
            var rootsK = IntFactorisation.AlgebraicRoots(P);
            var gal = GaloisTheory.GaloisGroup(rootsK, details: true);

            CheckChebotarev(P, 30, gal);
        }

        {
            var P = x.Pow(8) + 1;
            var rootsK = IntFactorisation.AlgebraicRoots(P);
            var gal = GaloisTheory.GaloisGroup(rootsK, details: true);

            CheckChebotarev(P, 30, gal);
        }

        {
            var P = x.Pow(4) + x + 1;
            CheckChebotarev(P, 60, new Symm(4), true); // slow
        }

        {
            var P = x.Pow(5) + 2;
            var s5 = new Sn(5);
            var gal = Group.Generate("C5 x: C4", s5, s5[(2, 3, 5, 4)], s5[(1, 2, 3, 4, 5)]);
            DisplayGroup.HeadElements(gal);
            CheckChebotarev(P, 60, gal, true); // slow
        }
    }

    public static void Example_Computing_Cos_2Pi_over_5()
    {
        Ring.DisplayPolynomial = MonomDisplay.Caret;
        var P = FG.CyclotomicPolynomial(5);
        var roots = new NthRootQ(5).Where(e => !e.Equals(e.One)).Order().ToList();

        var y = roots[0];
        var X = FG.KPoly('X', y);
        var cos = (y + 1 / y) / 2; // cos(2π/5)
        var isin = (y - 1 / y) / 2; // I * sin(2π/5)
        Console.WriteLine(new { cos });
        Console.WriteLine(new { isin });

        var croots = FG.NRoots(P.ToCPoly());
        var yc = croots[0];
        Console.WriteLine((yc, cos.Poly.Substitute(yc), isin.Poly.Substitute(yc)));
        croots.Select((c, k) => (c, (float)(2 * (k + 1)) / 5, ((float)c.Magnitude, (float)(c.Phase / Double.Pi))))
            .Println($"P = {croots.Aggregate(Cplx.X.One, (prod, r) => prod * (Cplx.X - r))}");

        var subFields = GaloisTheory.SubFields(roots).ToArray();
        var extTowers = GaloisApplications.ExtensionsTower(subFields);
        GaloisApplications.GaloisCorrespondence(extTowers);

        GaloisApplications.FindExtension(subFields, cos, $"Q(cos(2π/5))"); // Find containing subfield of cos(2π/5) 
        var P1 = RewriteQuadraticPolynomial(subFields[1].minPoly).newP;
        Console.WriteLine();

        var sqrt5 = IntFactorisation.AlgebraicRoots(P1.Substitute(X)).First(c => c.Poly.Substitute(yc).RealPart > 0);
        GaloisApplications.FindExtension(subFields, sqrt5, "Q(√5)");

        var sqrt5c = sqrt5.Poly.Substitute(yc);
        var rw = GaloisTheory.Rewrite(sqrt5, cos);
        Console.WriteLine($"{y} = {yc} with {y.F} = 0");
        Console.WriteLine($"a = √5 = {sqrt5} = {sqrt5c}");
        Console.WriteLine($"a^2 = {sqrt5.Pow(2)} = {sqrt5c.Pow(2)}");
        Console.WriteLine($"cos(2π/5) = {rw} = {rw.Substitute(sqrt5c)} = {Single.Cos(2 * Single.Pi / 5)}");
    }

    // Verifying expression from Daniel Guin book
    public static void Example_Computing_Cos_2Pi_over_17()
    {
        Ring.DisplayPolynomial = MonomDisplay.Caret;
        var n = 17;
        var P = FG.CyclotomicPolynomial(n).SubstituteChar('a');
        var roots = new NthRootQ(n).Where(e => !e.Equals(e.One)).Order().ToList();

        var y = roots[0];
        var cos = (y + 1 / y) / 2; // cos(2pi/n)
        var isin = (y - 1 / y) / 2; // I * sin(2pi/n)
        Console.WriteLine(new { cos });
        Console.WriteLine(new { isin });

        var croots = FG.NRoots(P.ToCPoly());
        var yc = croots[0];
        Console.WriteLine((yc, cos.Poly.Substitute(yc), isin.Poly.Substitute(yc)));
        croots.Select((c, k) => (c, (float)(2 * (k + 1)) / n, ((float)c.Magnitude, (float)(c.Phase / Double.Pi))))
            .Println($"P = {croots.Aggregate(Cplx.X.One, (prod, r) => prod * (Cplx.X - r))}");

        var subFields = GaloisTheory.SubFields(roots).ToArray();
        var extTowers = GaloisApplications.ExtensionsTower(subFields);
        GaloisApplications.GaloisCorrespondence(extTowers);

        GaloisApplications.FindExtension(subFields, cos, $"Q(cos(2π/{n}))"); // Find containing subfield of cos(2π/17)

        var sf3 = subFields[3];
        var (X3, y3) = FG.EPolyXc(sf3.minPoly, 'c');
        var yf = FG.NRoots(sf3.minPoly.ToCPoly())[0];
        var roots2 = IntFactorisation.AlgebraicRoots(sf3.minPoly.Substitute(X3), true);
        var minPoly = IntFactorisation.MinPolynomial(roots2.Aggregate(X3.One, (p0, r) => p0 * (X3 - r / 2)));
        Console.WriteLine(minPoly);

        Console.WriteLine(new { cos });
        Console.WriteLine(cos.Poly.Substitute(yc));
        Console.WriteLine(minPoly.Substitute(cos));

        var (cosf0, _) = roots2.Select(r => (r / 2, (r / 2).Poly.Substitute(yf).K))
            .Where(e => Complex.IsRealNumber(e.K) && e.K.Real > 0)
            .MaxBy(e => e.K.Real);

        Console.WriteLine(new { cos_actual = cosf0 });
        Console.WriteLine(cosf0.Poly.Substitute(yf));
        Console.WriteLine(minPoly.Substitute(cosf0));

        // Daniel Guin book, page 438
        // cos(2π/17) = 1/16 * (√17 - 1 + √2 * (√(17 - √17) + √(34 + 6 * √17 + √2 * (√17 - 1) * √(17 - √17) - 8 * √2 * √(17 + √17)))
        // cos(2π/17) = 1/16 * (√17 - 1 + √(34 - 2*√17) + √(68 + 12 * √17 + 2 * (√17 - 1) * √(34 - 2 * √17) - 16 * √(34 + 2 * √17))

        var sqrt_17 = IntFactorisation.AlgebraicRoots(X3.Pow(2) - 17)
            .First(e => e.Poly.Substitute(yf).IsReal() && e.Poly.Substitute(yf).RealPart > 0);

        var alpha = IntFactorisation.AlgebraicRoots(X3.Pow(2) - (34 - 2 * sqrt_17))
            .First(e => e.Poly.Substitute(yf).IsReal() && e.Poly.Substitute(yf).RealPart > 0);

        var beta = IntFactorisation.AlgebraicRoots(X3.Pow(2) - (34 + 2 * sqrt_17))
            .First(e => e.Poly.Substitute(yf).IsReal() && e.Poly.Substitute(yf).RealPart > 0);

        var c0 = 68 + 12 * sqrt_17 + 2 * (sqrt_17 - 1) * alpha - 16 * beta;
        var gamma = IntFactorisation.AlgebraicRoots(X3.Pow(2) - c0)
            .First(e => e.Poly.Substitute(yf).IsReal() && e.Poly.Substitute(yf).RealPart > 0);

        var cosf1 = (sqrt_17 - 1 + alpha + gamma) / 16;
        Console.WriteLine(new { cos_expected = cosf1 });
        Console.WriteLine(cosf1.Poly.Substitute(yf));
        Console.WriteLine(minPoly.Substitute(cosf1));
    }

    // Scraping the tower of subfields extensions
    public static void Example2_Computing_Cos_2Pi_over_17()
    {
        Ring.DisplayPolynomial = MonomDisplay.Caret;
        var n = 17;
        var P = FG.CyclotomicPolynomial(n).SubstituteChar('a');
        var roots = new NthRootQ(n).Where(e => !e.Equals(e.One)).Order().ToList();

        // Stage 0 : the Floor
        var y = roots[0];
        var X = FG.KPoly('X', y);
        var cos = (y + 1 / y) / 2; // cos(2pi/n)
        var isin = (y - 1 / y) / 2; // I * sin(2pi/n)
        Console.WriteLine(new { cos });
        Console.WriteLine(new { isin });

        var croots = FG.NRoots(P.ToCPoly());
        var yc = croots[0];
        Console.WriteLine((yc, cos.Poly.Substitute(yc), isin.Poly.Substitute(yc)));
        croots.Select((c, k) => (c, (float)(2 * (k + 1)) / n, ((float)c.Magnitude, (float)(c.Phase / Double.Pi))))
            .Println($"P = {croots.Aggregate(Cplx.X.One, (prod, r) => prod * (Cplx.X - r))}");

        var subFields = GaloisTheory.SubFields(roots).ToArray();
        var extTowers = GaloisApplications.ExtensionsTower(subFields);
        GaloisApplications.GaloisCorrespondence(extTowers);

        GaloisApplications.FindExtension(subFields, cos, $"Q(cos(2π/{n}))"); // Find containing subfield of cos(2π/17)

        var (_, sf1, sf2, sf3, sf4) = subFields.Deconstruct();
        Console.WriteLine(new { sf1 = "subField1", subGr = sf1.SubGr.ShortName, dim = sf1.roots.Length, sf1.minPoly, sf1.primElt });
        var P1 = RewriteQuadraticPolynomial(sf1.minPoly).newP;
        Console.WriteLine();

        var sqrt17 = IntFactorisation.AlgebraicRoots(P1.Substitute(X)).First(c => c.Poly.Substitute(yc).RealPart > 0);
        var sqrt17c = sqrt17.Poly.Substitute(yc);
        GaloisApplications.FindExtension(subFields, sqrt17, "Q(√17)");

        Console.WriteLine(new { sf2 = "subField2", subGr = sf2.SubGr.ShortName, dim = sf2.roots.Length, sf2.minPoly, sf2.primElt });

        // Stage 1
        var (X1, y1) = FG.EPolyXc(P1, 'a');
        var (f1, newP1, subs1, newPc1) = IntFactorisation.AlgebraicFactors(sf2.minPoly.Substitute(X1), true)
            .Select(f => (f, RewriteQuadraticPolynomial2(f, 'Y')))
            .Select(e => (e.f, e.Item2.newP, e.Item2.Y, e.Item2.newP.ToCPoly(sqrt17c)))
            .Where(e => e.Item4[0].RealPart < 0)
            .MinBy(e => -e.Item4[0]);

        var fc1 = f1.ToCPoly(sqrt17c);
        var subsc1 = subs1.ToCPoly(sqrt17c);

        Console.WriteLine(new { f1, fc1 });
        Console.WriteLine(new { newP1, subs1, newPc1, subsc1 });
        var e1 = FG.NRoots(newPc1).First(e => e.IsReal() && e.RealPart > 0);
        var xe1 = subsc1.Substitute(e1);
        Console.WriteLine(new { e = e1, e2 = e1.Pow(2), xe = xe1 });

        Console.WriteLine(newPc1.Substitute(e1));
        Console.WriteLine(fc1.Substitute(xe1));

        Console.WriteLine(new { sf3 = "subField3", dim = sf3.roots.Length, sf3.minPoly, sf3.primElt });

        // Stage 2
        var (r2, a2, b2) = IntFactorisation.PrimitiveElt(newP1);
        var (X2, y2) = FG.EPolyXc(r2, 'b');
        var sqrt17_1 = a2.Substitute(y2);
        var s2 = b2.Substitute(y2);
        Console.WriteLine((sqrt17_1, sqrt17_1.Pow(2)));
        Console.WriteLine((s2, s2.Pow(2), (17 - sqrt17_1) / 2));

        var (f2, newP2, subs2, newPc2) = IntFactorisation.AlgebraicFactors(sf3.minPoly.Substitute(X2), true)
            .Select(f => (f, RewriteQuadraticPolynomial2(f, 'Y')))
            .Select(e => (e.f, e.Item2.newP, e.Item2.Y, e.Item2.newP.ToCPoly(e1)))
            .Where(e => e.Item4[0].RealPart < 0)
            .MinBy(e => -e.Item4[0]);

        Console.WriteLine(new { f2, newP2, subs2, newPc2 });
        var fc2 = f2.ToCPoly(e1);
        var subsc2 = subs2.ToCPoly(e1);

        Console.WriteLine(new { f2, fc2 });
        Console.WriteLine(new { newP2, subs2, newPc2, subsc2 });
        var e2 = FG.NRoots(newPc2).First(e => e.IsReal() && e.RealPart > 0);
        var xe2 = subsc2.Substitute(e2);
        Console.WriteLine(new { e = e2, e2 = e2.Pow(2), xe = xe2 });

        Console.WriteLine(newPc2.Substitute(e2));
        Console.WriteLine(fc2.Substitute(xe2));

        // Stage 3
        var (r3, a3, b3) = IntFactorisation.PrimitiveElt(newP2);
        Console.WriteLine(new { r3, a3, b3 });
        var (X3, y3) = FG.EPolyXc(r3, 'c');
        var s2_1 = a3.Substitute(y3);
        var sqrt17_2 = sqrt17_1.Substitute(s2_1);
        var s3 = b3.Substitute(y3);
        Console.WriteLine((sqrt17_2, sqrt17_2.Pow(2)));
        Console.WriteLine((s2_1, s2_1.Pow(2), (17 - sqrt17_2) / 2));
        Console.WriteLine((s3, s3.Pow(2)));
        Console.WriteLine($"d = {s2_1} s3^2 = {GaloisTheory.Rewrite(s2_1, s3.Pow(2)).SubstituteChar('d')}");

        // Final Boss
        var mP3 = sf3.minPoly;
        var roots2 = IntFactorisation.AlgebraicRoots(mP3.Substitute(X3), true);

        var minPoly = IntFactorisation.MinPolynomial(roots2.Aggregate(X3.One, (p0, r) => p0 * (X3 - r / 2)));
        Console.WriteLine(minPoly);
        Console.WriteLine(mP3.Substitute(mP3.X * 2).Monic);

        var yf = FG.NRoots(r3.ToCPoly())[0];
        var (cosf, e3) = roots2.Select(r => (r / 2, (r / 2).Poly.Substitute(yf).K))
            .Where(e => Complex.IsRealNumber(e.K) && e.K.Real > 0)
            .MaxBy(e => e.K.Real);

        Console.WriteLine(new { cos });
        Console.WriteLine(cos.Poly.Substitute(yc));
        Console.WriteLine(minPoly.Substitute(cos));
        Console.WriteLine(new { cosf });
        Console.WriteLine(cosf.Poly.Substitute(yf));
        Console.WriteLine(minPoly.Substitute(cosf));
    }
}
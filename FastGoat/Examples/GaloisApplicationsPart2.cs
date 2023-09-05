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
    static IEnumerable<int[]> FpFactorsTypes(KPoly<Rational> f, bool details = false)
    {
        var discQ = BigInteger.Abs(Ring.Discriminant(f).Num);
        var discDecomp = IntExt.PrimesDec(discQ);
        if (details)
        {
            Console.WriteLine($"f = {f}");
            Console.WriteLine($"Disc(f) = {discQ} ~ {discDecomp.AscendingByKey().GlueMap(" * ", "{0}^{1}")}");
        }

        var k = 1;
        foreach (var p in IntExt.Primes10000)
        {
            if (BigInteger.Remainder(discQ, p) == 0)
                continue;

            var listIrr0 = IntFactorisation.FirrFp(f, p);
            var type = listIrr0.Select(g => g.Degree).Order().ToArray();
            if (details)
            {
                Console.WriteLine("#{2,-3} P = {0} shape ({1})", p, type.Glue(", "), k++);
            }

            yield return type;
        }
    }
    public static void CheckChebotarev(KPoly<Rational> P, ConcreteGroup<Perm> gal, bool detail = false)
    {
        var types = gal.Select(perm => IntExt.PermutationToCycles(perm.Sn.N, perm.Table).Select(l => l.Length).Order().Deconstruct())
            .GroupBy(e => e).ToDictionary(e => e.Key, e => e.Count());
        var d = gal.Count();

        var setTypes = types.Keys.Select(e => e.ToArray()).ToHashSet(new SequenceEquality<int>());
        var allTypes = new List<int[]>();
        foreach (var fpType in FpFactorsTypes(P, detail))
        {
            allTypes.Add(fpType);
            setTypes.Remove(fpType);
            if (setTypes.Count == 0)
                break;
        }

        var lt = allTypes.GroupBy(e => e.Deconstruct()).ToDictionary(e => e.Key, e => e.Count());

        var nb = lt.Sum(e => e.Value);
        types.Keys.OrderBy(e => e.ToString()).ToDictionary(e => e, e => types[e]).Println($"Expected types");
        lt.Keys.OrderBy(e => e.ToString()).ToDictionary(e => e, e => Double.Round((1.0 * d * lt[e]) / nb, 1)).Println("Actual types");

        Console.WriteLine("######################################");
        Console.WriteLine();
    }

    public static Dictionary<Array2Tuple<int>, int> GroupTypes(ConcreteGroup<Perm> gal)
    {
        var res = gal.Select(perm => IntExt.PermutationToCycles(perm.Sn.N, perm.Table).Select(l => l.Length).Order().Deconstruct())
            .GroupBy(e => e).ToDictionary(e => e.Key, e => e.Count());

        return res;
    }

    static double DistanceTypes(Dictionary<Array2Tuple<int>, int> actual, Dictionary<Array2Tuple<int>, int> expected)
    {
        var d = expected.Sum(e => e.Value);
        var nb = actual.Sum(e => e.Value);
        return expected.Select(e => actual.TryGetValue(e.Key, out var ae) ? Double.Abs((1.0 * d * ae) / nb - e.Value) : e.Value)
            .Average();
    }

    public static ConcreteGroup<Perm> GaloisGroupChebotarev(KPoly<Rational> P, bool details = false)
    {
        var deg = P.Degree;
        if (deg > 7)
            throw new("Deg(P) greater than 7 isnt supported");

        if (IntFactorisation.FirrZ2(P).Length != 1)
            throw new($"P = {P} isnt irreductible");

        var sym = new Symm(deg);
        var transSubGr = sym.TransitiveSubGroups().Select(e => Group.Generate(e.name, sym, e.gens))
            .Select(g => (g, GroupTypes(g))).Select(e => (e.Item2.Keys.OrderBy(a => a.ToString()).ToArray(), (e.Item2, e.g)))
            .ToArray();
        var transTypes = transSubGr.Select(e => e.Item1).ToArray();

        var types = new Dictionary<Array2Tuple<int>, int>();
        foreach (var ct in FpFactorsTypes(P, details: details))
        {
            var ctype = ct.Order().Deconstruct();
            if (!types.ContainsKey(ctype))
                types[ctype] = 1;
            else
                types[ctype] = types[ctype] + 1;

            if (types.Sum(e => e.Value) < 10)
                continue;

            var sel1 = transTypes.Where(tp => tp.ToHashSet().IsSupersetOf(types.Keys)).ToArray();
            var sel2 = sel1.Select(e => transSubGr.First(f => f.Item1.Equals(e)).Item2)
                .Select(e => (e, DistanceTypes(types, e.Item1))).OrderBy(e => e.Item2).ToArray();

            if (sel1.Length == 1 || sel2[1].Item2 / sel2[0].Item2 > 5) // testing precision
            {
                var g = sel2[0].e.g;
                if (details)
                {
                    types.OrderByDescending(e => e.Key.ToString()).Println("actual types");
                    sel2[0].e.Item1.OrderByDescending(e => e.Key.ToString()).Println("expected types");
                    sel2.Select(e => new { e.e.g.Name, order = e.e.g.Count(), dist = e.Item2 }).Println("Distances");
                    Console.WriteLine();
                
                    Console.WriteLine($"P = {P}");
                    Console.WriteLine($"Gal(P) = {g.Name}");
                    DisplayGroup.Head(g);
                }

                return g;
            }
        }

        throw new();
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
        Ring.DisplayPolynomial = MonomDisplay.StarCaret;
        var x = FG.QPoly();

        {
            var P = x.Pow(3) - 3 * x - 1;
            var rootsK = IntFactorisation.AlgebraicRoots(P);
            var gal = GaloisTheory.GaloisGroup(rootsK, details: true);

            CheckChebotarev(P, gal, detail: true);
        }
        
        {
            var P = x.Pow(4) + x.Pow(3) + x.Pow(2) + x + 1;
            var rootsK = IntFactorisation.AlgebraicRoots(P);
            var gal = GaloisTheory.GaloisGroup(rootsK, details: true);

            CheckChebotarev(P, gal, detail: true);
        }

        {
            var P = x.Pow(5) + x.Pow(4) - 4 * x.Pow(3) - 3 * x.Pow(2) + 3 * x + 1;
            var rootsK = IntFactorisation.AlgebraicRoots(P);
            var gal = GaloisTheory.GaloisGroup(rootsK, details: true);

            CheckChebotarev(P, gal, detail: true);
        }

        {
            var P = x.Pow(6) + 243;
            var rootsK = IntFactorisation.AlgebraicRoots(P);
            var gal = GaloisTheory.GaloisGroup(rootsK, details: true);

            CheckChebotarev(P, gal, detail: true);
        }

        {
            var P = x.Pow(8) + 1;
            var rootsK = IntFactorisation.AlgebraicRoots(P);
            var gal = GaloisTheory.GaloisGroup(rootsK, details: true);

            CheckChebotarev(P, gal, detail: true);
        }

        {
            var P = x.Pow(4) + x + 1;
            var gal = new Symm(4);
            DisplayGroup.HeadElements(gal);
            CheckChebotarev(P, gal, detail: true);
        }

        {
            var P = x.Pow(5) + 2;
            var s5 = new Sn(5);
            var gal = Group.Generate("C5 x: C4", s5, s5[(2, 3, 5, 4)], s5[(1, 2, 3, 4, 5)]);
            DisplayGroup.HeadElements(gal);
            CheckChebotarev(P, gal, detail: true);
        }

        {
            var P = x.Pow(7) - 14 * x.Pow(5) + 56 * x.Pow(3) - 56 * x + 22;
            var s7 = new Sn(7);
            var gal = Group.Generate("C7 x: C3", s7, s7[(1, 2, 3, 4, 5, 6, 7)], s7[(2, 3, 5), (4, 7, 6)]);
            DisplayGroup.HeadElements(gal);
            CheckChebotarev(P, gal, detail: true);
        }

        {
            var P = x.Pow(6) + 2 * x.Pow(5) + 3 * x.Pow(4) + 4 * x.Pow(3) + 5 * x.Pow(2) + 6 * x + 7;
            var s6 = new Sn(6);
            var gal = Group.Generate("PGL(2, 5)", s6, s6[(3, 6, 5, 4)], s6[(1, 2, 5), (3, 4, 6)]);
            DisplayGroup.Head(gal);
            CheckChebotarev(P, gal, detail: true);
        }
    }

    public static void GaloisGroupTransitivityExamples()
    {
        Ring.DisplayPolynomial = MonomDisplay.StarCaret;
        var x = FG.QPoly();
        GaloisGroupChebotarev(x.Pow(3) - 3 * x - 1, details: true);
        GaloisGroupChebotarev(x.Pow(4) + x.Pow(3) + x.Pow(2) + x + 1, details: true);
        GaloisGroupChebotarev(x.Pow(5) + x.Pow(4) - 4 * x.Pow(3) - 3 * x.Pow(2) + 3 * x + 1, details: true);
        GaloisGroupChebotarev(x.Pow(6) + 243, details: true);
        GaloisGroupChebotarev(x.Pow(4) + x + 1, details: true);
        GaloisGroupChebotarev(x.Pow(5) + 2, details: true);
        GaloisGroupChebotarev(x.Pow(7) - 14 * x.Pow(5) + 56 * x.Pow(3) - 56 * x + 22, details: true);

        // J.S. Milne exercice 4.6 page 60 without using the “polgalois” command in PARI
        GaloisGroupChebotarev(x.Pow(6) + 2 * x.Pow(5) + 3 * x.Pow(4) + 4 * x.Pow(3) + 5 * x.Pow(2) + 6 * x + 7, details: true);
    }
    /***
        f = x^6 + 2*x^5 + 3*x^4 + 4*x^3 + 5*x^2 + 6*x + 7
        Disc(f) = 157351936 ~ 2^16 * 7^4
        #1   P = 3 shape [6]
        #2   P = 5 shape [3, 3]
        #3   P = 11 shape [6]
        #4   P = 13 shape [1, 5]
        #5   P = 17 shape [1, 5]
        #6   P = 19 shape [1, 1, 4]
        #7   P = 23 shape [1, 1, 4]
        #8   P = 29 shape [1, 5]
        #9   P = 31 shape [1, 1, 4]
        #10  P = 37 shape [3, 3]
        types
            [(6), 2]
            [(3, 3), 2]
            [(1, 5), 3]
            [(1, 1, 4), 3]
        Distances
            { Name = L(6):2 = PGL(2,5) = S_5(6), dist = 7.428571428571429 }
            { Name = S6, dist = 59.27272727272727 }
        P = x^6 + 2*x^5 + 3*x^4 + 4*x^3 + 5*x^2 + 6*x + 7
        |L(6):2 = PGL(2,5) = S_5(6)| = 120
     */

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
        var rw = IntFactorisation.Rewrite(sqrt5, cos);
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
        var minPoly = IntFactorisation.GetBaseAndMinPolynomial(cos).Item2;
        var roots2 = IntFactorisation.AlgebraicRoots(minPoly.Substitute(X3 / 2).Monic, true).Select(r => r / 2).ToList();
        Console.WriteLine(minPoly);

        Console.WriteLine(new { cos });
        Console.WriteLine(cos.Poly.Substitute(yc));
        Console.WriteLine(minPoly.Substitute(cos));

        var (cosf0, _) = roots2.Select(r => (r, r.Poly.Substitute(yf).K))
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

        var croots = FG.NRoots(P.ToCPoly()).Order().ToArray();
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
        Console.WriteLine("({0})^2 =  {1} ", sqrt17, sqrt17.Pow(2));
        Console.WriteLine("{0}^2 = {1} ", sqrt17c, sqrt17c.Pow(2));
        Console.WriteLine();
        
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
        Console.WriteLine($"d = {s2_1} s3^2 = {IntFactorisation.Rewrite(s2_1, s3.Pow(2)).SubstituteChar('d')}");
        
        // Final Boss
        var minPoly = IntFactorisation.GetBaseAndMinPolynomial(cos).Item2;
        var roots2 = IntFactorisation.AlgebraicRoots(minPoly.Substitute(X3 / 2).Monic, true).Select(r => r / 2).ToList();;
        Console.WriteLine(minPoly);
        Console.WriteLine(sf3.minPoly);

        var yf = FG.NRoots(r3.ToCPoly())[0];
        var (cosf, e3) = roots2.Select(r => (r, (r).Poly.Substitute(yf).K))
            .Where(e => Complex.IsRealNumber(e.K) && e.K.Real > 0)
            .MaxBy(e => e.K.Real);

        Console.WriteLine(new { cos });
        Console.WriteLine(cos.Poly.Substitute(yc));
        Console.WriteLine(minPoly.Substitute(cos));
        Console.WriteLine(new { cosf });
        Console.WriteLine(cosf.Poly.Substitute(yf));
        Console.WriteLine(minPoly.Substitute(cosf));
    }
    
    static void IsConstructible(int n, bool fast = true)
    {
        var a = Cnf.Nth(2 * n); // e(i*π/n)
        Console.WriteLine($"a = e(i*π/{n}) and {a.E.F} = 0");
        var a0 = FG.NRoots(a.E.F.ToBcPoly()).Where(e => e.ImaginaryPart.CompareTo(e.Zero.RealPart) == 1).MaxBy(e => e.RealPart);
        Console.WriteLine("[ cos(π/{0}); sin(π/{0})] = {1}", n, a0.ToFixForm());
        Console.WriteLine("cos(π/{0}) = {1}", n, a.Re.E);
        Console.WriteLine();

        var (_, P0) = IntFactorisation.GetBaseAndMinPolynomial(a.Re.E, 'x');
        Console.WriteLine($"[Q( cos(π/{n}) ) / Q] = {P0.Degree}");
        Console.WriteLine($"    MinPolynomial P = {P0}");
        Console.WriteLine($"    and P( cos(π/{n}) ) = 0");

        if (fast)
        {
            var deg = P0.Degree;
            if((new BigInteger(deg).IsPowerOfTwo))
                Console.WriteLine($"########### cos(π/{n}) CONSTRUCTIBLE");
            else
                Console.WriteLine($"########### cos(π/{n}) NOT CONSTRUCTIBLE");
    
            Console.WriteLine();
            return;
        }

        var P1 = P0.Substitute(P0.X / 2).Monic;
        var roots = IntFactorisation.SplittingField(P1).ToList();
        var subFields = GaloisTheory.SubFields(roots).ToArray();
        var extTowers = GaloisApplications.ExtensionsTower(subFields);
        GaloisApplications.GaloisCorrespondence(extTowers);
        if (extTowers.Any(e => e.Select((c0, i) => (c0, i)).All(f => 2.Pow(f.i) == f.c0.SubGr.Count())))
            Console.WriteLine($"########### cos(π/{n}) CONSTRUCTIBLE");
        else
            Console.WriteLine($"########### cos(π/{n}) NOT CONSTRUCTIBLE");
    
        Console.WriteLine();
    }

    public static void ConstructibilityFast()
    {
        // Fast checking Constructibility of cos(π/n) for n from 3 to 33
        for (int n = 3; n < 34; ++n)
            IsConstructible(n);
    }
    
    public static void ConstructibilityDetails()
    {
        // Detailed Galois Group and Constructibility of cos(π/n) for n from 5 to 20
        for (int n = 5; n < 20; ++n)
            IsConstructible(n, fast: false);
    }

    public static void GaloisGroupPolynomialDegree7()
    {
        Ring.DisplayPolynomial = MonomDisplay.StarCaret;
        var j = 1;
        for (int i = 1; i < 8; i++)
        {
            Console.WriteLine($"Degree {i}");
            var si = new Symm(i);
            var subGr = si.TransitiveSubGroups().Select(e => Group.Generate(e.name, si, e.gens));
            var isos = new Queue<ConcreteGroup<Perm>>(subGr);
            foreach (var P in FG.GaloisGroupPolynomialsList(i))
            {
                var iso = isos.Dequeue();
                var galP = GaloisGroupChebotarev(P);
                Console.WriteLine($"P{j,-2} = {P}");
                if (!iso.IsIsomorphicTo(galP))
                    throw new();
            
                var d = Ring.Discriminant(P).IsSquare ? "\"+\"" : "\"-\"";
                Console.WriteLine($"    Gal(P{j++}) = [{galP.Name}, {d}, {galP.Count()}, {iso.GetGenerators().Glue(" ")}]");
            }

            if (isos.Count != 0)
                throw new();
        
            Console.WriteLine();
        }
    }
}
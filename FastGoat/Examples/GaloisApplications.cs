using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Perms;
using FastGoat.UserGroup.Polynoms;

namespace FastGoat.Examples;

public static class GaloisApplications
{
    static List<ConcreteGroup<Perm>> MaxSubGr(List<ConcreteGroup<Perm>> allSubGr, ConcreteGroup<Perm> g)
    {
        var allMax = new List<ConcreteGroup<Perm>>();
        foreach (var h in allSubGr)
        {
            if (g.Count() <= h.Count() || !h.SubSetOf(g))
                continue;

            allMax = allMax.Except(allMax.Where(h0 => h0.SubSetOf(h))).ToList();
            if (allMax.All(h0 => !h.SubSetOf(h0)))
                allMax.Add(h);
        }

        return allMax;
    }

    static IEnumerable<List<ConcreteGroup<Perm>>> SubGrLattice(List<ConcreteGroup<Perm>> allSubGr)
    {
        var g0 = allSubGr.MaxBy(g => g.Count());
        if (g0 is not null)
        {
            var all = new List<List<ConcreteGroup<Perm>>>() { new() { g0 } };
            while (all.Count != 0)
            {
                var tmp = all.ToList();
                all.Clear();
                foreach (var lt in tmp)
                {
                    var g = lt.Last();
                    if (g.Count() == 1)
                    {
                        yield return lt;
                        continue;
                    }

                    foreach (var h in MaxSubGr(allSubGr, g))
                    {
                        var lt2 = lt.Append(h).ToList();
                        all.Add(lt2);
                    }
                }
            }
        }
    }

    public static GaloisCorrespondence[][] ExtensionsTower(GaloisCorrespondence[] subFields)
    {
        var gr2field = subFields.ToDictionary(e => e.SubGr, e => e);
        return SubGrLattice(gr2field.Keys.ToList()).Select(lt =>
            lt.Reverse<ConcreteGroup<Perm>>().Select(sg => gr2field[sg]).ToArray()
        ).ToArray();
    }

    public static void FindExtension(GaloisCorrespondence[] galCor, EPoly<Rational> a, string name = "Q(r)")
    {
        var sel = galCor.Select(e => (e, pol: IntFactorisation.Rewrite(e.primElt, a)))
            .Where(e => !e.pol.IsZero()).OrderBy(e => e.e.SubGr.Count()).ToArray();

        if (sel.Length == 0)
        {
            Console.WriteLine("Not Found");
            Console.WriteLine();
        }
        else
        {
            var minPol = IntFactorisation.GetBaseAndMinPolynomial(a).Item2;
            var rows = sel.Select(e => new[]
                    { e.e.SubGr.ShortName, $"[{e.e.FieldName}:Q] = {e.e.minPoly.Degree}", $"z = {e.pol.SubstituteChar(e.e.primEltName[0])}" })
                .ToArray();

            var sf = sel.Last().e;
            Console.WriteLine($"{name} = Q(z) is subfield of {sf.FieldName} and [{name}:Q] = {sf.minPoly.Degree}");
            Console.WriteLine($"With {minPol.SubstituteChar('z')} = 0");
            Console.WriteLine($"With {a.F} = 0");
            Ring.DisplayMatrix(Ring.Matrix(rows.Length, rows.SelectMany(e => e).ToArray()), Ring.MatrixDisplay.TableLeft, "  \t");
            Console.WriteLine();
        }
    }

    public static GaloisCorrespondence[] GaloisCorrespondence(KPoly<Rational> P, int nbGens = 2)
    {
        var subFields = GaloisTheory.SubFields(P, nbGens).ToArray();
        var extTowers = ExtensionsTower(subFields);
        GaloisCorrespondence(extTowers);
        return subFields;
    }
    
    public static GaloisCorrespondence[] GaloisCorrespondence(List<EPoly<Rational>> roots, int nbGens = 2)
    {
        var subFields = GaloisTheory.SubFields(roots, nbGens).ToArray();
        var extTowers = ExtensionsTower(subFields);
        GaloisCorrespondence(extTowers);
        return subFields;
    }

    public static GaloisCorrespondence[] GaloisCorrespondence(List<KAut<Rational>> roots, int nbGens = 2)
    {
        return GaloisCorrespondence(roots.Select(r => r.E).ToList(), nbGens);
    }

    public static void GaloisCorrespondence(GaloisCorrespondence[][] exts)
    {
        var galCor = exts[0].Last();
        if (galCor is null)
            throw new();

        var P = galCor.primElt.F;
        
        var (X, y) = FG.EPolyXc(galCor.primElt.F, galCor.primElt.F.x);
        Console.WriteLine($"With P = {P.Substitute(X)}");
        galCor.roots.Order().Select(r => X - r).Println($"Factorization in Q({y})/Q");
        Console.WriteLine();
        Console.WriteLine("Galois Group");
        DisplayGroup.HeadElements(galCor.SubGr);
        var i = 1;
        foreach (var lt in exts.OrderBy(t => t.Select(gc => gc.SubGr.Count()).ToArray(),
                     Comparer<int[]>.Create((l1, l2) => l1.SequenceCompareTo(l2)))
                     .ThenBy(t => t.Select(gc => gc.primEltName).ToArray(),
                         Comparer<string[]>.Create((l1, l2) => l1.SequenceCompareTo(l2))))
        {
            Console.WriteLine($"Tower {i++}");
            foreach (var gc in lt)
            {
                var name = $"[{gc.FieldName}:Q]";
                var details = gc.primEltName == "1" || gc.primEltName == $"{y}"
                    ? ""
                    : $" with {gc.primEltName} = {gc.primElt}";
                Console.WriteLine($"  {gc.SubGr.ShortName,-10} => {name,-8} = {gc.minPoly.Degree}{details}");
            }
        }

        Console.WriteLine();
    }

    public static void Example1()
    {
        Ring.DisplayPolynomial = MonomDisplay.Caret;
        var x = FG.QPoly('X');
        var (minPoly, X, sqrt2, sqrt3) = IntFactorisation.PrimitiveElt(x.Pow(2) - 2, x.Pow(2) - 3).First();
        Console.WriteLine("Q(√2, √3) = Q(α)");
        Console.WriteLine($"MinPoly P = {minPoly} with P(α) = 0");

        var subFields = GaloisTheory.SubFields(minPoly).ToArray();
        var extTowers = ExtensionsTower(subFields);
        GaloisCorrespondence(extTowers);

        FindExtension(subFields, sqrt2.One, "Q");

        FindExtension(subFields, sqrt2, "Q(√2)");
        FindExtension(subFields, sqrt3, "Q(√3)");
        FindExtension(subFields, sqrt2 * sqrt3, "Q(√6)");
        FindExtension(subFields, sqrt2 + sqrt3, "Q(√2 + √3)");
    }

    public static void Example2()
    {
        Ring.DisplayPolynomial = MonomDisplay.Caret;
        var x = FG.QPoly('X');
        var (minPoly, X, j, cbrt2) = IntFactorisation.PrimitiveElt(x.Pow(2) + x + 1, x.Pow(3) - 2).First();
        Console.WriteLine("Q(j, ∛2) with j^2+j+1=0");
        Console.WriteLine($"MinPoly P = {minPoly} with P(α) = 0");

        var subFields = GaloisTheory.SubFields(minPoly).ToArray();
        var extTowers = ExtensionsTower(subFields);
        GaloisCorrespondence(extTowers);

        FindExtension(subFields, j.One, "Q");

        FindExtension(subFields, j, "Q(j)");
        FindExtension(subFields, cbrt2, "Q(∛2)");
        FindExtension(subFields, j * cbrt2, "Q(j * ∛2)");
        FindExtension(subFields, j.Pow(2) * cbrt2, "Q(j^2 * ∛2)");

        FindExtension(subFields, j.X, "Q(j, ∛2)");
    }

    public static void Example3()
    {
        Ring.DisplayPolynomial = MonomDisplay.Caret;
        var x = FG.QPoly('X');
        var (minPoly, X, i, a) = IntFactorisation.PrimitiveElt(x.Pow(2) + 1, x.Pow(4) - 2).First();
        Console.WriteLine("Q(i, α) with α^4 - 2 = 0");

        var subFields = GaloisTheory.SubFields(minPoly).ToArray();
        var extTowers = ExtensionsTower(subFields);
        GaloisCorrespondence(extTowers);

        var sqrt2 = IntFactorisation.AlgebraicRoots(X.Pow(2) - 2).Last();
        var primElt_i_sqrt2 = GaloisTheory.PrimitiveEltComb(i, sqrt2).W;

        // Daniel Guin book example at pages 350-353
        FindExtension(subFields, a.One, "Q");

        FindExtension(subFields, i, "Q(i)");
        FindExtension(subFields, sqrt2, "Q(√2)");
        FindExtension(subFields, sqrt2 * i, "Q(i√2)");

        FindExtension(subFields, primElt_i_sqrt2, "Q(i, √2)");
        FindExtension(subFields, a, "Q(a)");
        FindExtension(subFields, i * a, "Q(i*a)");
        FindExtension(subFields, (1 + i) * a, "Q((1+i)a)");
        FindExtension(subFields, (1 - i) * a, "Q((1-i)a)");

        FindExtension(subFields, i + a, "Q(i + a)");
    }

    // Quaternion Galois Group With assistance of the New Bing 
    public static void Example4()
    {
        Ring.DisplayPolynomial = MonomDisplay.StarCaret;
        var x = FG.QPoly('X');
        var (minPoly1, X, sqrt2, sqrt3) = IntFactorisation.PrimitiveElt(x.Pow(2) - 2, x.Pow(2) - 3).First();
        var (minPoly2, _, e0) = IntFactorisation.PrimitiveElt(X.Pow(2) - ((2 + sqrt2) * (3 + sqrt3)));

        var subFields = GaloisTheory.SubFields(minPoly2).ToArray();
        var extTowers = ExtensionsTower(subFields);
        GaloisCorrespondence(extTowers);

        var X0 = FG.KPoly('X', subFields.First().primElt.X);
        var (a, _) = IntFactorisation.AlgebraicRoots(X0.Pow(2) - 2).Deconstruct();
        var (b, _) = IntFactorisation.AlgebraicRoots(X0.Pow(2) - 3).Deconstruct();
        var e = e0.Substitute(a.X);

        FindExtension(subFields, a.One, "Q");
        FindExtension(subFields, a + b, "Q(α+β)");
        FindExtension(subFields, a, "Q(α)");
        FindExtension(subFields, b, "Q(β)");
        FindExtension(subFields, a * b, "Q(αβ)");
        FindExtension(subFields, e, "K");
    }
    
    /***
        With P = X^8 - 24*X^6 + 144*X^4 - 288*X^2 + 144
        Factorization in Q(y)/Q
            X - y
            X + 1/12*y^5 - 3/2*y^3 + 3*y
            X + y
            X - 1/12*y^5 + 3/2*y^3 - 3*y
            X + 1/12*y^7 - 11/6*y^5 + 17/2*y^3 - 10*y
            X + 1/24*y^7 - 5/6*y^5 + 5/2*y^3 + y
            X - 1/12*y^7 + 11/6*y^5 - 17/2*y^3 + 10*y
            X - 1/24*y^7 + 5/6*y^5 - 5/2*y^3 - y

        Galois Group
        |G1| = 8
        Type        NonAbelianGroup
        BaseGroup   S8
        SuperGroup  |Gal( Q(y)/Q )| = 8

        Elements
        (1)[1] = []
        (2)[2] = [(1 2)(3 4)(5 8)(6 7)]
        (3)[4] = [(1 3 2 4)(5 7 8 6)]
        (4)[4] = [(1 4 2 3)(5 6 8 7)]
        (5)[4] = [(1 5 2 8)(3 6 4 7)]
        (6)[4] = [(1 6 2 7)(3 8 4 5)]
        (7)[4] = [(1 7 2 6)(3 5 4 8)]
        (8)[4] = [(1 8 2 5)(3 7 4 6)]

        Tower 1
          |G6| = 1   => [Q(y):Q] = 8
          |G5| = 2   => [Q(d):Q] = 4 with d = -1/2*y^2
          |G2| = 4   => [Q(a):Q] = 2 with a = 1/24*y^6 - 5/6*y^4 + 5/2*y^2
          |G1| = 8   => [Q:Q]    = 1
        Tower 2
          |G6| = 1   => [Q(y):Q] = 8
          |G5| = 2   => [Q(d):Q] = 4 with d = -1/2*y^2
          |G3| = 4   => [Q(b):Q] = 2 with b = 1/24*y^6 - y^4 + 11/2*y^2
          |G1| = 8   => [Q:Q]    = 1
        Tower 3
          |G6| = 1   => [Q(y):Q] = 8
          |G5| = 2   => [Q(d):Q] = 4 with d = -1/2*y^2
          |G4| = 4   => [Q(c):Q] = 2 with c = 1/12*y^6 - 7/4*y^4 + 7*y^2
          |G1| = 8   => [Q:Q]    = 1

        Q = Q(z) is subfield of Q and [Q:Q] = 1
        With z - 1 = 0
        With y^8 - 24*y^6 + 144*y^4 - 288*y^2 + 144 = 0
        [|G6| = 1       [Q(y):Q] = 8    z = 1]
        [|G5| = 2       [Q(d):Q] = 4    z = 1]
        [|G2| = 4       [Q(a):Q] = 2    z = 1]
        [|G3| = 4       [Q(b):Q] = 2    z = 1]
        [|G4| = 4       [Q(c):Q] = 2    z = 1]
        [|G1| = 8       [Q:Q] = 1       z = 1]

        Q(a+b) = Q(z) is subfield of Q(d) and [Q(a+b):Q] = 4
        With z^4 - 10*z^2 + 1 = 0
        With y^8 - 24*y^6 + 144*y^4 - 288*y^2 + 144 = 0
        [|G6| = 1       [Q(y):Q] = 8    z = -1/8*y^6 + 31/12*y^4 - 19/2*y^2 + 6]
        [|G5| = 2       [Q(d):Q] = 4    z = d^3 + 31/3*d^2 + 19*d + 6          ]

        Q(a) = Q(z) is subfield of Q(a) and [Q(a):Q] = 2
        With z^2 - 2 = 0
        With y^8 - 24*y^6 + 144*y^4 - 288*y^2 + 144 = 0
        [|G6| = 1       [Q(y):Q] = 8    z = -1/24*y^6 + 5/6*y^4 - 5/2*y^2]
        [|G5| = 2       [Q(d):Q] = 4    z = 1/3*d^3 + 10/3*d^2 + 5*d     ]
        [|G2| = 4       [Q(a):Q] = 2    z = -a                           ]

        Q(b) = Q(z) is subfield of Q(c) and [Q(b):Q] = 2
        With z^2 - 3 = 0
        With y^8 - 24*y^6 + 144*y^4 - 288*y^2 + 144 = 0
        [|G6| = 1       [Q(y):Q] = 8    z = -1/12*y^6 + 7/4*y^4 - 7*y^2 + 6]
        [|G5| = 2       [Q(d):Q] = 4    z = 2/3*d^3 + 7*d^2 + 14*d + 6     ]
        [|G4| = 4       [Q(c):Q] = 2    z = -c + 6                         ]

        Q(ab) = Q(z) is subfield of Q(b) and [Q(ab):Q] = 2
        With z^2 - 6 = 0
        With y^8 - 24*y^6 + 144*y^4 - 288*y^2 + 144 = 0
        [|G6| = 1       [Q(y):Q] = 8    z = 1/24*y^6 - y^4 + 11/2*y^2 - 6]
        [|G5| = 2       [Q(d):Q] = 4    z = -1/3*d^3 - 4*d^2 - 11*d - 6  ]
        [|G3| = 4       [Q(b):Q] = 2    z = b - 6                        ]

        K = Q(z) is subfield of Q(y) and [K:Q] = 8
        With z^8 - 24*z^6 + 144*z^4 - 288*z^2 + 144 = 0
        With y^8 - 24*y^6 + 144*y^4 - 288*y^2 + 144 = 0
        [|G6| = 1       [Q(y):Q] = 8    z = y]
     */
}
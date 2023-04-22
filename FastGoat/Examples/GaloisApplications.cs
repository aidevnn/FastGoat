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
        var sel = galCor.Select(e => (e, pol: GaloisTheory.Rewrite(e.primElt, a)))
            .Where(e => !e.pol.IsZero()).OrderBy(e => e.e.SubGr.Count()).ToArray();

        if (sel.Length == 0)
        {
            Console.WriteLine("Not Found");
            Console.WriteLine();
        }
        else
        {
            var minPol = IntFactorisation.CharacPoly2(a);
            var rows = sel.Select(e => new[] { e.e.SubGr.ShortName, $"w = {e.e.primElt}", $"z = {e.pol.SubstituteChar('w')}" })
                .ToArray();

            var dim = GaloisTheory.GetBase(sel.Last().e.primElt).Length;
            Console.WriteLine($"{name} = Q(z) is subfield of Q(w) and [{name}:Q] = {dim}");
            Console.WriteLine($"With {minPol.SubstituteChar('z')} = 0");
            Console.WriteLine($"With {a.F} = 0");
            Ring.DisplayMatrix(Ring.Matrix(rows.Length, rows.SelectMany(e => e).ToArray()), Ring.MatrixDisplay.TableLeft, "  \t");
            Console.WriteLine();
        }
    }

    public static void GaloisCorrespondence(KPoly<Rational> P, int nbGens = 2)
    {
        var subFields = GaloisTheory.SubFields(P).ToArray();
        var extTowers = ExtensionsTower(subFields);
        GaloisCorrespondence(extTowers);
    }
    
    public static void GaloisCorrespondence(GaloisCorrespondence[][] exts)
    {
        var subFields = exts.SelectMany(e => e.Select(gc => gc.primElt)).Distinct().ToArray();

        var galCor = exts[0].Last();
        if (galCor is null)
            throw new();

        var P = galCor.primElt.F;
        var field2name = subFields.Where(e => !e.Equals(e.X) && e.Degree > 0).Select((e, k) => (e, $"{(char)('a' + k)}"))
            .ToDictionary(e => e.e, e => e.Item2);
        field2name[galCor.primElt] = "1";
        field2name[galCor.primElt.X] = "y";
        
        var (X, _) = FG.EPolyXc(galCor.primElt.F, 'y');
        Console.WriteLine($"With P = {P.Substitute(X)}");
        galCor.roots.Select(r => X - r).Println("Factorization in Q(y)/Q");
        Console.WriteLine();
        Console.WriteLine("Galois Group");
        DisplayGroup.HeadElements(galCor.SubGr);
        var i = 1;
        foreach (var lt in exts.OrderBy(t => t.Select(gc => gc.SubGr.Count()).ToArray(),
                     Comparer<int[]>.Create((l1, l2) => l1.SequenceCompareTo(l2)))
                     .ThenBy(t => t.Select(gc => field2name[gc.primElt]).ToArray(),
                         Comparer<string[]>.Create((l1, l2) => l1.SequenceCompareTo(l2))))
        {
            Console.WriteLine($"Tower {i++}");
            foreach (var gc in lt)
            {
                var a = field2name[gc.primElt];
                if (a == "1")
                    Console.WriteLine($"  {gc.SubGr.ShortName,-10} => [Q:Q]    = 1");
                else if (a == "y")
                    Console.WriteLine($"  {gc.SubGr.ShortName,-10} => [Q(y):Q] = {P.Degree}");
                else
                    Console.WriteLine($"  {gc.SubGr.ShortName,-10} => [Q({a}):Q] = {GaloisTheory.GetBase(gc.primElt).Length} with {a}={gc.primElt}");
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
        FindExtension(subFields, a + b, "Q(a+b)");
        FindExtension(subFields, a, "Q(a)");
        FindExtension(subFields, b, "Q(b)");
        FindExtension(subFields, a * b, "Q(ab)");
        FindExtension(subFields, e, "K");
    }
    
    /***
        With P = X^8 + -24*X^6 + 144*X^4 + -288*X^2 + 144
        Factorization in Q(y)/Q
            X + -y
            X + 1/12*y^5 + -3/2*y^3 + 3*y
            X + y
            X + -1/12*y^5 + 3/2*y^3 + -3*y
            X + 1/12*y^7 + -11/6*y^5 + 17/2*y^3 + -10*y
            X + 1/24*y^7 + -5/6*y^5 + 5/2*y^3 + y
            X + -1/12*y^7 + 11/6*y^5 + -17/2*y^3 + 10*y
            X + -1/24*y^7 + 5/6*y^5 + -5/2*y^3 + -y

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
          |G5| = 2   => [Q(a):Q] = 4 with a=-y^2
          |G2| = 4   => [Q(b):Q] = 2 with b=y^6 + -20*y^4 + 60*y^2
          |G1| = 8   => [Q:Q]    = 1
        Tower 2
          |G6| = 1   => [Q(y):Q] = 8
          |G5| = 2   => [Q(a):Q] = 4 with a=-y^2
          |G3| = 4   => [Q(c):Q] = 2 with c=y^6 + -24*y^4 + 132*y^2
          |G1| = 8   => [Q:Q]    = 1
        Tower 3
          |G6| = 1   => [Q(y):Q] = 8
          |G5| = 2   => [Q(a):Q] = 4 with a=-y^2
          |G4| = 4   => [Q(d):Q] = 2 with d=y^6 + -21*y^4 + 84*y^2
          |G1| = 8   => [Q:Q]    = 1

        Q = Q(z) is subfield of Q(w) and [Q:Q] = 1
        With z + -1 = 0
        With y^8 + -24*y^6 + 144*y^4 + -288*y^2 + 144 = 0
        [|G6| = 1       w = y                           z = 1]
        [|G5| = 2       w = -y^2                        z = 1]
        [|G2| = 4       w = y^6 + -20*y^4 + 60*y^2      z = 1]
        [|G3| = 4       w = y^6 + -24*y^4 + 132*y^2     z = 1]
        [|G4| = 4       w = y^6 + -21*y^4 + 84*y^2      z = 1]
        [|G1| = 8       w = 1                           z = 1]

        Q(a+b) = Q(z) is subfield of Q(w) and [Q(a+b):Q] = 4
        With z^4 + -10*z^2 + 1 = 0
        With y^8 + -24*y^6 + 144*y^4 + -288*y^2 + 144 = 0
        [|G6| = 1       w = y           z = -1/8*w^6 + 31/12*w^4 + -19/2*w^2 + 6]
        [|G5| = 2       w = -y^2        z = 1/8*w^3 + 31/12*w^2 + 19/2*w + 6    ]

        Q(a) = Q(z) is subfield of Q(w) and [Q(a):Q] = 2
        With z^2 + -2 = 0
        With y^8 + -24*y^6 + 144*y^4 + -288*y^2 + 144 = 0
        [|G6| = 1       w = y                           z = -1/24*w^6 + 5/6*w^4 + -5/2*w^2]
        [|G5| = 2       w = -y^2                        z = 1/24*w^3 + 5/6*w^2 + 5/2*w    ]
        [|G2| = 4       w = y^6 + -20*y^4 + 60*y^2      z = -1/24*w                       ]

        Q(b) = Q(z) is subfield of Q(w) and [Q(b):Q] = 2
        With z^2 + -3 = 0
        With y^8 + -24*y^6 + 144*y^4 + -288*y^2 + 144 = 0
        [|G6| = 1       w = y                           z = -1/12*w^6 + 7/4*w^4 + -7*w^2 + 6]
        [|G5| = 2       w = -y^2                        z = 1/12*w^3 + 7/4*w^2 + 7*w + 6    ]
        [|G4| = 4       w = y^6 + -21*y^4 + 84*y^2      z = -1/12*w + 6                     ]

        Q(ab) = Q(z) is subfield of Q(w) and [Q(ab):Q] = 2
        With z^2 + -6 = 0
        With y^8 + -24*y^6 + 144*y^4 + -288*y^2 + 144 = 0
        [|G6| = 1       w = y                           z = 1/24*w^6 + -w^4 + 11/2*w^2 + -6]
        [|G5| = 2       w = -y^2                        z = -1/24*w^3 + -w^2 + -11/2*w + -6]
        [|G3| = 4       w = y^6 + -24*y^4 + 132*y^2     z = 1/24*w + -6                    ]

        K = Q(z) is subfield of Q(w) and [K:Q] = 8
        With z^8 + -24*z^6 + 144*z^4 + -288*z^2 + 144 = 0
        With y^8 + -24*y^6 + 144*y^4 + -288*y^2 + 144 = 0
        [|G6| = 1       w = y   z = w]
     */
}
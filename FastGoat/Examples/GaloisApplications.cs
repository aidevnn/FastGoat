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
        With P = X^8 + -44·X^6 + -96·X^5 + 222·X^4 + 672·X^3 + -92·X^2 + -1152·X + -647
        Factorization in Q(y)/Q
            X + -y
            X + 15/92·y^7 + -39/92·y^6 + -285/46·y^5 + 15/23·y^4 + 3605/92·y^3 + 1139/92·y^2 + -2923/46·y + -916/23
            X + -5/46·y^7 + 27/92·y^6 + 94/23·y^5 + -77/92·y^4 + -51/2·y^3 + -403/92·y^2 + 956/23·y + 1645/92
            X + -5/92·y^7 + 3/23·y^6 + 97/46·y^5 + 17/92·y^4 + -1259/92·y^3 + -8·y^2 + 1057/46·y + 2019/92
            X + 21/184·y^7 + -37/184·y^6 + -827/184·y^5 + -581/184·y^4 + 4443/184·y^3 + 4045/184·y^2 + -6045/184·y + -5827/184
            X + 1/46·y^7 + -1/23·y^6 + -45/46·y^5 + -9/46·y^4 + 224/23·y^3 + 371/46·y^2 + -568/23·y + -537/23
            X + -9/184·y^7 + 7/184·y^6 + 403/184·y^5 + 511/184·y^4 + -2815/184·y^3 + -3663/184·y^2 + 5549/184·y + 5913/184
            X + -2/23·y^7 + 19/92·y^6 + 151/46·y^5 + 53/92·y^4 + -855/46·y^3 + -933/92·y^2 + 630/23·y + 2105/92

        Galois Group
        |G1| = 8
        Type        NonAbelianGroup
        BaseGroup   S8
        SuperGroup  |Gal( Q(y)/Q )| = 8

        Elements
        (1)[1] = []
        (2)[2] = [(1 8)(2 6)(3 5)(4 7)]
        (3)[4] = [(1 2 8 6)(3 4 5 7)]
        (4)[4] = [(1 3 8 5)(2 7 6 4)]
        (5)[4] = [(1 4 8 7)(2 3 6 5)]
        (6)[4] = [(1 5 8 3)(2 4 6 7)]
        (7)[4] = [(1 6 8 2)(3 7 5 4)]
        (8)[4] = [(1 7 8 4)(2 5 6 3)]

        Tower 1
          |G6| = 1   => [Q(y):Q] = 8
          |G5| = 2   => [Q(a):Q] = 4 with a=y^5 + -35/16·y^4 + -143/4·y^3 + -143/8·y^2 + 475/4·y
          |G2| = 4   => [Q(b):Q] = 2 with b=y^7 + -53/26·y^6 + -514/13·y^5 + -427/26·y^4 + 3181/13·y^3 + 4705/26·y^2 + -396·y
          |G1| = 8   => [Q:Q]    = 1
        Tower 2
          |G6| = 1   => [Q(y):Q] = 8
          |G5| = 2   => [Q(a):Q] = 4 with a=y^5 + -35/16·y^4 + -143/4·y^3 + -143/8·y^2 + 475/4·y
          |G3| = 4   => [Q(c):Q] = 2 with c=y^7 + -3·y^6 + -41·y^5 + 28·y^4 + 383·y^3 + 53·y^2 + -871·y
          |G1| = 8   => [Q:Q]    = 1
        Tower 3
          |G6| = 1   => [Q(y):Q] = 8
          |G5| = 2   => [Q(a):Q] = 4 with a=y^5 + -35/16·y^4 + -143/4·y^3 + -143/8·y^2 + 475/4·y
          |G4| = 4   => [Q(d):Q] = 2 with d=y^7 + -21/8·y^6 + -147/4·y^5 + 21/8·y^4 + 395/2·y^3 + 297/8·y^2 + -995/4·y
          |G1| = 8   => [Q:Q]    = 1

        Q = Q(z) is subfield of Q(w) and [Q:Q] = 1
        With z + -1 = 0
        With y^8 + -44·y^6 + -96·y^5 + 222·y^4 + 672·y^3 + -92·y^2 + -1152·y + -647 = 0
        [|G6| = 1       w = y                                                                                   z = 1]
        [|G5| = 2       w = y^5 + -35/16·y^4 + -143/4·y^3 + -143/8·y^2 + 475/4·y                                z = 1]
        [|G2| = 4       w = y^7 + -53/26·y^6 + -514/13·y^5 + -427/26·y^4 + 3181/13·y^3 + 4705/26·y^2 + -396·y   z = 1]
        [|G3| = 4       w = y^7 + -3·y^6 + -41·y^5 + 28·y^4 + 383·y^3 + 53·y^2 + -871·y                         z = 1]
        [|G4| = 4       w = y^7 + -21/8·y^6 + -147/4·y^5 + 21/8·y^4 + 395/2·y^3 + 297/8·y^2 + -995/4·y          z = 1]
        [|G1| = 8       w = -1                                                                                  z = 1]

        Q(a+b) = Q(z) is subfield of Q(w) and [Q(a+b):Q] = 4
        With z^4 + -10·z^2 + 1 = 0
        With y^8 + -44·y^6 + -96·y^5 + 222·y^4 + 672·y^3 + -92·y^2 + -1152·y + -647 = 0
        [|G6| = 1       w = y                                                           z = -5/92·w^7 + 27/184·w^6 + 47/23·w^5 + -77/184·w^4 + -51/4·w^3 + -403/184·w^2 + 933/46·w + 1645/184       ]
        [|G5| = 2       w = y^5 + -35/16·y^4 + -143/4·y^3 + -143/8·y^2 + 475/4·y        z = -89776/477201907·w^3 + -21797565/477201907·w^2 + -26371453865/7635230512·w + -9807902104215/122163688192]

        Q(a) = Q(z) is subfield of Q(w) and [Q(a):Q] = 2
        With z^2 + -2 = 0
        With y^8 + -44·y^6 + -96·y^5 + 222·y^4 + 672·y^3 + -92·y^2 + -1152·y + -647 = 0
        [|G6| = 1       w = y                                                                                   z = -1/23·w^7 + 21/184·w^6 + 147/92·w^5 + -21/184·w^4 + -395/46·w^3 + -297/184·w^2 + 995/92·w + 801/184]
        [|G5| = 2       w = y^5 + -35/16·y^4 + -143/4·y^3 + -143/8·y^2 + 475/4·y                                z = -400/2963987·w^3 + -97671/2963987·w^2 + -121143891/47423792·w + -47336300197/758780672             ]
        [|G4| = 4       w = y^7 + -21/8·y^6 + -147/4·y^5 + 21/8·y^4 + 395/2·y^3 + 297/8·y^2 + -995/4·y          z = -1/23·w + 801/184                                                                                  ]

        Q(b) = Q(z) is subfield of Q(w) and [Q(b):Q] = 2
        With z^2 + -3 = 0
        With y^8 + -44·y^6 + -96·y^5 + 222·y^4 + 672·y^3 + -92·y^2 + -1152·y + -647 = 0
        [|G6| = 1       w = y                                                                   z = -1/92·w^7 + 3/92·w^6 + 41/92·w^5 + -7/23·w^4 + -383/92·w^3 + -53/92·w^2 + 871/92·w + 211/46 ]
        [|G5| = 2       w = y^5 + -35/16·y^4 + -143/4·y^3 + -143/8·y^2 + 475/4·y                z = -1952/36707839·w^3 + -467118/36707839·w^2 + -264126439/293662712·w + -84106068173/4698603392]
        [|G3| = 4       w = y^7 + -3·y^6 + -41·y^5 + 28·y^4 + 383·y^3 + 53·y^2 + -871·y         z = -1/92·w + 211/46                                                                            ]

        Q(ab) = Q(z) is subfield of Q(w) and [Q(ab):Q] = 2
        With z^2 + -6 = 0
        With y^8 + -44·y^6 + -96·y^5 + 222·y^4 + 672·y^3 + -92·y^2 + -1152·y + -647 = 0
        [|G6| = 1       w = y                                                                                   z = -13/92·w^7 + 53/184·w^6 + 257/46·w^5 + 427/184·w^4 + -3181/92·w^3 + -4705/184·w^2 + 1287/23·w + 9033/184]
        [|G5| = 2       w = y^5 + -35/16·y^4 + -143/4·y^3 + -143/8·y^2 + 475/4·y                                z = -48/1107197·w^3 + -17557/1107197·w^2 + -29511761/17715152·w + -14120487551/283442432                    ]
        [|G2| = 4       w = y^7 + -53/26·y^6 + -514/13·y^5 + -427/26·y^4 + 3181/13·y^3 + 4705/26·y^2 + -396·y   z = -13/92·w + 9033/184                                                                                     ]

        K = Q(z) is subfield of Q(w) and [K:Q] = 8
        With z^8 + -24·z^6 + 144·z^4 + -288·z^2 + 144 = 0
        With y^8 + -44·y^6 + -96·y^5 + 222·y^4 + 672·y^3 + -92·y^2 + -1152·y + -647 = 0
        [|G6| = 1       w = y   z = -5/92·w^7 + 27/184·w^6 + 47/23·w^5 + -77/184·w^4 + -51/4·w^3 + -403/184·w^2 + 979/46·w + 1645/184]
     */
}
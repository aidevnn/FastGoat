using FastGoat.Commons;
using FastGoat.Examples;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Characters;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Matrix;
using FastGoat.UserGroup.Perms;
using FastGoat.UserGroup.Polynoms;
using FastGoat.UserGroup.Words;
using FastGoat.UserGroup.Words.ToddCoxeter;
using static FastGoat.Commons.IntExt;
using static FastGoat.Commons.EnumerableExt;
using FastGoat.UserGroup.Padic;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");

RngSeed(1234);

KPoly<K> M1<K>(KPoly<K> g, int q) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
{
    var h = g.Zero;
    var g0 = g.One;
    if (Int32.IsPow2(q))
    {
        for (int i = 0; i < q; i++)
        {
            if (i == 0 || Int32.IsPow2(i))
                h += g0;

            g0 *= g;
        }
    }
    else
    {
        // for (int i = 0; i < (q - 1) / 2; i++)
        //     g0 *= g;
        // h = g0 - 1;
        h = g.Pow((q - 1) / 2) - 1;
    }

    return h;
}

// A Computational Introduction to Number Theory and Algebra
// Victor Shoup
// 20.5 Factoring polynomials: Berlekamp’s algorithm page 541
IEnumerable<KPoly<K>> Berlekamp2<K>(KPoly<K> f, K a0) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
{
    var acc = a0.One;
    var allF = new List<K>() { a0.Zero };
    var q = 0;
    do
    {
        allF.Add(acc);
        acc *= a0;
        ++q;
    } while (!acc.Equals(a0.One));

    var polys = PolynomialFactorization.FrobeniusKernel(f, q + 1);
    var r = polys.Length;

    var H0 = new List<KPoly<K>>() { f };
    var H1 = new List<KPoly<K>>();
    while (H0.Count < r)
    {
        var g = polys.Aggregate(f.Zero, (sum, gi) => sum + allF[Rng.Next(q + 1)] * gi);
        H1.Clear();
        foreach (var h in H0)
        {
            var b = g.Div(h).rem;
            var d = Ring.Gcd(M1(b, q + 1), h).Monic;
            if (d.Equals(d.One) || d.Equals(h))
                H1.Add(h.Monic);
            else
                H1.AddRange(new[] { d, (h / d).Monic });
        }

        H0.Clear();
        H0.AddRange(H1);
    }

    return H0.Order();
}

// AECF Algorithme de Berlekamp 353
IEnumerable<KPoly<K>> Berlekamp3<K>(KPoly<K> f, K a0) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
{
    // var irrs = PolynomialFactorization.Firr(f, a0).Order().ToArray();
    // var r = irrs.Length.Range();
    // foreach (var (i, j) in r.Grid2D(r))
    // {
    //     var w = FG.EPoly(irrs[j], 'w');
    //     var df = f.Derivative.Substitute(w);
    //     var b0 = f / irrs[i] * irrs[i].Derivative;
    //     var bi = b0.Substitute(w) / df;
    //     Console.WriteLine(((i, j), bi)); // (Bi mod Fj) = (i==j ? 1 : 0)
    // }

    var acc = a0.One;
    var allF = new List<K>() { a0.Zero };
    do
    {
        allF.Add(acc);
        acc *= a0;
    } while (!acc.Equals(a0.One));

    var Gi = PolynomialFactorization.FrobeniusKernel(f, allF.Count);
    return BerlekampRec(f, Gi, allF);
}

// AECF Algorithme de Berlekamp 351
IEnumerable<KPoly<K>> BerlekampRec<K>(KPoly<K> F, KPoly<K>[] Gi, List<K> allF) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
{
    if (Gi.Length == 1)
    {
        yield return F.Monic;
        yield break;
    }

    KPoly<K> H1;
    var q = allF.Count;
    var degF = F.Degree;
    while (true)
    {
        var G = Gi.Aggregate(F.Zero, (sum, gi) => sum + allF[Rng.Next(q)] * gi);
        H1 = Ring.Gcd(F, G);
        if (H1.Degree > 0 && H1.Degree != degF)
            break;

        var H = M1(G, q);
        H1 = Ring.Gcd(F, H);
        if (H1.Degree > 0 && H1.Degree != degF)
            break;
    }

    var H2 = F / H1;
    
    var Gi1 = new HashSet<KPoly<K>>() { F.One };
    var Gi2 = new HashSet<KPoly<K>>() { F.One };
    
    var df = F.Derivative;
    var h1 = FG.EPoly(H1);
    var h2 = FG.EPoly(H2);
    var dh1h2 = (H1.Derivative * H2).Substitute(h1).Inv();
    var h1dh2 = (H1 * H2.Derivative).Substitute(h2).Inv();
    foreach (var gi in Gi)
    {
        var gidf = gi * df;
        
        var gi1 = (gidf * dh1h2).Poly;
        if (!gi1.IsZero())
            Gi1.Add(gi1.Monic);

        var gi2 = (gidf * h1dh2).Poly;
        if (!gi2.IsZero())
            Gi2.Add(gi2.Monic);
    }

    foreach (var f in BerlekampRec(H1, Gi1.ToArray(), allF))
        yield return f;

    foreach (var f in BerlekampRec(H2, Gi2.ToArray(), allF))
        yield return f;
}

{
    Monom.Display = MonomDisplay.StarCaret;
    for (int i = 0; i < 10; ++i)
    {
        var p = Primes10000[Rng.Next(5)]; // 2, 3, 5, 7, 11
        var d = Rng.Next((int)(Math.Log(50) / Math.Log(p))) + 1; // p^d < 50 => 4, 8, 9, 16, 25, 27, 32, 49
        var fq = new Fq(p.Pow(d), 'a');
        var gf = FG.Galois(p.Pow(d), 'a');
        var a0 = gf.GetGenerators().First();
        var n = 2 + Rng.Next(6);
        var f0 = PolynomialFactorization.RandPolySep(fq.One, p, n);
        var f1 = PolynomialFactorization.RandPolySep(fq.One, p, n);
        var f2 = PolynomialFactorization.RandPolySep(fq.One, p, n);
        var f = f0 * f1 * f2;
        if (Ring.Discriminant(f).IsZero())
        {
            --i;
            continue;
        }

        Console.WriteLine($"{fq} with {fq.F} = 0");
        Console.WriteLine($"f = {f} mod ({p})");
        var firr0 = PolynomialFactorization.Firr(f, a0).Order().ToArray();
        Console.WriteLine($"Disc(f) = {Ring.Discriminant(f)} mod ({p})");
        Console.WriteLine($"Fact1(f) = {firr0.Glue("*", "({0})")} mod ({p})");

        var firr1 = Berlekamp2(f, a0).Order().ToArray();
        Console.WriteLine($"Fact2(f) = {firr1.Glue("*", "({0})")} mod ({p})");

        var firr2 = Berlekamp3(f, a0).Order().ToArray();
        Console.WriteLine($"Fact3(f) = {firr2.Glue("*", "({0})")} mod ({p})");
        var check2 = firr0.SequenceEqual(firr1);
        var check3 = firr0.SequenceEqual(firr2);
        Console.WriteLine($"Check2 : {check2}");
        Console.WriteLine($"Check3 : {check3}");

        if (!check2 || !check3)
            throw new();

        var nb = 5;
        GlobalStopWatch.Bench(nb, "B1", () => PolynomialFactorization.Firr(f, a0).Order().ToArray());
        GlobalStopWatch.Bench(nb, "B2", () => Berlekamp2(f, a0).ToArray());
        GlobalStopWatch.Bench(nb, "B3", () => Berlekamp3(f, a0).ToArray());

        Console.WriteLine();
    }
}

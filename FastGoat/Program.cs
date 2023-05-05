using System.Diagnostics.Tracing;
using System.Globalization;
using System.Numerics;
using System.Security.Principal;
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


BigReal GenBR((int M0, int M1) m, (int E0, int E1) e, int O)
{
    var a0 = Rng.Next(m.M0, m.M1 + 1);
    var e0 = Rng.Next(e.E0, e.E1 + 1);
    return BigReal.FromBigInteger(a0, O).Mul10PowN(e0);
}

string fmt(double a, int O) => String.Format($"{{0:E{O}}}", a);

void BigRealAdd()
{
    RngSeed(1231545);
    var k = 20000;
    var M = 1000;
    var E = 6;
    var O = 14;
    var lt = k.Range()
        .Select(_ => (GenBR((-M * M, M * M), (-2 * E, 2 * E), O), GenBR((-M, M), (-E, E), O)))
        .Select(e => Rng.NextDouble() < 0.025 ? (e.Item1, e.Item1) : Rng.NextDouble() < 0.5 ? e : (e.Item2, e.Item1))
        .ToArray();
    GlobalStopWatch.AddLap();
    var err = Double.Pow(10, -O + 2);
    foreach (var (a, b) in lt)
    {
        var a_b = a + b;
        var a_b2 = a.ToDouble + b.ToDouble;
        var diff = Double.Abs(a_b.ToDouble - a_b2) * Double.Pow(10, -a_b.V);
        if ((!a_b.IsZero() && diff > err) || (a_b.IsZero() && Double.Abs(a_b2) > err))
        {
            Console.WriteLine(new
                { a, ad = $"{fmt(a.ToDouble, O)}", a0 = a.Details, b, bd = $"{fmt(b.ToDouble, O)}", b0 = b.Details });
            Console.WriteLine(
                new
                {
                    a_b, a_bd = $"{fmt(a_b.ToDouble, O)}", a_bd2 = $"{fmt(a_b2, O)}", diff,
                    epsilon = err
                });
            Console.WriteLine();
        }
    }

    GlobalStopWatch.Show("END BigRealAdd");
}

void BigRealMul()
{
    // RngSeed(1231545);
    var k = 20000;
    var M = 1000;
    var E = 6;
    var O = 14;
    var lt = k.Range()
        .Select(_ => (GenBR((-M * M, M * M), (-2 * E, 2 * E), O), GenBR((-M, M), (-E, E), O)))
        .Select(e => Rng.NextDouble() < 0.025 ? (e.Item1, e.Item1) : Rng.NextDouble() < 0.5 ? e : (e.Item2, e.Item1))
        .ToArray();
    GlobalStopWatch.AddLap();
    var err = Double.Pow(10, -O + 2);
    foreach (var (a, b) in lt)
    {
        var a_b = a * b;
        var a_b2 = a.ToDouble * b.ToDouble;
        var diff = Double.Abs(a_b.ToDouble - a_b2) * Double.Pow(10, -a_b.V);
        if ((!a_b.IsZero() && diff > err) || (a_b.IsZero() && Double.Abs(a_b2) > err))
        {
            Console.WriteLine(new
                { a, ad = $"{fmt(a.ToDouble, O)}", a0 = a.Details, b, bd = $"{fmt(b.ToDouble, O)}", b0 = b.Details });
            Console.WriteLine(
                new
                {
                    a_b, a_bd = $"{fmt(a_b.ToDouble, O)}", a_bd2 = $"{fmt(a_b2, O)}", diff,
                    epsilon = err
                });
            Console.WriteLine();
        }
    }

    GlobalStopWatch.Show("END BigRealMul");
}

void BigRealKMul()
{
    // RngSeed(1231545);
    var k = 20000;
    var M = 1000;
    var E = 6;
    var O = 14;
    var lt = k.Range().Select(_ => (GenBR((-M * M, M * M), (-2 * E, 2 * E), O), Rng.Next(-M, M))).ToArray();
    GlobalStopWatch.AddLap();
    var err = Double.Pow(10, -O + 2);
    foreach (var (a, b) in lt)
    {
        var a_b = a * b;
        var a_b2 = a.ToDouble * b;
        var diff = Double.Abs(a_b.ToDouble - a_b2) * Double.Pow(10, -a_b.V);
        if ((!a_b.IsZero() && diff > err) || (a_b.IsZero() && Double.Abs(a_b2) > err))
        {
            Console.WriteLine(new
                { a, ad = $"{fmt(a.ToDouble, O)}", a0 = a.Details, b, b0 = BigReal.FromBigInteger(b, O) });
            Console.WriteLine(
                new
                {
                    a_b, a_bd = $"{fmt(a_b.ToDouble, O)}", a_bd2 = $"{fmt(a_b2, O)}", diff,
                    epsilon = err
                });
            Console.WriteLine();
        }
    }

    GlobalStopWatch.Show("END BigRealKMul");
}

void BigRealDiv()
{
    // RngSeed(1231545);
    var k = 20000;
    var M = 1000;
    var E = 6;
    var O = 14;
    var lt = k.Range()
        .Select(_ => (GenBR((-M * M, M * M), (-2 * E, 2 * E), O), GenBR((-M, M), (-E, E), O)))
        .Select(e => Rng.NextDouble() < 0.025 ? (e.Item1, e.Item1) : Rng.NextDouble() < 0.5 ? e : (e.Item2, e.Item1))
        .ToArray();
    GlobalStopWatch.AddLap();
    var err = Double.Pow(10, -O + 2);
    foreach (var (a, b) in lt)
    {
        if (b.IsZero())
            continue;

        var a_b = a / b;
        var a_b2 = a.ToDouble / b.ToDouble;
        var diff = Double.Abs(a_b.ToDouble - a_b2) * Double.Pow(10, -a_b.V);
        if ((!a_b.IsZero() && diff > err) || (a_b.IsZero() && Double.Abs(a_b2) > err))
        {
            Console.WriteLine(new
                { a, ad = $"{fmt(a.ToDouble, O)}", a0 = a.Details, b });
            Console.WriteLine(
                new
                {
                    a_b, a_bd = $"{fmt(a_b.ToDouble, O)}", a_bd2 = $"{fmt(a_b2, O)}", diff,
                    epsilon = err
                });
            Console.WriteLine();
        }
    }

    GlobalStopWatch.Show("END BigRealDiv");
}

void BigRealStr()
{
    // RngSeed(1231545);
    var k = 20000;
    var M = 1000;
    var E = 6;
    var O = 14;
    var lt = k.Range().Select(_ => GenBR((-M * M, M * M), (-2 * E, 2 * E), O)).ToArray();
    GlobalStopWatch.AddLap();
    var err = Double.Pow(10, -O + 2);
    foreach (var a in lt)
    {
        var a_b = a;
        var a_b2 = Double.Parse($"{a}");
        var diff = Double.Abs(a_b.ToDouble - a_b2) * Double.Pow(10, -a_b.V);
        if (diff > err)
        {
            Console.WriteLine(new
                { a, ad = $"{fmt(a.ToDouble, O)}", a0 = a.Details });
            Console.WriteLine(
                new
                {
                    a_b, a_bd = $"{fmt(a_b.ToDouble, O)}", a_bd2 = $"{fmt(a_b2, O)}", diff,
                    epsilon = err
                });
            Console.WriteLine();
        }
    }

    GlobalStopWatch.Show("END BigRealStr");
}

void BigRealBInt()
{
    // RngSeed(1231545);
    var k = 20000;
    var M = 1000;
    var E = 6;
    var O = 14;
    var lt = k.Range().Select(_ => Rng.Next(-M * M, M * M) * BigInteger.Pow(10, Rng.Next(0, E))).ToArray();
    GlobalStopWatch.AddLap();
    var err = Double.Pow(10, -O + 2);
    foreach (var a in lt)
    {
        var a_b = BigReal.FromBigInteger(a, O);
        var a_b2 = Double.Parse($"{a}");
        var diff = Double.Abs(a_b.ToDouble - a_b2) * Double.Pow(10, -a_b.V);
        if (diff > err)
        {
            Console.WriteLine(
                new
                {
                    a, a_b.Details, a_b, a_bd = $"{fmt(a_b.ToDouble, O)}", a_bd2 = $"{fmt(a_b2, O)}", diff,
                    epsilon = err
                });
            Console.WriteLine();
        }
    }

    GlobalStopWatch.Show("END BigRealBInt");
}

void BigRealRat()
{
    // RngSeed(1231545);
    var k = 20000;
    var M = 1000;
    var E = 6;
    var O = 14;
    var lt = k.Range().Select(_ => GenBR((-M * M, M * M), (-2 * E, 2 * E), O)).ToArray();
    GlobalStopWatch.AddLap();
    var err = Double.Pow(10, -O + 2);
    foreach (var a in lt)
    {
        var a_b = a;
        var a_b2 = (double)a.ToRational;
        var diff = Double.Abs(a_b.ToDouble - a_b2) * Double.Pow(10, -a_b.V);
        if ((!a_b.IsZero() && diff > err) || (a_b.IsZero() && Double.Abs(a_b2) > err))
        {
            Console.WriteLine(
                new
                {
                    a, a_b.Details, a_b, a_bd = $"{fmt(a_b.ToDouble, O)}", a_bd2 = $"{fmt(a_b2, O)}", diff,
                    epsilon = err
                });
            Console.WriteLine();
        }
    }

    GlobalStopWatch.Show("END BigRealRat");
}

void BigRealInv()
{
    // RngSeed(1231545);
    var k = 20000;
    var M = 1000;
    var E = 6;
    var O = 14;
    var lt = k.Range().Select(_ => GenBR((-M * M, M * M), (-2 * E, 2 * E), O)).ToArray();
    GlobalStopWatch.AddLap();
    var err = Double.Pow(10, -O + 3);
    foreach (var a in lt)
    {
        var a_b = a.Inv();
        var a_b2 = 1.0 / a.ToDouble;
        var diff = Double.Abs(a_b.ToDouble - a_b2) * Double.Pow(10, -a_b.V);
        if ((!a_b.IsZero() && diff > err) || (a_b.IsZero() && Double.Abs(a_b2) > err))
        {
            Console.WriteLine(
                new
                {
                    a, a_b.Details, a_b, a_bd = $"{fmt(a_b.ToDouble, O)}", a_bd2 = $"{fmt(a_b2, O)}", diff,
                    epsilon = err
                });
            Console.WriteLine();
        }
    }

    GlobalStopWatch.Show("END BigRealInv");
}

void BigRealRound()
{
    // RngSeed(1231545);
    var k = 2000;
    var M = 1000;
    var E = 6;
    var O = 14;
    var lt = k.Range().Select(_ => GenBR((-M * M, M * M), (-2 * E, 2 * E), O))
        .Concat(k.Range().Select(_ => GenBR((-M * M, M * M), (0, 0), O))).ToArray();
    GlobalStopWatch.AddLap();
    foreach (var a in lt)
    {
        for (int i = 0; i <= O / 2; i++)
        {
            var err = Double.Pow(10, -i);

            var a_b = BigReal.Round(a, i);
            var a_b2 = Double.Round(a.ToDouble, i, MidpointRounding.ToEven);
            var diff = Double.Abs(a_b.ToDouble - a_b2) * Double.Pow(10, -a_b.V);
            if ((!a_b.IsZero() && diff > err) || (a_b.IsZero() && Double.Abs(a_b2) > err))
            {
                Console.WriteLine(
                    new
                    {
                        i, a = a.Details, ad = $"{fmt(a.ToDouble, O)}", a_b.Details, a_b,
                        a_bd = $"{fmt(a_b.ToDouble, O)}", a_bd2 = $"{fmt(a_b2, O)}", diff,
                        epsilon = err
                    });
                Console.WriteLine();
            }
        }
    }

    GlobalStopWatch.Show("END BigRealRound");
}

// From AECF page 362
void LLL_Application_Pi()
{
    var n = 20;
    var N = new Rational(BigInteger.Pow(10, n));
    // pi=3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117067
    var pi = new Rational(BigInteger.Parse("314159265358979323846264338327950"), BigInteger.Pow(10, 32));

    var d = 8;
    var Nalpha = new Rational(BigInteger.Parse("-1669947371922907049619"));
    var alpha = Nalpha / N;

    Console.WriteLine(pi);
    Console.WriteLine($"{Nalpha.Num / Nalpha.Denom}");

    var mat = new KMatrix<Rational>(Rational.KZero(), d, d).One;
    for (int i = 0; i < d; i++)
    {
        var a = pi.Pow(i) * N;
        mat.Coefs[i, d - 1] = new Rational(a.Num / a.Denom);
    }

    mat.Coefs[d - 1, d - 1] = new Rational(Nalpha.Num / Nalpha.Denom);
    // mat.Coefs[6, d - 1] = new Rational(BigInteger.Parse("96138919357530443703022"));

    Console.WriteLine(mat);
    var lll = IntFactorisation.LLL(mat.T);
    Console.WriteLine();
    Console.WriteLine(lll);
    var l0 = lll.GetCol(0).ToArray();

    var N2 = N * N;
    var cols = lll.Cols.Where(l => l.Aggregate(Rational.KZero(), (acc, v) => acc + v.Pow(2)).CompareTo(N.One * 200 * 200) == -1)
        .ToArray();
    cols.Select(l => l.T).Println();

    Console.WriteLine(l0.Glue("; "));
    var sum = l0.SkipLast(1).Select((v, i) => (i, v)).Aggregate(pi.Zero, (acc, c) => acc + c.v * pi.Pow(c.i));
    Console.WriteLine(sum);
    Console.WriteLine((0.0 + sum) / (0.0 + alpha));
    Console.WriteLine(sum / alpha);
    Console.WriteLine($"{(double)alpha}");

    var alpha0 = (5 * pi.Pow(2) / 24) * (3 * pi.Pow(4) - 28 * pi.Pow(2) - 24);
    var Nalpha0 = N * alpha0;

    Console.WriteLine();
    Console.WriteLine(alpha0);
    Console.WriteLine($"{Nalpha0.Num / Nalpha0.Denom}");
}

// {
//     // LLL_Application_Pi();
//
//     GlobalStopWatch.Restart();
//     GlobalStopWatch.AddLap();
//
//     BigRealBInt();
//     BigRealStr();
//     BigRealAdd();
//     BigRealMul();
//     BigRealKMul();
//     BigRealDiv();
//     BigRealInv();
//     BigRealRat();
//     BigRealRound();
//
//     GlobalStopWatch.Show("END");
// }

void test()
{
    LLL_Application_Pi();
    var x = FG.QPoly();
    var P = x.Pow(4) - 7 * x.Pow(2) + 10;
    FG.NRoots(P.ToCPoly()).Println();

    var Px = P.ToBcPoly();
    var roots = FG.NRoots(Px);
    roots.Order().Println();
    var prod = roots.Aggregate(Px.One, (acc, r) => acc * (Px.X - r));
    Console.WriteLine(prod);
    Console.WriteLine(prod.ToRatPoly());

    var sqrt2 = BigCplx.NthRoot(Px.KOne * 2, 2);
    var sqrt5 = BigCplx.NthRoot(Px.KOne * 5, 2);
    var x150 = BigCplx.NthRoot(Px.KOne * 150, 4);
    Console.WriteLine(sqrt2.ToString());
    Console.WriteLine(sqrt5.ToString());
    Console.WriteLine(x150.ToString());
    Console.WriteLine(x150.Pow(4).ToString());

    var i = BigCplx.BgI();
    var a = 1 + i;
    Console.WriteLine(a);
    Console.WriteLine(a + a);
    Console.WriteLine(a * 5);
    Console.WriteLine(BigCplx.MagnitudeBigReal(a));
    var b = BigCplx.NthRoot(-i.One * 150, 4);
    Console.WriteLine(b);
    Console.WriteLine(b.Pow(4));
}

void bench()
{
    // GlobalStopWatch.Bench(5, "Example3 S3", AlgebraicIntegerRelationLLL.Example3);
    // GlobalStopWatch.Bench(5, "Example3 S3", AlgebraicIntegerRelationLLL.Example3);
    GlobalStopWatch.Bench(5, "Example4 D8", AlgebraicIntegerRelationLLL.Example4);
    // GlobalStopWatch.Bench(5, "Example4 D8", AlgebraicIntegerRelationLLL.Example4);
    // GlobalStopWatch.Bench(5, "Example5 D10", AlgebraicIntegerRelationLLL.Example5);
    // GlobalStopWatch.Bench(5, "Example5 D10", AlgebraicIntegerRelationLLL.Example5);
    // GlobalStopWatch.Bench(5, "Example6 A4", AlgebraicIntegerRelationLLL.Example6);
    // GlobalStopWatch.Bench(5, "Example6 A4", AlgebraicIntegerRelationLLL.Example6);
}

void Quintic1()
{
    var x = FG.QPoly();
    var O1 = 120;
    var O2 = 140;
    var P = x.Pow(5) + x.Pow(4) + 1;
    var (f0, f1) = IntFactorisation.FirrZ(P, details: true).Order().Deconstruct();
    var (_, X0, a0, a1) = IntFactorisation.PrimitiveElt(f0, f1)[0];
    Console.WriteLine(P.Substitute(a0));
    Console.WriteLine(P.Substitute(a1));
    var facts = IntFactorisation.AlgebraicFactors(P.Substitute(X0), details: true);
    Console.WriteLine(facts.Count(e => e.Degree > 1));
    var (minPoly1, _, _) = IntFactorisation.PrimitiveElt(facts.First(e => e.Degree > 1));
    IntFactorisation.FirrZ(minPoly1, details: true);
    
    var gal = AlgebraicIntegerRelationLLL.GaloisGroupLLL(minPoly1, O1, O2);
    DisplayGroup.HeadElements(gal);
    DisplayGroup.AreIsomorphics(gal, Product.Generate(new Cn(2), new Symm(3)));
    
    var subFields = GaloisTheory.SubFields(gal.Select(ka => ka.E).ToList(), nbGens: 3).ToArray();
    var extTowers = GaloisApplications.ExtensionsTower(subFields);
    GaloisApplications.GaloisCorrespondence(extTowers);

    var (X1, y1) = FG.EPolyXc(minPoly1, 'y');
    var (r0, r1) = IntFactorisation.AlgebraicRoots(f0.Substitute(X1), details: true).Deconstruct();
    var (r2, r3, r4) = IntFactorisation.AlgebraicRoots(f1.Substitute(X1), details: true).Deconstruct();
    
    GaloisApplications.FindExtension(subFields, r0);
    GaloisApplications.FindExtension(subFields, r1);
    GaloisApplications.FindExtension(subFields, r2);
    GaloisApplications.FindExtension(subFields, r3);
    GaloisApplications.FindExtension(subFields, r4);
    
    Console.WriteLine(gal.Aggregate(X1.One, (acc, r) => acc * (X1 - r)));
}

void Cyclo11_1()
{
    Ring.DisplayPolynomial = MonomDisplay.StarCaret;
    var x = FG.QPoly();
    var P0 = FG.CyclotomicPolynomial(11);
    var (X0, e0) = FG.EPolyXc(P0, 'e');
    var (minPol, e1, i0) = IntFactorisation.PrimitiveElt(X0.Pow(2) + 1);
    Console.WriteLine(minPol);
    IntFactorisation.FirrZ2(minPol, true);
    var (X, y) = FG.EPolyXc(minPol, 'y');
    var e = e1.Substitute(y);
    var i = i0.Substitute(y);
    Console.WriteLine("###################");
    Console.WriteLine((X - i) * (X + i));
    
    var p0Roots = 10.Range(1).Select(k => e.Pow(k)).ToArray();
    Console.WriteLine(P0);
    Console.WriteLine(p0Roots.Aggregate(X.One, (acc, ri) => acc * (X - ri)));
    
    var y0k = p0Roots.Select(ek => -ek + i).ToArray();
    var y1k = p0Roots.Select(ek => -ek - i).ToArray();
    var roots = y0k.Union(y1k).ToList();
    Console.WriteLine(minPol);
    Console.WriteLine(roots.Aggregate(X.One, (acc, ri) => acc * (X - ri)));
    
    var gal = GaloisTheory.GaloisGroup(roots, details: true);
    DisplayGroup.AreIsomorphics(FG.Abelian(2, 10), gal);
    DisplayGroup.AreIsomorphics(FG.Abelian(2, 2, 5), gal);
    
    var cos = (e + 1 / e) / 2;
    var sin = (e - 1 / e) / (2 * i);
    Console.WriteLine(cos);
    Console.WriteLine(sin);
    Console.WriteLine(cos.Pow(2) + sin.Pow(2));
    var P1 = IntFactorisation.MinPolynomial(X - cos).Substitute(x / 2).Monic;
    var P2 = IntFactorisation.MinPolynomial(X - sin).Substitute(x / 2).Monic;
    Console.WriteLine(P1);
    Console.WriteLine(P2);
    // IntFactorisation.SplittingField(P1, true);
}

void Cyclo11_2()
{
    Ring.DisplayPolynomial = MonomDisplay.StarCaret;
    var pRoots = new NthRootQ(44).PrimitivesRoots();
    var gal = Group.KAut(pRoots);
    DisplayGroup.HeadElements(gal);
    DisplayGroup.AreIsomorphics(FG.Abelian(2, 10), gal);
    var y = gal.Neutral().E;
    var X = FG.KPoly('X', y);
    
    var i = y.Pow(11);
    var ξ11 = y.Pow(4);
    
    var subFields = GaloisApplications.GaloisCorrespondence(pRoots.ToList());
    GaloisApplications.FindExtension(subFields, i.Substitute(gal.Neutral()), "Q(i)");
    GaloisApplications.FindExtension(subFields, ξ11.Substitute(gal.Neutral()), "Q(ξ11)");
    
    Console.WriteLine(ξ11);
    Console.WriteLine(i);
    Console.WriteLine((X - i) * (X + i));
    Console.WriteLine(10.Range(1).Aggregate(X.One, (acc, k) => acc * (X - ξ11.Pow(k))));
    
    var cos = (ξ11 + 1 / ξ11) / 2;
    var sin = (ξ11 - 1 / ξ11) / (2 * i);
    Console.WriteLine(cos);
    Console.WriteLine(GaloisTheory.Rewrite(ξ11, cos));
    Console.WriteLine(sin);
    Console.WriteLine(cos.Pow(2) + sin.Pow(2));
    var exp_2ipi_11 = new Cnf(11);
    Console.WriteLine(exp_2ipi_11.Re);
    Console.WriteLine(exp_2ipi_11.Im);
    Console.WriteLine(exp_2ipi_11.Re.Pow(2) + exp_2ipi_11.Im.Pow(2));
}

{
    // Quintic1();
    // Cyclo11_1();
    // Cyclo11_2();
    AlgebraicIntegerRelationLLL.Example6();
}
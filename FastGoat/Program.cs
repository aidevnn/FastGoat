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

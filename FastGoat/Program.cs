using System.Numerics;
using System.Reflection;
using FastGoat.Commons;
using FastGoat.Examples;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.LWE;
using FastGoat.UserGroup.Matrix;
using FastGoat.UserGroup.Padic;
using FastGoat.UserGroup.Polynoms;
using Rq = FastGoat.Structures.VecSpace.KPoly<FastGoat.UserGroup.Integers.Rational>;
using static FastGoat.Commons.IntExt;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");
Ring.DisplayPolynomial = MonomDisplay.StarCaret;

int Phi2(int n)
{
    var (num, denom) = (1, 1);
    foreach (var pi in PrimesDec(n).Keys)
    {
        num *= pi - 1;
        denom *= pi;
    }

    return n / denom * num;
}

int Phi3(int n) => (int)double.Round(n * PrimesDec(n).Keys.Aggregate(1.0, (acc, pi) => acc * (1.0 - 1.0 / pi)));

IEnumerable<int> Dividors2(int n)
{
    for (int i = 1; i <= double.Sqrt(n); i++)
    {
        if (n % i == 0)
        {
            yield return i;
            var i2 = n / i;
            if (i != i2)
                yield return i2;
        }
    }
}

IEnumerable<int> Dividors3(int n)
{
    var decomp = PrimesDecomposition(int.Abs(n)).GroupBy(e => e)
        .ToDictionary(e => e.Key, e => e.Count());

    return decomp.Select(e => (e.Value + 1).Range().Select(k => (int)double.Round(double.Pow(e.Key, k))))
        .MultiLoop()
        .Select(l => l.Aggregate(1, (acc, e) => acc * e));
}

// {
//     for (int n = 1000000; n <= 1100000; n += 10000)
//     {
//         var phi1 = Phi(n);
//         var phi2 = Phi2(n);
//         var phi3 = Phi3(n);
//         Console.WriteLine($"n={n,3} Phi1={phi1,3} Phi2={phi2,3} Phi3={phi3,3}");
//         if (phi1 != phi2 || phi1 != phi3)
//             throw new();
//
//         var n0 = n;
//         GlobalStopWatch.Bench(50, "Phi1", () => Phi(n0));
//         GlobalStopWatch.Bench(50, "Phi2", () => Phi2(n0));
//         GlobalStopWatch.Bench(50, "Phi3", () => Phi3(n0));
//         Console.WriteLine();
//     }
// }
//
// {
//     for (int n = 10000000; n <= 11000000; n += 100000)
//     {
//         var divs1 = Dividors(n).Append(n).ToArray();
//         var divs2 = Dividors2(n).Order().ToArray();
//         var divs3 = Dividors3(n).Order().ToArray();
//         Console.WriteLine($"n={n,3}");
//         Console.WriteLine($"Divs1={divs1.Length} Divs2={divs2.Length} Divs3={divs3.Length}");
//         if (!divs1.SequenceEqual(divs2) || !divs1.SequenceEqual(divs3))
//             throw new();
//
//         var n0 = n;
//         GlobalStopWatch.Bench(50, "Divs1", () => Dividors(n0).Append(n0).ToArray());
//         GlobalStopWatch.Bench(50, "Divs2", () => Dividors2(n0).Append(n0).ToArray());
//         GlobalStopWatch.Bench(50, "Divs3", () => Dividors3(n0).Append(n0).ToArray());
//         Console.WriteLine();
//     }
// }

// {
//     EllipticCurves.Example5FromLMFDB();
//     EllipticCurves.Example5FromLMFDB();
//     EllipticCurves.Example5FromLMFDB();
//     EllipticCurves.Example5FromLMFDB();
// }

void testPowFast(int p, int n, int nbTests = 100)
{
    var q = p.Pow(n);
    var x = FG.FqX(q);
    for (int i = 0; i < nbTests; i++)
    {
        var k = Rng.Next(2, q);
        var xk = x.Pow(k);
        var m = Rng.Next(2, q);
        var xkm1 = xk.Pow(m);
        var xkm2 = xk.FastPow(m);
        if (!xkm1.Equals(xkm2))
        {
            Console.WriteLine(new { i, n, p, k, m, xkm1, xkm2 });
            throw new();
        }
    }

    Console.WriteLine($"n={n} p={p} q={q} Pass {nbTests} Tests");
}

void benchPowFast(int p, int n, int m, int nbTests = 3)
{
    var q = p.Pow(n);
    var x = FG.FqX(q);
    for (int i = 0; i < nbTests; i++)
    {
        var k = Rng.Next(2, q);
        var xk = x.Pow(k);
        GlobalStopWatch.Bench(50, $"Pow  p={p} n={n} m={m}", () => xk.Pow(m));
        GlobalStopWatch.Bench(50, $"Pow2 p={p} n={n} m={m}", () => xk.FastPow(m));
    }

    Console.WriteLine();
}

void runPowFast()
{
    var max_q = 1 << 10;
    foreach (var p in IntExt.Primes10000.Where(p => p * p < max_q))
    {
        var max_n = (int)(Double.Log(max_q) / Double.Log(p));
        for (int n = 1; n <= max_n; n++)
        {
            testPowFast(p, n);
            benchPowFast(p, n, 500 * (12 - n));
            Console.WriteLine();
        }
    }
}

{
    var q = 3.Pow(4);
    var x = FG.FqX(q);
    var set = new List<EPoly<ZnInt>>();
    var divs = Dividors(q - 1).ToArray();
    for (int i = 1; i < q - 1; i++)
    {
        var a = x.FastPow(i);
        if (divs.All(d => !a.FastPow(d).IsOne()))
            set.Add(a);
    }
    
    var gf = new GFp($"GF({x.F})", x);
    var primElts = Group.Generate(gf).ElementsOrders.Where(e => e.Value == q - 1).Select(e => e.Key).Order().ToArray();
    set = set.Order().ToList();
    set.Order().Println($"Seqs equal {set.SequenceEqual(primElts)} Count={set.Count}");
}
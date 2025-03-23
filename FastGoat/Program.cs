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

{
    for (int n = 1000000; n <= 1100000; n += 10000)
    {
        var phi1 = Phi(n);
        var phi2 = Phi2(n);
        var phi3 = Phi3(n);
        Console.WriteLine($"n={n,3} Phi1={phi1,3} Phi2={phi2,3} Phi3={phi3,3}");
        if (phi1 != phi2 || phi1 != phi3)
            throw new();

        var n0 = n;
        GlobalStopWatch.Bench(50, "Phi1", () => Phi(n0));
        GlobalStopWatch.Bench(50, "Phi2", () => Phi2(n0));
        GlobalStopWatch.Bench(50, "Phi3", () => Phi3(n0));
        Console.WriteLine();
    }
}

{
    for (int n = 10000000; n <= 11000000; n += 100000)
    {
        var divs1 = Dividors(n).Append(n).ToArray();
        var divs2 = Dividors2(n).Order().ToArray();
        var divs3 = Dividors3(n).Order().ToArray();
        Console.WriteLine($"n={n,3}");
        Console.WriteLine($"Divs1={divs1.Length} Divs2={divs2.Length} Divs3={divs3.Length}");
        if (!divs1.SequenceEqual(divs2) || !divs1.SequenceEqual(divs3))
            throw new();

        var n0 = n;
        GlobalStopWatch.Bench(50, "Divs1", () => Dividors(n0).Append(n0).ToArray());
        GlobalStopWatch.Bench(50, "Divs2", () => Dividors2(n0).Append(n0).ToArray());
        GlobalStopWatch.Bench(50, "Divs3", () => Dividors3(n0).Append(n0).ToArray());
        Console.WriteLine();
    }
}
using System.Numerics;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Lattice;
using Rq = FastGoat.Structures.VecSpace.KPoly<FastGoat.UserGroup.Integers.Rational>;
using static FastGoat.Commons.IntExt;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");

Dictionary<int, double> DGS_CDF(double sigma, double tau)
{
    var bound = (int)double.Floor(sigma * tau);
    var tlength = 2 * bound + 1;
    var sigma2 = 2 * sigma * sigma;
    var norm = tlength.SeqLazy(-bound).Sum(x => double.Exp(-x * x / sigma2));
    var pdf = tlength.SeqLazy(-bound).ToDictionary(x => x, x => double.Exp(-x * x / sigma2) / norm);
    return tlength.SeqLazy(1).ToDictionary(i => i, i => i.SeqLazy(-bound).Sum(e => pdf[e]));
}

void BatchTest(int sizeI, double sigma = 3.2, double tau = 1.0)
{
    var bound = (int)double.Floor(sigma * tau);
    var tlength = 2 * bound + 1;
    var digits = $"{4 * sizeI / (1 + bound)}".Length;
    var fmt = $"[{{0,2}}]:{{1,{digits}}}";
    var cdf = DGS_CDF(sigma, tau);
    var size = 0;
    var table = tlength.SeqLazy(-bound).ToDictionary(i => i, _ => 0);
    double max, test;
    var fix = false;
    do
    {
        size += sizeI;
        foreach (var i in DistributionExt.DiscreteGaussianSample(sizeI, sigma, tau))
            table[i]++;

        var cdfTest = tlength.SeqLazy(1).ToDictionary(i => i, i => i.SeqLazy(-bound).Sum(e => table[e]) / (1.0 * size));
        max = cdfTest.Max(e => double.Abs(e.Value - cdf[e.Key]));
        test = 1.358101 / double.Sqrt(size);

        fix |= max > test;
    } while (max > test);

    var sizeInfo = fix ? $"size:{sizeI} => {size} (Fixed x{size / sizeI})" : $"size:{size}";
    Console.WriteLine($"DiscreteGaussian Sample {sizeInfo} Sigma:{sigma} Tau:{tau} Bounds:[{-bound}, {bound}]");
    Console.WriteLine(
        $"results    : {{ {table.OrderBy(e => e.Key).Select(e => string.Format(fmt, e.Key, e.Value)).Glue(", ")} }}");
    Console.WriteLine($"Discrete Kolmogorov-Smirnov:{max:F9} Test:{test:F6} => Accepted alpha=0.05");
    Console.WriteLine();
}

{
    var size = 5000;
    var sigmas = new[] { 0.4, 0.5, 1.0, 1.5, 2, 3, 5, 7 };
    var taus = new[] { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 };
    foreach (var (sigma, tau) in sigmas.Grid2D(taus))
        BatchTest(size, sigma, tau);
}
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
// IntExt.RngSeed(7532159);
Ring.DisplayPolynomial = MonomDisplay.StarCaret;

{
    var size = 1000000;
    var sigma = 7.0;
    var sample = DistributionExt.DiscreteGaussianSample(size, sigma);
    var table = sample.GroupBy(e => e).ToDictionary(e => e.Key, e => e.Count());

    // D = DiscreteGaussianDistributionIntegerSampler(sigma=7,tau=1)
    // defaultdict(<class 'sage.rings.integer.Integer'>, {0: 79153, 7: 47886, 5: 61306, 4: 67843, -4: 67731, -1: 79087, 6: 55064, -3: 72199, -2: 76615, -5: 61623, 3: 73368, 2: 76455, -7: 47507, 1: 78678, -6: 55485})
    var sageCounter =
        "0: 79153, 7: 47886, 5: 61306, 4: 67843, -4: 67731, -1: 79087, 6: 55064, -3: 72199, -2: 76615, -5: 61623, 3: 73368, 2: 76455, -7: 47507, 1: 78678, -6: 55485";
    var table2 = sageCounter.Split(',').Select(e => e.Split(":").Select(i => int.Parse(i)).ToArray())
        .ToDictionary(e => e[0], e => e[1]);

    var min = table.Keys.Min();
    var tlength = table.Count;
    var cdf1 = tlength.SeqLazy(1).ToDictionary(i => i, i => i.SeqLazy(min).Sum(e => table[e]) / (1.0 * size));
    var cdf2 = tlength.SeqLazy(1).ToDictionary(i => i, i => i.SeqLazy(min).Sum(e => table2[e]) / (1.0 * size));
    var max = cdf1.Max(e => double.Abs(e.Value - cdf2[e.Key]));
    var test = 1.36 / double.Sqrt(size);

    Console.WriteLine($"DiscreteGaussian Sample size:{size} Sigma:{sigma} Tau:1");
    Console.WriteLine(
        $"results    : {{ {table.OrderBy(e => e.Key).Select(e => $"{e.Key,2} =>{e.Value,6}").Glue(", ")} }}");
    Console.WriteLine($"sage       :{sageCounter}");
    Console.WriteLine($"max diff   :{max:E3} test:{test:E3}");
    Console.WriteLine($"Accepted 5%:{test > max}");
    // DiscreteGaussian Sample size:1000000 Sigma:7 Tau:1
    // results    : { -7 => 48238, -6 => 55333, -5 => 61523, -4 => 67234, -3 => 72717, -2 => 76066, -1 => 78844,  0 => 79200,  1 => 78609,  2 => 76282,  3 => 72884,  4 => 67815,  5 => 61950,  6 => 54689,  7 => 48616 }
    // sage       :0: 79153, 7: 47886, 5: 61306, 4: 67843, -4: 67731, -1: 79087, 6: 55064, -3: 72199, -2: 76615, -5: 61623, 3: 73368, 2: 76455, -7: 47507, 1: 78678, -6: 55485
    // max diff   :9.990E-004 test:1.360E-003
    // Accepted 5%:True
}
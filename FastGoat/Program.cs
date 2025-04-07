using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Integers;
using static FastGoat.Commons.IntExt;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");

void BSGSonAbelianGroup(int nbTest = 50)
{
    var z = FG.Abelian([Primes10000[Rng.Next(5)], Primes10000[Rng.Next(5)], Primes10000[Rng.Next(5, 10)]]);
    var n = z.Count();
    DisplayGroup.HeadOrders(z);
    for (int i = 0; i < nbTest; i++)
    {
        var a = z.ElementAt(Rng.Next(1, n));
        var oa = z.ElementsOrders[a];
        var k0 = Rng.Next(oa / 4, oa) + 1;
        var b = z.Times(a, k0);
        var k1 = Group.BSGS(z, a, b, n);
        Console.WriteLine("{0,5} x {1} = {2}", k1, a, b);
        if (k0 != k1)
            throw new($"k0={k0} k1={k1} b = {b} != {z.Times(a, k1)}");
    }

    Console.WriteLine();
}

void BenchBSGS<T>(ConcreteGroup<T> g) where T : struct, IElt<T>
{
    DisplayGroup.HeadOrders(g);
    var g1g2Arr = g.ToArray();
    var eltsOrdBSGS = Group.ElementsOrdersBSGS(g, g1g2Arr, g1g2Arr.Length);
    Console.WriteLine("Elements Orders : {0} BSGS",
        eltsOrdBSGS.Values.GroupBy(a => a).ToDictionary(a => a.Key, a => a.Count()).AscendingByKey()
            .GlueMap(fmt: "[{0}]:{1}"));
    
    if (g.ElementsOrders.Any(e => eltsOrdBSGS[e.Key] != e.Value))
        throw new();

    Console.WriteLine();

    GlobalStopWatch.Bench(5, "Old ", () => Group.ElementsOrders(g, g1g2Arr));
    GlobalStopWatch.Bench(5, "BSGS", () => Group.ElementsOrdersBSGS(g, g1g2Arr, g1g2Arr.Length));
    Console.WriteLine();
}

{
    BSGSonAbelianGroup();
    BSGSonAbelianGroup();
    
    BenchBSGS(new Cn(908));
    BenchBSGS(Product.Generate(new Cn(23), new Cn(61)));
    BenchBSGS(Product.Generate(FG.Dihedral(4), FG.Dihedral(8), FG.Dihedral(10)));
    BenchBSGS(Product.Generate(new Cn(41), FG.GL2p(5)));
}
// |C23 x C61| = 1403
// Type        AbelianGroup
// BaseGroup   C23 x C61
// 
// Elements Orders : [1]:1, [23]:22, [61]:60, [1403]:1320
// 
// Elements Orders : [1]:1, [23]:22, [61]:60, [1403]:1320 BSGS
// 
// # Old  Avg Time:169 ms Dev:17.938
// # BSGS Avg Time:9 ms Dev:1.720
// 
// |C41 x GL(2,5)| = 19680
// Type        NonAbelianGroup
// BaseGroup   C41 x GL(2,5)
// 
// Elements Orders : [1]:1, [2]:31, [3]:20, [4]:152, [5]:24, [6]:20, [8]:40, [10]:24, [12]:40, [20]:48, [24]:80, [41]:40, [82]:1240, [123]:800, [164]:6080, [205]:960, [246]:800, [328]:1600, [410]:960, [492]:1600, [820]:1920, [984]:3200
// 
// Elements Orders : [1]:1, [2]:31, [3]:20, [4]:152, [5]:24, [6]:20, [8]:40, [10]:24, [12]:40, [20]:48, [24]:80, [41]:40, [82]:1240, [123]:800, [164]:6080, [205]:960, [246]:800, [328]:1600, [410]:960, [492]:1600, [820]:1920, [984]:3200 BSGS
// 
// # Old  Avg Time:2278 ms Dev:24.244
// # BSGS Avg Time:458 ms Dev:7.414
// 
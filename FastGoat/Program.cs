using FastGoat.Commons;
using FastGoat.Examples;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.EllCurve;
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

void testsBSGS()
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

// Elliptic Curve Discrete Logarithm Problem
void ECDLP(int p, (int x, int y) P0, (int x, int y) Q0, int[] curve)
{
    var (Poly, trans, rev) = EllipticCurves.MinimizedFormZnInt(p, curve);
    var (A, B) = (Poly[1], Poly[0]);
    var E = new EllGroup<ZnInt>(A, B);
    var nb = EllipticCurves.SchoofEllPtsCount(A.K, B.K, p);
    Console.WriteLine($"|{E}| = {nb}");

    var P = new EllPt<ZnInt>(new(p, P0.x), new(p, P0.y));
    var Q = new EllPt<ZnInt>(new(p, Q0.x), new(p, Q0.y));

    var (P1, Q1) = (trans(P), trans(Q));
    var k = Group.BSGS(E, P1, Q1, nb);
    var Q2 = E.Times(P1, k);
    Console.WriteLine($"P1={P1} Q1={Q1} {k}xP1=Q1");
    Console.WriteLine($"{k} x {P} = {Q}");
    Console.WriteLine();
    if (!Q2.Equals(Q1))
        throw new();
}

{
    Ring.DisplayPolynomial = MonomDisplay.StarCaret;
    Logger.Level = LogLevel.Level1;
    
    // C : y^2 = x^3 + x^2 + x + 1 p = 97, P = (7, 20), Q = (17, 46)
    // sage: E=EllipticCurve(GF(97),[0,1,0,1,1]);P=E([7, 20]);Q=Ea([17, 46]);Q.log(P)
    // pari: E=ellinit([0,1,0,1,1],97);elllog(E,[17, 46],[7, 20])
    ECDLP(97, (7, 20), (17, 46), [0, 1, 0, 1, 1]);
    
    // (a) C : y^2 = x^3 + x^2 + x + 3, p = 103, P = (7, 14), Q = (8, 22).
    ECDLP(103, (7, 14), (8, 22), [0, 1, 0, 1, 3]);

    // (b) C : y^2 = x^3 − 2x^2 + 5x + 6, p = 149, P = (11, 16), Q = (110, 46).
    ECDLP(149, (11, 16), (110, 46), [0, -2, 0, 5, 6]);

    // (c) C : y^2 = x^3 + x^2 + x + 2, p = 10037, P = (8, 7358), Q = (2057, 5437).
    ECDLP(10037, (8, 7358), (2057, 5437), [0, 1, 0, 1, 2]);
}
// Elliptic curve      y^2 = x^3 + x^2 + x + 1 in Z/97Z
// Simplified form     y^2 = x^3 + 33*x + 69
// 
// |Ell[33,69](Z/97Z)| = 100
// P1=(72,77) Q1=(82,51) 47xP1=Q1
// 47 x ( 7,20) = (17,46)
// 
// Elliptic curve      y^2 = x^3 + x^2 + x +   3 in Z/103Z
// Simplified form     y^2 = x^3 +  35*x +  18
// 
// |Ell[ 35, 18](Z/103Z)| = 109
// P1=( 76, 89) Q1=( 77, 81) 76xP1=Q1
// 76 x (  7, 14) = (  8, 22)
// 
// Elliptic curve      y^2 = x^3 + 147*x^2 +   5*x +   6 in Z/149Z
// Simplified form     y^2 = x^3 + 103*x +  86
// 
// |Ell[103, 86](Z/149Z)| = 169
// P1=( 60,133) Q1=( 10,103) 87xP1=Q1
// 87 x ( 11, 16) = (110, 46)
// 
// Elliptic curve      y^2 = x^3 + x^2 + x +     2 in Z/10037Z
// Simplified form     y^2 = x^3 +  6692*x +  9667
// 
// |Ell[ 6692, 9667](Z/10037Z)| = 10151
// P1=( 3354, 2679) Q1=( 5403, 4600) 1277xP1=Q1
// 1277 x (    8, 7358) = ( 2057, 5437)
// 
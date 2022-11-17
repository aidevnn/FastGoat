using System.Collections;
using System.Diagnostics;
using FastGoat.Commons;
using FastGoat.Examples;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Matrix;
using FastGoat.UserGroup.Perms;
using FastGoat.UserGroup.Polynoms;
using FastGoat.UserGroup.Words;
using FastGoat.UserGroup.Words.ToddCoxeter;
using static FastGoat.Commons.IntExt;
using static FastGoat.Commons.EnumerableExt;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");

void Possibilities(int n)
{
    int[][] Part(int p, int pow)
    {
        var parts = Partitions32[pow];
        return parts.Select(i => i.Select(i0 => (int)Math.Pow(p, i0)).ToArray()).ToArray();
    }

    var decomp = PrimesDecomposition(n).GroupBy(i => i).Select(e => Part(e.Key, e.Count())).ToArray();
    List<int[]> all = new();
    foreach (var l in decomp.MultiLoop())
    {
        var l1 = l.SelectMany(e => e).Ascending().ToArray();
        all.Add(l1);
    }

    foreach (var l3 in all.OrderBy(l=>l.Length).ThenBy(l => l, Comparer<int[]>.Create((e1, e2) => e1.SequenceCompareTo(e2))))
        Console.WriteLine($"{n,-3} = {l3.Glue(" x ")}");
}

void Test128()
{
    for (int i = 2; i <= 128; i++)
    {
        Possibilities(i);
        Console.WriteLine();
    }
}

void ES1(int p)
{
    if (!Primes10000.Contains(p)) 
        throw new GroupException(GroupExceptionType.GroupDef);
    
    var wg = new WordGroup($"a{p*p}, b{p}, ab = ba{p+1}");
    DisplayGroup.HeadElements(wg);
}

void ES2(int p)
{
    if (!Primes10000.Contains(p)) 
        throw new GroupException(GroupExceptionType.GroupDef);
    
    var wg = new WordGroup($"a{p}, b{p}, c{p}, c = a-1b-1ab, ac = ca, bc = cb");
    DisplayGroup.HeadElements(wg);
}

void Suzuki8()
{
    var gl48 = new GLnq(4, 8);
    var x = gl48.Fq['x'];
    var a = gl48[
        1, 0, 0, 0,
        1, 1, 0, 0,
        x + 1, 1, 1, 0,
        x * x + 1, x, 1, 1
    ];
    var b = gl48[
        0, 0, 0, 1,
        0, 0, 1, 0,
        0, 1, 0, 0,
        1, 0, 0, 0
    ];

    var sz8 = Group.Generate("Sz(8)", gl48, a, b);
    DisplayGroup.Head(sz8);

    var n = sz8.Count();
    Console.WriteLine(PrimesDecomposition(n).GroupBy(i => i).ToDictionary(i => i.Key, i => i.Count())
        .AscendingByKey()
        .GlueMap(" . ", "{0}^{1}"));

    GlobalStopWatch.Time("Sz(8)", () => Group.Generate("Sz(8)", sz8, a, b));
}

void DetGenerators()
{
    var n = 3; // 2,3,4,5
    var alphabet = "abcdefghijklmnopqrstvwxyz";
    Monom.Display = MonomDisplay.StarPowFct;
    var coefs = Ring.Polynomial(ZnInt.KZero(), alphabet.Take(n * n).ToArray());
    var z0 = Ring.PolynomialZero(ZnInt.KZero());
    var mat = Ring.Matrix(n, coefs);
    var cmat =Ring.Transpose(Ring.CoMatrix(mat, z0));

    Console.WriteLine("Matrix");
    Ring.DisplayMatrix(mat);

    for (int i = 0; i < n; i++)
    {
        var i0 = i;
        var rg0 = n.Range().Select(j => mat[i0, j]).Glue(", ");
        var rg1 = n.Range().Select(j => $"mat[{n * i0 + j}]").Glue(", ");
        Console.WriteLine($"var ({rg0}) = ({rg1});");
    }

    Console.WriteLine();
    Console.WriteLine("var det = {0};", Ring.Determinant(mat, z0));
    Console.WriteLine("var idet = det.Inv();");
    Console.WriteLine();
    Console.WriteLine($"var inv = new EPoly<ZnInt>[{n * n}];");
    
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            Console.WriteLine("inv[{0}] = ({1}) * idet;", n * i + j, cmat[i, j]);
        }
    }

    Console.WriteLine();
    Console.WriteLine("return inv.Aggregate(0, (acc, a0) => a0.GetHashCode() + Fq.P * acc);");
}

{

}
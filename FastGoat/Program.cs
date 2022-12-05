using System.ComponentModel.DataAnnotations;
using System.Diagnostics;
using System.Linq.Expressions;
using FastGoat.Commons;
using FastGoat.Examples;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Matrix;
using FastGoat.UserGroup.Perms;
using FastGoat.UserGroup.Polynoms;
using FastGoat.UserGroup.Words;
using FastGoat.UserGroup.Words.ToddCoxeter;
using static FastGoat.Commons.IntExt;
using static FastGoat.Commons.EnumerableExt;
using System.Numerics;
using System.Threading.Channels;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");

Random rnd = new Random();

KPoly<K> ProdIrr<K>(KPoly<K> scalar, int p, int d)
    where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
{
    var x = scalar.X;
    return x.Pow(p.Pow(d)) - x;
}

bool IsIrreductibleFp<K>(KPoly<K> f)
    where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
{
    var p = f.P;
    if (p == 0)
        throw new ArgumentException();

    if (f.Degree < 1)
        return true;

    var n = f.Degree;
    var divs = IntExt.Dividors(n);
    var x = f.X;
    return ProdIrr(x, p, n).Div(f).rem.IsZero() && divs.All(d => !ProdIrr(x, p, d).Div(f).rem.IsZero());
}

KMatrix<K> RandMatrix<K>(int amp, int offset, int m, int n, K scalar)
    where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
{
    var A = new K[m, n];
    foreach (var (i, j) in m.Range().Grid2D(n.Range()))
        A[i, j] = scalar * rnd.Next(-amp, amp + 1) + rnd.Next(-offset, offset + 1);

    var mA = new KMatrix<K>(A);
    return mA;
}

KPoly<K> RandPolySep<K>(K scalar, int p, int n) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
{
    while (true)
    {
        var coefs = (n + 1).Range().Select(i => rnd.Next(-p, p + 1) * scalar.One).TrimSeq().ToArray();
        var f = new KPoly<K>('x', scalar, coefs);
        if (f.Degree > 1 && !Ring.Discriminant(f).IsZero())
            return f.Monic;
    }
}

EPoly<ZnInt> RandPoly(KPoly<ZnInt> f, int n)
{
    var p = f.P;
    var coefs = (n + 1).Range().Select(i => new ZnInt(p, rnd.Next(p))).ToArray();
    var g = new KPoly<ZnInt>(f.x, ZnInt.KZero(p), coefs);
    return new EPoly<ZnInt>(f, g.Div(f).rem);
}

List<(KPoly<K> g, int q, int i)> MusserSFF<K>(KPoly<K> f)
    where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
{
    var L = new List<(KPoly<K>, int, int)>();
    var c = Ring.Gcd(f, f.Derivative);
    var i = 1;
    var g = f / c;
    while (g.Degree >= 1)
    {
        var p = Ring.Gcd(c, g);
        c = c / p;
        if (g.Degree > p.Degree)
            L.Add(((g / p).Monic, 1, i));

        g = p;
        ++i;
    }

    return L;
}
//
// {
//     var x = FG.QFracPoly();
//     var amp = 0;
//     var offset = 2;
//     var m = 3;
//     var n = 3;
//     var A = RandMatrix(amp, offset, m, n, x);
//     Ring.MatrixDisplayForm = Ring.MatrixDisplay.Bracket;
//     Console.WriteLine(A);
//     
//     Ring.MatrixDisplayForm = Ring.MatrixDisplay.Table;
//     Console.WriteLine(A);
//     var det = A.Det;
//     Console.WriteLine(det);
//     
//     var ai = A.Inv();
//     Console.WriteLine(ai);
//
//     var M = A - x * A.One;
//     Console.WriteLine(M);
//     Console.WriteLine(M.Det);
//     Console.WriteLine(M.Inv());
//     Console.WriteLine(M.Pow(3));
// }

// {
//     var x = Rational.KZero();
//     var amp = 0;
//     var offset = 5;
//     var m = 5;
//     var n = 3;
//     // var A = RandMatrix(amp, offset, m, n, x);
//     var a = Ring.Matrix(3, Rational.KZero(),
//         1, 1, 2,
//         2, 2, 4,
//         3, 3, 6
//     );
//     var A = new KMatrix<Rational>(a);
//     Console.WriteLine(A);
//     var det = A.Det;
//     Console.WriteLine(det);
//     var (nt, ns) = A.NullSpace();
//     Console.WriteLine($"Nullity {nt}");
//     Console.WriteLine("NullSpace");
//     Console.WriteLine(ns.Select(e => e.Glue("\n", "[{0}]")).Glue("\n\n"));
// }

EPoly<K>[] CanonicalBase<K>(KPoly<K> f) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
{
    var x = new EPoly<K>(f);
    return f.Degree.Range().Select(i => x.Pow(i)).ToArray();
}

KMatrix<K> BerlekampMatrix<K>(KPoly<K> f, int q)
    where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
{
    var baseCan = CanonicalBase(f);
    var n = baseCan.Length;
    var M = new K[n, n];
    var polys = baseCan.Select(g => g.Pow(q) - g).ToArray();
    foreach (var (i, j) in n.Range().Grid2D(n.Range()))
    {
        M[i, j] = polys[j][i];
    }

    return new(M);
}

void TrivialFactor<K>(KPoly<K> f, int q)
    where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
{
    var baseCan = CanonicalBase(f);
    var bm = BerlekampMatrix(f, q);
    var (nt, ns) = bm.NullSpace();
    if (nt == 1)
    {
        Console.WriteLine("Irreductible");
        return;
    }

    var polys = ns.Select(c => c.Select((k, i) => k * baseCan[i]).Aggregate((a, b) => a + b)).ToArray();
    var fqa = (q - 1).Range(1).Select(i => f.One * i).ToArray();

    foreach (var g in polys.Select(g => g.Poly).Where(g => g.Degree > 0))
    {
        foreach (var a in fqa)
        {
            var g_a = g - a;
            if (Ring.Gcd(f, g_a).Degree != 0)
            {
                var gcd = Ring.Gcd(f, g_a).Monic;
                Console.WriteLine($"g - {a} = {g_a}");
                Console.WriteLine($"gcd = {gcd}");
                Console.WriteLine($"f/gcd = {(f / gcd)}");
                Console.WriteLine();
            }
        }
    }
}

{
    var p = 17;
    var n = 6;
    var f = RandPolySep(ZnInt.KZero(p), p, n);
    Console.WriteLine(f);
    Console.WriteLine($"Disc {Ring.Discriminant(f)}");
    Console.WriteLine(MusserSFF(f).Glue("\n"));
    var bm = BerlekampMatrix(f, p);
    Console.WriteLine(bm);

    var (nt, ns) = bm.NullSpace();
    Console.WriteLine($"Nullity {nt}");
    Console.WriteLine("NullSpace");
    Console.WriteLine(ns.Select(e => e.Glue("\n", "[{0}]")).Glue("\n\n"));
    TrivialFactor(f, p);
}
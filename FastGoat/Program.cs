using System.ComponentModel.DataAnnotations;
using System.Diagnostics;
using System.IO.IsolatedStorage;
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

int FirstPrime(KPoly<Rational> f)
{
    var disc = Ring.Discriminant(f);
    return Primes10000.First(p => disc.Num % p != 0);
}

double Nu(KPoly<Rational> f)
{
    var n = f.Degree;
    var norm = f.Coefs.Select(e => Double.Abs(e)).Sum();
    return Double.Sqrt(n + 1) * Double.Pow(2, n) * norm;
}

(int p, int s) PSigma(KPoly<Rational> f)
{
    var p = FirstPrime(f);
    var nu = Nu(f);
    var s = 1;
    while (double.Pow(p, s) <= 2 * nu)
    {
        ++s;
    }

    return (p, s);
}

KPoly<Padic> PadicPoly(int p, int o, char x = 'x') => new KPoly<Padic>(x, new Padic(p, o)).X;

KPoly<Rational> Padic2QPoly(KPoly<Padic> f) =>
    new(f.x, Rational.KZero(), f.Coefs.Select(e => new Rational(e.ToBigInt)).ToArray());

KPoly<Padic> QPoly2Padic(KPoly<Rational> f, int p, int o)
{
    var coefs = f.Coefs.Select(e => Padic.Convert(p, o, e.Num)).ToArray();
    return new(f.x, new Padic(p, o), coefs);
}

KPoly<ZnInt> QPoly2ZPoly(KPoly<Rational> f, int p)
{
    var coefs = f.Coefs.Select(e => new ZnInt(p, (int)e.Num)).ToArray();
    return new(f.x, ZnInt.KZero(p), coefs);
}

KPoly<Padic> Resize(KPoly<Padic> f, int o0) =>
    new(f.x, f.KZero.Resize(o0), f.Coefs.Select(e => e.Resize(o0)).ToArray());

{
    var x = FG.QPoly();
    // var f = x.Pow(4) - 1; // AECF example 21.2, page 387
    var f = x.Pow(4) + -x.Pow(3) + -x.Pow(2) + 2*x + 4;
    // var f = x.Pow(2) - 1;
    // var f = x.Pow(12) - 50 * x.Pow(10) + 753 * x.Pow(8) - 4520 * x.Pow(6) + 10528 * x.Pow(4) - 6720 * x.Pow(2) + 576;
    // var f = x.Pow(8) - 40 * x.Pow(6) + 352 * x.Pow(4) - 960 * x.Pow(2) + 576; // bug 
    var (p, o) = PSigma(f);
    
    var padicZero = new Padic(p, 1);
    var a0 = (new Un(p)).GetGenerators().First()[new(p, 1)];
    var f0 = QPoly2Padic(f, p, 1);
    Console.WriteLine(f);
    Console.WriteLine(f0);
    Console.WriteLine($"Prime P = {p}; Disc={Ring.Discriminant(f)}; Sigma = {o}");
    var firr0 = PolynomialFactorization.Firr(f0, a0 + padicZero);
    var all = new List<KPoly<Padic>>(firr0);
    var o0 = 1;
    while (o0 <= o)
    {
        Console.WriteLine();
        Console.WriteLine(all.Order().Glue("\n"));
        var tmp = new List<KPoly<Padic>>();
        o0 *= 2;
        var fa = QPoly2Padic(f, p, o0);
        foreach (var g in all)
        {
            var gi = Resize(g, o0);
            var dgi = new EPoly<Padic>(gi, gi.Derivative);
            var fi = new EPoly<Padic>(gi, fa.Div(gi).rem);
            var dfi = new EPoly<Padic>(gi, fa.Derivative.Div(gi).rem);
            var ri = (dgi * fi / dfi).Poly;
            tmp.Add(gi + ri);
        }

        all.Clear();
        all = tmp.ToList();
    }
}
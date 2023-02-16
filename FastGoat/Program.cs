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

Dictionary<int, Dictionary<int, Polynomial<ZnInt, Xi>>> allCnPolys = new();

Polynomial<ZnInt, Xi> ConwayPoly(int p, int n)
{
    var x = Ring.Polynomial("x", ZnInt.KZero(p));
    var xPow = new Dictionary<int, Polynomial<ZnInt, Xi>>();
    var xp = x;
    for (int k = 1; k <= n; k++)
    {
        var xp0 = xp.One;
        for (int j = 0; j < p; j++)
            xp0 *= xp;

        xp = xp0;
        xPow[k] = xp - x;
    }

    Polynomial<ZnInt, Xi> Poly_PpowD(int d) => xPow[d]; // x.Pow(p.Pow(d)) - x;
    var seq = EnumerableExt.MultiLoop(n.Range().OrderDescending()
        .Select(i => p.Range().Select(j => j == 0 ? x.Zero : x.Pow(i) * ((n - i) % 2 == 0 ? j : p - j))));

    var Px = Poly_PpowD(n);
    var lPolys = IntExt.Dividors(n).Select(m => (m, l: Poly_PpowD(m))).ToList();

    var cnPoly = x.Zero;
    var mnm = new Monom<Xi>(x.Indeterminates);
    var t = 0;
    foreach (var lt in seq)
    {
        ++t;
        var lx = x.Pow(n) + lt.Aggregate(x.Zero, (sum, xi) => xi + sum);
        if (lx.IsZero() || lx.Equals(x) || lx[mnm].IsZero() || !Px.Div(lx).rem.IsZero())
            continue;

        var a = new EPolynomial<ZnInt>(x, new PolynomialBasis<ZnInt,Xi>(x.Indeterminates, lx));
        var pn = p.Pow(n) - 1;
        var aPow = new Dictionary<int, EPolynomial<ZnInt>>() { [0] = a.One };
        var acc = a.One;
        var i = 0;
        do
        {
            ++i;
            acc *= a;
            aPow[i] = acc;
        } while (!acc.Equals(a.One));

        if (i != pn)
            continue;

        if (lPolys.All(pi =>
                !pi.l.Div(lx).rem.IsZero() &&
                (n == 1 || allCnPolys[p][pi.m].Substitute(aPow[pn / (p.Pow(pi.m) - 1)],x.ExtractIndeterminate).IsZero())))
        {
            cnPoly = lx;
            break;
        }
    }

    if (!allCnPolys.ContainsKey(p))
        allCnPolys[p] = new() { [1] = cnPoly };
    else
        allCnPolys[p][n] = cnPoly;

    var cnPoly0 = FG.FqX(p.Pow(n)).F.ToPolynomial(x);
    var check = cnPoly.Equals(cnPoly0);

    Console.WriteLine(new { p, n, cnPoly, cnPoly0, t, check });

    return cnPoly;
}

{
    var nb = 4000;
    foreach (var p in Primes10000.Where(p => p * p <= nb))
    {
        var mx = (int)(Double.Log(nb) / Double.Log(p));
        for (int n = 1; n <= mx; n++)
        {
            var poly = ConwayPoly(p, n);
        }
    }
}
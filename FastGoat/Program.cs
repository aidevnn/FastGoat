using System.ComponentModel.DataAnnotations;
using System.Diagnostics;
using System.Numerics;
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
Ring.DisplayPolynomial = MonomDisplay.StarCaret;
ZnInt.Display = ZnDisplay.Signed;

Polynomial<K, Xi> ApplyF<K>(Polynomial<K, Xi> f, KMatrix<K> A, params Polynomial<K, Xi>[] xi)
    where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
{
    if (xi.Length != A.M || xi.Distinct().Count() != xi.Length)
        throw new();

    var nb = xi.Length;
    var monoms = xi.Select(x => x.ExtractIndeterminate).ToArray();

    var kone = f.One;
    var A0 = A.Select(a => a * kone).ToKMatrix(A.M);
    var X = monoms.Select(m => m.ToPolynomial(f.One)).ToKMatrix(A.M);
    var A0X = A0 * X;
    var dicoSubs = nb.Range().ToDictionary(i => monoms[i], i => A0X[i, 0]);
    return f.Substitute(dicoSubs);
}

Polynomial<K, Xi>[] Reynolds<K>(KMatrix<K>[] G, params Polynomial<K, Xi>[] xi)
    where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
{
    var (M, N) = G[0].Dim;
    if (M != N || G.Select(g => g.Dim).Distinct().Count() != 1)
        throw new();

    if (xi.Length != M || xi.Distinct().Count() != xi.Length)
        throw new();

    var f0 = xi[0].One;
    var coefs = G.Length.Range(1).Select(i => xi.Aggregate((a0, a1) => a0 + a1).Pow(i)).Aggregate((a0, a1) => a0 + a1);

    var mn = coefs.Coefs.Keys.Select(m => m.ToPolynomial(f0)).Order().ToArray();
    var facts = mn.Select(f => (f, G.Aggregate(f0.Zero, (acc, g) => acc + ApplyF(f, g, xi)))).ToArray();
    facts.Println(
        $"Reynolds Table {xi.Select((xk, k) => $"{xk}^i{k}").Glue(" * ")} when {xi.Select((xk, k) => $"i{k}").Glue(" + ")}<={G.Length}");
    Console.WriteLine(mn.Length);
    var res = facts.Where(e => !e.Item2.IsZero()).Select(e => e.Item2 / e.Item2.LeadingDetails.lc).Distinct().Order().ToArray();
    return res;
}

(Polynomial<K, Xi>[] inv, Polynomial<K, Xi>[] rfs) InvariantGLnK<K>(ConcreteGroup<KMatrix<K>> G, MonomOrder order = MonomOrder.GrLex)
    where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
{
    var A = G.Neutral();
    var (M, N) = A.Dim;
    if (M != N)
        throw new();

    var xi0 = Ring.Polynomial(A.KZero, order, M.Range().Select(i => $"x{i}").ToArray());
    var Rfs0 = Reynolds(G.ToArray(), xi0);
    xi0[0].Indeterminates.Extend(Rfs0.Length.Range().Select(i => new Xi($"u{i}")).ToArray());
    var xi = xi0[0].Indeterminates.Select(u => u.ToPolynomial(xi0[0].One)).ToArray();
    var Rfs1 = Rfs0.Select((f, i) => f - xi[i + M]).ToArray();
    Rfs1.Println("System");
    var red = Ring.ReducedGrobnerBasis(Rfs1);
    red.Println("Reduced System");
    Console.WriteLine();
    var red2 = SimplifyLoop(red);
    red2.Println("Simplifyed generators");
    Console.WriteLine();
    var ui = red2.SelectMany(s => s.ExtractAllIndeterminates).Distinct().Where(e => e.xi.Contains('u')).Order().ToArray();
    var idl = red2.Where(f => f.ExtractAllIndeterminates.All(u => ui.Contains(u))).ToArray();
    if (idl.Length != 0)
    {
        var uis = idl.SelectMany(f => f.ExtractAllIndeterminates).Distinct().Order().ToArray();
        var inv = Rfs1.Where(f => f.ExtractAllIndeterminates.Distinct().Any(u => uis.Contains(u))).Concat(idl).ToArray();
        inv.Println("Invariant generators");
        Console.WriteLine();
    }

    return (red2, Rfs1);
}

Polynomial<K, Xi>[] Simplify<K>(Polynomial<K, Xi>[] sys) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
{
    var p0 = sys[0];
    var ind = sys[0].Indeterminates;
    var ui = sys.SelectMany(s => s.ExtractAllIndeterminates).Distinct().Where(e => e.xi.Contains('u')).Order().ToArray();

    var xi = ind.Except(ui).ToArray();
    var subs = sys.Grid2D(ui).Where(e =>
            e.t1.Coefs.Keys.Count(e1 => e1[e.t2] == 1 && e1.Degree == 1) == 1 &&
            xi.All(x => !e.t1.ExtractAllIndeterminates.Contains(x)))
        .Select(e => (e.t2, -e.t1 / e.t1.Coefs[new Monom<Xi>(ind, e.t2, 1)] + e.t2.ToPolynomial(p0.One)))
        .GroupBy(e => e.t2).Select(e => (e.Key, e.First().Item2)).ToArray();

    if (subs.Length == 0)
        return sys;

    var sub = subs.MaxBy(e => e.Key);
    Console.WriteLine($"Sub : {sub}");
    var sys2 = sys.Select(s => s.Substitute(sub.Item2, sub.Key)).Where(e => !e.IsZero()).Distinct().ToArray();
    var sys3 = Ring.ReducedGrobnerBasis(sys2);
    return sys3;
}

Polynomial<K, Xi>[] SimplifyLoop<K>(Polynomial<K, Xi>[] sys) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
{
    var sys0 = sys.ToArray();
    var sz = 0;
    while (sys0.Length != sz)
    {
        sz = sys0.Length;
        sys0 = Simplify(sys0);
    }

    return sys0;
}

{
     var b = FG.EPoly(Rational.KZero(), 'b', -3, 0, 1);
     var r1 = -b.One / 2;
     var gl = FG.GLnK($"Q({b})", 2, b);
     var A = gl[r1, b / 2, -b / 2, r1];
     var G = Group.Generate(gl, A);
     DisplayGroup.HeadElements(G);
     InvariantGLnK(G);
}

{
    var (x, y, u, v, w) = Ring.Polynomial(Rational.KZero(), "x", "y", "u", "v", "w").Deconstruct();
    var p = u.Pow(3) - 9 * v.Pow(2) - w.Pow(2);
    var uvw = new[] { u, v, w };
    var (u0, v0, w0) = uvw.Select(e => e.ExtractIndeterminate).Deconstruct();
    var dicoSubs = new Dictionary<Xi, Polynomial<Rational, Xi>>()
    {
        [u0] = x.Pow(2) + y.Pow(2),
        [v0] = x.Pow(2) * y - y.Pow(3) / 3,
        [w0] = x.Pow(3) - 3 * x * y.Pow(2)
    };

    dicoSubs.Println();

    for (int k = 0; k < 5; k++)
    {
        var r = 4.Range().Select(i => (Rng.Next(0, 5), Rng.Next(0, 5), Rng.Next(0, 5), Rng.Next(-4, 5))).ToArray();
        var Fuvw = r.Select(e => (u.Pow(e.Item1) - v.Pow(e.Item2) + w.Pow(e.Item3)) * e.Item4).Aggregate((a0, a1) => a0 + a1)
            .Div(p).rem;

        var Fxyz = Fuvw.Substitute(dicoSubs);
        Console.WriteLine($"{Fxyz} = 0");
        Console.WriteLine("g=lambda x,y:{0}", Fxyz.ToString().Replace("^", "**"));
        Console.WriteLine();
    }
    /*
     Ideals, Varieties, and Algorithms. chap 7 Invariant Theory of Finite Groups page 345, 385
     
    Invariant generators
    x^2 + y^2 - u
    x^2*y - 1/3*y^3 - v
    x^3 - 3*x*y^2 - w
    u^3 - 9*v^2 - w^2

    k[x,y]C3 = k[] ~ k[u,v,w]/(u^3 - 9*v^2 - w^2)

    Random gen Example
    x^8*y^4 - 4/3*x^6*y^6 + 2/3*x^4*y^8 - 4/27*x^2*y^10 
    + 1/81*y^12 - 3*x^8 - 12*x^6*y^2 - 18*x^4*y^4 - 12*x^2*y^6 
    - 3*y^8 - 2*x^6 + 3*x^4*y^2 - 12*x^2*y^4 - y^6 + 2*x^3 
    - 2*x^2*y - 6*x*y^2 + 2/3*y^3 + 3*x^2 + 3*y^2 + 2 = 0

    python code

    def draw(delta,xm,ym,f):
       xr=np.arange(-xm,xm,delta)
       yr=np.arange(-ym,ym,delta)
       x,y=np.meshgrid(xr,yr)
       plt.contour(f(x,y),[0])
       plt.show()

    g=lambda x,y:x**8*y**4 - 4/3*x**6*y**6 + 2/3*x**4*y**8 - 4/27*x**2*y**10 + 1/81*y**12 - 3*x**8 - 12*x**6*y**2 - 18*x**4*y**4 - 12*x**2*y**6 - 3*y**8 - 2*x**6 + 3*x**4*y**2 - 12*x**2*y**4 - y**6 + 2*x**3 - 2*x**2*y - 6*x*y**2 + 2/3*y**3 + 3*x**2 + 3*y**2 + 2
    draw(0.025,10,10,g)

     */
}


// {
//     var gl = FG.GLnK("Q", 2, Rational.KZero());
//     
//     var A0 = gl[0, -1, 1, 0];
//     var G0 = Group.Generate("C4", gl, A0);
//     DisplayGroup.HeadElements(G0);
//     InvariantGLnK(G0);
//     
//     var A1 = gl[0, -1, 1, -1];
//     var G1 = Group.Generate("C3", gl, A1);
//     DisplayGroup.HeadElements(G1);
//     InvariantGLnK(G1);
//     
//     var A2 = gl[-1, 0, 0, 1];
//     var A3 = gl[1, 0, 0, -1];
//     var G3 = Group.Generate("V4", gl, A2, A3);
//     DisplayGroup.HeadElements(G3);
//     InvariantGLnK(G3);
//     
//     var G4 = Group.Generate("C2", gl, -gl.Neutral());
//     DisplayGroup.HeadElements(G4);
//     InvariantGLnK(G4);
//     
//     var A4 = gl[0, 1, 1, 0];
//     var G5 = Group.Generate("D8", gl, A0, A4);
//     DisplayGroup.HeadElements(G5);
//     InvariantGLnK(G5);
//     
//     var G2 = Group.Generate("C6", gl, -A1);
//     DisplayGroup.HeadElements(G2);
//     InvariantGLnK(G2);
// }
//
// {
//     var gl = FG.GLnK("Q", 3, Rational.KZero());
//
//     var A0 = gl[
//         -1, 0, 0, 
//         0, 1, 0, 
//         0, 0, 1];
//     var A1 = gl[
//         1, 0, 0, 
//         0, -1, 0, 
//         0, 0, 1];
//     var A2 = gl[
//         1, 0, 0, 
//         0, 1, 0, 
//         0, 0, -1];
//     var G = Group.Generate(gl, A0, A1, A2);
//     DisplayGroup.HeadElements(G);
//     InvariantGLnK(G);
// }

// {
//     var a = FG.EPoly(Rational.KZero(), 'a', -2, 0, 1);
//     Console.WriteLine(a);
//     Console.WriteLine(a.F);
//     var gl = FG.GLnK($"Q({a})", 2, a);
//     var A = gl[1 / a, -1 / a, 1 / a, 1 / a];
//     var G = Group.Generate(gl, A);
//     InvariantGLnK(G);
// }

// {
//     var gl = FG.GLnK("Q", 2, Rational.KZero());
//
//     var A0 = gl[0, -1, 1, 0];
//     var A1 = gl[0, 1, 1, 0];
//     // var G0 = Group.Generate("D8", gl, A0, A1);
//     // DisplayGroup.HeadElements(G0);
//     // InvariantGLnK(G0);
//
//     // f (x, y) = x8 + 2x6 y2 − x5 y3 + 2x4 y4 + x3 y5 + 2x2 y6 + y8
//
//     var G1 = Group.Generate("C4", gl, A0);
//     DisplayGroup.HeadElements(G1);
//     var (invC4, rfs) = InvariantGLnK(G1);
//     var p0 = invC4[0];
//     p0.Indeterminates.Extend(new[] { new Xi("f") });
//     var f = invC4[0].Indeterminates.Last().ToPolynomial(p0.One);
//     var (x, y) = invC4[0].Indeterminates.Take(2).Select(e => e.ToPolynomial(p0.One)).Deconstruct();
//     var fxy = x.Pow(8) + 2 * x.Pow(6) * y.Pow(2) - x.Pow(5) * y.Pow(3) + 2 * x.Pow(4) * y.Pow(4) + x.Pow(3) * y.Pow(5) +
//               2 * x.Pow(2) * y.Pow(6) + y.Pow(8);
//     
//     var sys = invC4.Prepend(fxy - f).ToArray();
//     sys.Println();
//     Ring.ReducedGrobnerBasis(sys).Println();
// }
//
// {
//     var gl = FG.GLnK("F3", 2, ZnInt.ZnZero(3));
//     var a = gl[1, 1, 0, 1];
//     var b = gl[0, 1, 2, 0];
//
//     var g1 = Group.Generate("SL(2, F3)", gl, a, b);
//     DisplayGroup.HeadElements(g1); // 24 elements
//     InvariantGLnK(g1);
// }
//
// {
//     var gl = FG.GLnK("Q", 2, Rational.KZero());
//
//     var m0 = gl[0, 1, -1, 0];
//     var m1 = gl[0, 1, 1, 0];
//
//     var g1 = Group.Generate("D8a", gl, m0, m1);
//     DisplayGroup.HeadElements(g1);
//     InvariantGLnK(g1);
// }
//
// {
//
//     var x = FG.QPoly();
//     var i = FG.NumberFieldQ(x.Pow(2) + 1, "i");
//     var gl = FG.GLnK($"Q({i})", 2, i);
//     
//     var m0 = gl[i, 0, 0, -i];
//     var m1 = gl[0, 1, -1, 0];
//     
//     var g1 = Group.Generate("Q8a", gl, m0, m1);
//     DisplayGroup.HeadElements(g1);
//     InvariantGLnK(g1);
// }

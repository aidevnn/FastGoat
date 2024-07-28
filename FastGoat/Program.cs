using System.CodeDom;
using System.Collections;
using System.Diagnostics;
using System.Numerics;
using FastGoat.Commons;
using FastGoat.Examples;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.Naming;
using FastGoat.Structures.Subgroups;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Characters;
using FastGoat.UserGroup.Floats;
using FastGoat.UserGroup.GModuleN;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Matrix;
using FastGoat.UserGroup.Perms;
using FastGoat.UserGroup.Polynoms;
using FastGoat.UserGroup.Words;
using static FastGoat.Commons.IntExt;
using static FastGoat.Commons.EnumerableExt;
using FastGoat.UserGroup.Padic;
using FastGoat.UserGroup.Words.Tools;
using GroupRegX = System.Text.RegularExpressions;

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
    return f.Substitute(dicoSubs.Select(e => (e.Key, e.Value)).ToList());
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
    xi0[0].Indeterminates.ExtendAppend(Rfs0.Length.Range().Select(i => new Xi($"u{i}")).ToArray());
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

void InvariantC3ContourPlot()
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

        var Fxyz = Fuvw.Substitute(dicoSubs.Select(e => (e.Key, e.Value)).ToList());
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
    
    import numpy as np
    import matplotlib.pyplot as plt

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

(GLn<EPoly<Rational>> gl, KMatrix<EPoly<Rational>> A) RealRotation_2Pi_over_k(int k)
{
    var c = Cnf.Nth(k);
    var n0 = Lcm(c.Re.N, c.Im.N);
    var (cos, sin) = (c.Re.ToCnfN(n0).E, c.Im.ToCnfN(n0).E);
    var gl = FG.GLnK($"Q(ξ{c.N})", 2, cos);
    return (gl, gl[cos, sin, -sin, cos]);
}

void InvariantCn(int n)
{
    var (gl, A) = RealRotation_2Pi_over_k(n);
    var G = Group.Generate(gl, A);
    DisplayGroup.HeadElements(G);
    InvariantGLnK(G);
}

{
    InvariantCn(2);
    InvariantCn(3);
    InvariantCn(4);
    InvariantCn(5);
    InvariantCn(6);
    InvariantCn(7);
    InvariantCn(8);
}

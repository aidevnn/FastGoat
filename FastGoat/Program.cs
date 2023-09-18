using System.Diagnostics;
using System.IO.IsolatedStorage;
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

int GLnpOrder(int n, int p) => p.Pow(n * (n - 1) / 2) * n.Range().Select(k => p.Pow(n - k) - 1).Aggregate(1, (acc, pk) => acc * pk);

IEnumerable<(int n, int p, int o)> CandidatesGLnp(int maxOrder)
{
    foreach (var p in Primes10000)
    {
        var oGL2p = GLnpOrder(2, p);
        if (oGL2p > maxOrder)
            break;
        else
            yield return (2, p, oGL2p);
        for (int n = 3;; n++)
        {
            var oGLnp = GLnpOrder(n, p);
            if (oGLnp > maxOrder)
                break;
            else
                yield return (n, p, oGLnp);
        }
    }
}

void GLnp(int n, int p)
{
    if (!Primes10000.Contains(p))
        throw new();

    var G = FG.ElementaryAbelian(p.Pow(n));
    Console.WriteLine($"###################### {G} ######################");
    var og = GLnpOrder(n, p);
    Console.WriteLine($"|GL{n}({p})| = {og}");
    GlobalStopWatch.AddLap();
    var AutG = Group.AutomorphismGroup(G);
    Console.WriteLine($"|Aut({G})| = {AutG.Count()}");
    GlobalStopWatch.Show($"Aut({G})");

    var one = ZnInt.ZnZero(p);
    var rg = n.Range();
    var ind = rg.Grid2D(rg).Select(e => $"x{e.t1}{e.t2}").ToArray();
    var xis = Ring.Polynomial(one, MonomOrder.GrevLex, ind);
    var dicoInds = xis.Select((e, k) => (k, e.ExtractIndeterminate))
        .ToDictionary(a => a.ExtractIndeterminate, a => a.k);
    var X = xis[0].X(xis[0].ExtractIndeterminate);
    var matInd = xis.ToKMatrix(n);

    var GLnp = FG.GLnp(n, p);

    Polynomial<ZnInt, Xi>[] Eq(Ep<ZnInt> x, Ep<ZnInt> y)
    {
        if (x.Ei.Length != y.Ei.Length)
            throw new();

        var xMat = x.Ei.Select(c => c * X.One).ToKMatrix(n);
        var yMat = y.Ei.Select(c => c * X.One).ToKMatrix(n);
        var r0 = matInd * xMat - yMat;

        return r0.ToArray();
    }

    var gens = new List<Mat>();
    GlobalStopWatch.AddLap();
    foreach (var gen in AutG.GetGenerators())
    {
        var og1 = AutG.ElementsOrders[gen];
        var system = gen.AutMap.SelectMany(kv => Eq(kv.Key, kv.Value)).Where(e => e.NbIndeterminates > 0).Distinct().Order().ToArray();
        var sols = Ring.ReducedGrobnerBasis(system);
        var mat = new int[n * n];
        foreach (var sol in sols)
        {
            var xi = sol.ExtractIndeterminate;
            var poly = sol.ToKPoly(sol.ExtractIndeterminate);
            if (poly.Degree != 1)
                throw new();
            
            var s0 = -poly[0] / poly[1];
            // Console.WriteLine($"{sol}; {poly}; {xi} = {s0};");
            mat[dicoInds[xi]] = s0.K;
        }

        var gMat = GLnp.Create(mat);
        var og2 = Group.Cycle(GLnp, gMat).Count;
        if (og1 != og2)
            throw new();
        
        gens.Add(gMat);
    }

    GlobalStopWatch.Show($"Solve Groebner : {gens.Count}");

    var G1 = Group.Generate(GLnp.Name, GLnp, gens.ToArray());
    var o0 = G1.ElementsOrders.Values.Where(k => k != 1).Min();
    var o1 = G1.ElementsOrders.Values.Where(k => k != 1).Max();
    var gen0 = G1.Where(e => G1.ElementsOrders[e] == o0).Min();
    var gen1 = G1.Where(e => G1.ElementsOrders[e] == o1).Min();
    var G2 = Group.Generate(GLnp.Name, GLnp, gen0, gen1);
    Console.WriteLine($"Gen0 Order {o0}");
    Console.WriteLine(gen0);
    Console.WriteLine($"Gen1 Order {o1}");
    Console.WriteLine(gen1);
    DisplayGroup.Head(G);
    DisplayGroup.Head(AutG);
    DisplayGroup.Head(G2);
    
    GlobalStopWatch.AddLap();
    DisplayGroup.AreIsomorphics(G2, AutG);
    GlobalStopWatch.Show("Isomorphism");
    
    Console.WriteLine();
}

{
    var candidates = CandidatesGLnp(2500).OrderBy(e => e.o).ToArray();
    candidates.Println();
    foreach (var (n, p, o) in candidates)
    {
        GLnp(n, p);
    }
    
    GLnp(4, 2);
}

/*###################### C2 x C2 x C2 x C2 ######################
   |GL4(2)| = 20160
   |Aut(C2 x C2 x C2 x C2)| = 20160
   # Aut(C2 x C2 x C2 x C2) Time:8023 ms
   # Solve Groebner : 9 Time:4209 ms
   Gen0 Order 2
   [0, 0, 0, 1]
   [0, 0, 1, 0]
   [0, 1, 0, 0]
   [1, 0, 0, 0]
   Gen1 Order 15
   [0, 0, 0, 1]
   [0, 0, 1, 0]
   [0, 1, 0, 1]
   [1, 1, 0, 0]
   |C2 x C2 x C2 x C2| = 16
   Type        AbelianGroup
   BaseGroup   C2 x C2 x C2 x C2
   
   |Aut[C2 x C2 x C2 x C2]| = 20160
   Type        NonAbelianGroup
   BaseGroup   Aut(C2 x C2 x C2 x C2)
   
   |GL(4,2)| = 20160
   Type        NonAbelianGroup
   BaseGroup   GL(4,2)
   
   GL(4,2) IsIsomorphicTo Aut[C2 x C2 x C2 x C2] : True
   # Isomorphism Time:56487 ms
*/
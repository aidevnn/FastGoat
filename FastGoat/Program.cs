using System.Collections;
using System.ComponentModel;
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
using FastGoat.UserGroup.GModuleN;
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

HashSet<ZNElt<Tn, Tg>> TwoCocycleCondition<Tn, Tg>(ConcreteGroup<Tn> N, ConcreteGroup<Tg> G, Dictionary<Ep<Tg>, ZNElt<Tn, Tg>> map)
    where Tg : struct, IElt<Tg>
    where Tn : struct, IElt<Tn>
{
    var sys = new HashSet<ZNElt<Tn, Tg>>();
    var G0 = map.Keys.SelectMany(e => e.Ei).Distinct().Order().ToArray();
    foreach (var (r, s, t) in G0.Grid3D(G0, G0))
    {
        var eq = r * map[new(s, t)] + map[new(r, G.Op(s, t))] - map[new(r, s)] - map[new(G.Op(r, s), t)];
        if (!eq.IsZero())
            sys.Add(eq);
    }

    sys.Order().Println($"N:{N.ShortName} by G:{G.ShortName} Sys:{sys.Count()}");
    return sys;
}

void Test<Tn, Tg>(ConcreteGroup<Tn> N, ConcreteGroup<Tg> G) where Tg : struct, IElt<Tg> where Tn : struct, IElt<Tn>
{
    var autN = Group.AutomorphismGroup(N);
    var ops = Group.AllHomomorphisms(G, autN);
    var GxG = Product.GpGenerate(G, G);
    var og = G.Count();
    var Nab = new AbelianDirectSum<Tn>(N);
    var p = Gcd(Nab.Decomp.Select(e => e.o).Distinct().ToArray());
    var x = FG.FqX(p, 'a').X;
    Console.WriteLine(new { p, x.F });
    var di = $"{(og * (og - 2))}".Length;
    var dj = $"{Nab.Decomp.Count - 1}".Length;
    var zi = Enumerable.Repeat('0', di).Glue();
    var zj = Enumerable.Repeat('0', dj).Glue();
    var fmt = $"x{{0:{zi}}}{{1:{zj}}}";
    var xis = (og * (og - 2) + 1).Range().Select(i => Nab.Decomp.Count.Range().Select(j => string.Format(fmt, i, j)).ToArray())
        .ToArray();
    var ind = Ring.Indeterminates(xis.SelectMany(e => e).ToArray());
    foreach (var op in ops)
    {
        var gXis = new Queue<string[]>(xis);
        var L = op.ToMapElt(autN);
        var zero = new ZNElt<Tn, Tg>(x, ind, Nab, L);
        var map = new Dictionary<Ep<Tg>, ZNElt<Tn, Tg>>();
        foreach (var ep in GxG.OrderByDescending(ep => ep.Ei.Count(ei => ei.Equals(G.Neutral()))).ThenBy(ep => ep))
        {
            if (ep.Ei.Count(ei => ei.Equals(G.Neutral())) != 0)
                map[ep] = zero;
            else
            {
                var xi = gXis.Dequeue().Select(e => new Xi(e)).ToArray();
                map[ep] = new(x, ind, Nab, L, xi);
            }
        }

        L.map.OrderBy(e => e.Key).Println("L");
        map.OrderKeys(G).Println("Map");
        var sys = TwoCocycleCondition(N, G, map).SelectMany(e => e.Coefs.Values).ToArray();
        sys.Println($"System:{sys.Length}");
        var list = new List<KMatrix<EPoly<ZnInt>>>();
        var dico = xis.SelectMany(e => e).Select((e, i) => (i, new Monom<Xi>(ind, new Xi(e)))).ToDictionary(e => e.Item2, e => e.i);
        var x0 = new Polynomial<EPoly<ZnInt>, Xi>(ind, x.Zero).Zero;
        var mat = KMatrix<EPoly<ZnInt>>.MergeSameRows(sys.Select(eq0 => dico.Select(e => eq0[e.Key]).ToKMatrix(dico.Count)).ToArray());
        var m0 = mat.T.GetCol(0).Zero;
        var mat0 = KMatrix<EPoly<ZnInt>>.MergeSameRows(mat.T.Cols.Append(m0).ToArray());
        Console.WriteLine(mat.T);
        Console.WriteLine();
        var red = Ring.ReducedRowsEchelonForm(mat.T).A0;
        var red0 = KMatrix<EPoly<ZnInt>>.MergeSameCols(red.Rows.Where(e => e.Any(k => !k.IsZero())).ToArray());
        Console.WriteLine(red);
        Console.WriteLine();
        Console.WriteLine(red0);
        Console.WriteLine();
        var indep = KMatrix<EPoly<ZnInt>>.MergeSameRows(red0.Cols.Where(m => m.Count(k => k.Equals(k.One)) == 1).ToArray());
        Console.WriteLine(indep);
        Console.WriteLine($"Nb X:{dico.Count} => Nb Indep:{indep.N}");
        Console.WriteLine($"Max Solutions:{BigInteger.Pow(p, indep.N)}");
        Console.WriteLine();
        return;
    }
}

{
    var (N, G) = (FG.Abelian(2), FG.Abelian(2, 2));
    Test(N, G);
}
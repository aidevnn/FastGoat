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

void RCochain<Tn, Tg>(ConcreteGroup<Tn> N, ConcreteGroup<Tg> G, int r = 1, bool details = false)
    where Tg : struct, IElt<Tg>
    where Tn : struct, IElt<Tn>
{
    var autN = Group.AutomorphismGroup(N);
    var ops = Group.AllHomomorphisms(G, autN);
    foreach (var (op, i) in ops.Select((l, i) => (l, i + 1)))
    {
        var L0 = op.ToMapElt(autN);
        LRCochain(N, G, L0, r, $"Lbl{i}/{ops.Count}", details);
    }
}

Dictionary<Ep<Tg>, ZNElt<Tn, Tg>> Dr<Tn, Tg>(Dictionary<Ep<Tg>, ZNElt<Tn, Tg>> map)
    where Tg : struct, IElt<Tg>
    where Tn : struct, IElt<Tn>
{
    if (map.Count == 0 || map.Keys.Select(e => e.Ei.Length).Distinct().Count() != 1)
        throw new();

    var r = map.Keys.Select(e => e.Ei.Length).First();
    if (r == 0 && map.Count != 1)
        throw new();

    var v = map.First().Value;
    var G = v.L.Domain;
    var Gr_next = Product.GpGenerate($"G{r + 1}", Enumerable.Repeat(G, r + 1).Cast<IGroup<Tg>>().ToArray());

    var mapNext = new Dictionary<Ep<Tg>, ZNElt<Tn, Tg>>();

    Ep<Tg> Chg(Ep<Tg> ep, int j)
    {
        var arr = ep.Ei;
        var arr0 = arr.SkipAt(j).ToArray();
        arr0[j] = G.Op(arr[j], arr[j + 1]);
        return new(arr0);
    }

    foreach (var ep in Gr_next)
    {
        var (e0, ep0) = (ep[0], ep.SkipAt(0));
        // Console.WriteLine($"r:{r} ep:{ep} => {e0} / {ep0}");
        var first = map[ep0].Act(e0);
        var last = map[ep.SkipAt(r)].Act((-1).Pow(r + 1));
        var other = r.Range().Select(j => map[Chg(ep, j)].Act((-1).Pow(j + 1))).Aggregate(v.ZNZero, (acc, z0) => acc + z0);
        mapNext[ep] = (first + other + last).Simplify();
    }

    return mapNext;
}

Dictionary<Ep<Tg>, ZNElt<Tn, Tg>> Cr<Tn, Tg>(ConcreteGroup<Tn> N, ConcreteGroup<Tg> G, MapElt<Tg, Automorphism<Tn>> L, int r)
    where Tg : struct, IElt<Tg>
    where Tn : struct, IElt<Tn>
{
    var Nab = new AbelianDirectSum<Tn>(N);
    var od = Nab.Decomp.Count;
    var Gr = Product.GpGenerate($"G{r}", Enumerable.Repeat(G, r).Cast<IGroup<Tg>>().ToArray());
    var eps = Gr.Where(ep => ep.Ei.All(e => !e.Equals(G.Neutral()))).OrderBy(ep => ep).ToArray();
    var xis = od.Range().Select(i => eps.Length.Range().Select(j => ($"a{i}{j}", Nab.Decomp[i].o)).ToArray()).ToArray();
    var ind = new Indeterminates<Xi>(xis.SelectMany(e => e.Select(ei => new Xi(ei.Item1))).ToArray());
    var xis2 = eps.Length.Range().Select(j => od.Range().Select(i => ind[i * eps.Length + j]).ToArray()).ToArray();
    var z0 = new ZNElt<Tn, Tg>(ind, Nab, L);
    var zi = xis2.Select(xi => new ZNElt<Tn, Tg>(ind, Nab, L, xi)).ToArray();
    var q = new Queue<ZNElt<Tn, Tg>>(zi);
    return Gr.OrderBy(ep => ep.Ei.Count(ei => ei.Equals(G.Neutral()))).ThenBy(ep => ep)
        .ToDictionary(e => e, e => eps.Contains(e) ? q.Dequeue() : z0.ZNZero);
}

void LRCochain<Tn, Tg>(ConcreteGroup<Tn> N, ConcreteGroup<Tg> G, MapElt<Tg, Automorphism<Tn>> L, int r, string lbl = "test",
    bool details = false)
    where Tg : struct, IElt<Tg>
    where Tn : struct, IElt<Tn>
{
    var Cp = Cr(N, G, L, r);
    var Cnext = Dr(Cp);
    // var Cnextnext = Dr(Cnext); // always zero
    
    Cp.OrderKeys(G).Println($"C{r}");
    Cnext.OrderKeys(G).Println($"C{r + 1}");

    var mapC2 = Cnext.ToDictionary(e => e.Key, e => e.Value);
    var prod = 1;
    while (mapC2.Values.Any(z => !z.IsZero()))
    {
        var infos = SolveSystem(mapC2);
        prod *= infos.nb;
        Console.WriteLine($"Variable:{infos.p0} Nb Possibilities:{infos.nb} Total:{prod}");
        mapC2 = infos.map0.ToDictionary(e => e.Key, e => e.Value);
    }

    if (details && r == 1)
    {
        var all = CocyclesDFS.TwoCocycles(N, G, L, lbl);
        all.AllTwoCocycles();
        if (all.AllCoboundaries.Count != prod)
            throw new("###############");
    }

    Console.WriteLine();
}

(Polynomial<ZnInt, Xi> p0, int nb, Dictionary<Ep<Tg>, ZNElt<Tn, Tg>> map0) SolveSystem<Tn, Tg>(
    Dictionary<Ep<Tg>, ZNElt<Tn, Tg>> map)
    where Tg : struct, IElt<Tg>
    where Tn : struct, IElt<Tn>
{
    if (map.Count == 0 || map.Keys.Select(e => e.Ei.Length).Distinct().Count() != 1)
        throw new();

    var zero = map.First().Value.ZNZero;
    var Nab = zero.Nab;
    var gens = Nab.Decomp.ToDictionary(e => e.g, e => e);
    var sys0 = map.Values.SelectMany(v => gens.Select(e => (n: v[e.Key], e.Value.g, e.Value.o))).Where(e => !e.n.IsZero())
        .OrderBy(e => e.n.NbIndeterminates).ThenByDescending(e => e.n).ToArray();
    if (sys0.Length == 0)
        return (zero.Zero, 0, map);

    var sys1 = sys0.GroupBy(e => e.n).ToDictionary(e => e.Key, e => e.Select(e0 => (e0.g, e0.o)).ToArray());
    var sys2 = sys1.Select(eq => (eq, invs: UnInvertible(eq.Value.Max(f => f.o))))
        .Select(e => (e.invs, p: e.eq.Key, go: e.eq.Value, e.eq.Key.ExtractAllIndeterminates
            .Select(xi => (xi, new Monom<Xi>(e.eq.Key.Indeterminates, xi, 1)))
            .Select(m => (m.xi, m.Item2, e.eq.Key[m.Item2]))
            .OrderBy(m => e.invs.ContainsKey(m.Item3.K) ? 0 : 1)
            .ThenBy(m => m.Item3.K)
            .ThenBy(m => m.xi)
            .ToArray()))
        .OrderBy(e => e.Item4.Length)
        .ThenBy(e => e.Item4.Any(f => e.invs.ContainsKey(f.Item3.K)) ? 0 : 1)
        .ThenBy(e => e.Item4.Min(f => f.xi))
        .ThenBy(e => e.p)
        .ToArray();

    var (inv, eq0, go, unknowns) = sys2.First();
    var ord = go.Max(e => e.o);

    var m0 = unknowns[0];
    var gcd = Gcd(eq0.Coefs.Values.Where(e => !e.IsZero()).Select(e => e.K).ToArray());
    var eq1 = new Polynomial<ZnInt, Xi>(eq0.Indeterminates, eq0.KZero,
        new(eq0.Coefs.ToDictionary(e => e.Key, e => e.Value.One * (e.Value.K / gcd))));
    var p0 = eq1.X(m0.xi) * (m0.Item3.K / gcd);
    var p1 = ZNElt<Tn, Tg>.Mod(p0 - eq1, ord);
    if (!p1.IsZero() && !inv.ContainsKey(m0.Item3.K / gcd))
    {
        Console.WriteLine($"Eq0:{eq0} = 0 <=> Eq1:{eq1} = 0");
        map.OrderKeys(zero.L.Domain).Println($"m0:{m0} p1:{p1}");
        throw new("Solve Problem");
    }

    var map0 = map.ToDictionary(e => e.Key, e => e.Value.Substitue(p1, m0.xi).Simplify());

    var zo = Group.Generate(new Zn(ord), new ZnInt(ord, 1));
    var nb = zo.Select(z => m0.Item3.K * z).Distinct().Count();
    Console.WriteLine($"Eq0:{eq0} = 0<= {p0}={p1}");
    return (p0, nb, map0);
}

{
    var (N, G) = (FG.Abelian(4), FG.Abelian(2, 2));
    RCochain(N, G, r: 0, details: true);
    RCochain(N, G, r: 1, details: true);
    // RCochain(N, G, r: 2, details: false);
}
using FastGoat.Commons;
using FastGoat.Examples;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Integers;

namespace FastGoat.UserGroup.GModuleN;

public static class ZNSolver
{
    public static CrMap<Tn, Tg> Dr<Tn, Tg>(CrMap<Tn, Tg> map)
        where Tg : struct, IElt<Tg>
        where Tn : struct, IElt<Tn>
    {
        if (map.Count == 0 || map.Keys.Select(e => e.Ei.Length).Distinct().Count() != 1)
            throw new();

        var r = map.R;
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
            var first = map[ep0].Act(e0);
            var last = map[ep.SkipAt(r)].Act((-1).Pow(r + 1));
            var other = r.Range().Select(j => map[Chg(ep, j)].Act((-1).Pow(j + 1))).Aggregate(v.ZNZero, (acc, z0) => acc + z0);
            mapNext[ep] = (first + other + last).Simplify();
        }

        return new(mapNext);
    }

    public static CrMap<Tn, Tg> Cr<Tn, Tg>(ConcreteGroup<Tn> N, ConcreteGroup<Tg> G,
        MapElt<Tg, Automorphism<Tn>> L, int r, char c = 'a')
        where Tg : struct, IElt<Tg>
        where Tn : struct, IElt<Tn>
    {
        var Nab = new AbelianDirectSum<Tn>(N);
        var od = Nab.DecompElementary.Count;
        var Gr = Product.GpGenerate($"G{r}", Enumerable.Repeat(G, r).Cast<IGroup<Tg>>().ToArray());
        var eps = Gr.Where(ep => ep.Ei.All(e => !e.Equals(G.Neutral()))).OrderBy(ep => ep).ToArray();
        var xis = od.Range().Select(i => eps.Length.Range().Select(j => $"{c}{i}{j}").ToArray()).ToArray();
        var ind = new Indeterminates<Xi>(xis.SelectMany(e => e.Select(ei => new Xi(ei))).ToArray());
        ind.ExtendAppend(xis.Length.Range().Select(i => new Xi($"q{i}")).ToArray());
        var xis2 = eps.Length.Range().Select(j => od.Range().Select(i => ind[i * eps.Length + j]).ToArray()).ToArray();
        var z0 = new ZNElt<Tn, Tg>(ind, Nab, L);
        var z1 = Enumerable.Repeat(z0.ZNZero, Gr.Count() - eps.Length).ToArray();
        var zi = xis2.Select(xi => new ZNElt<Tn, Tg>(ind, Nab, L, xi)).ToArray();
        var map = Gr.OrderBy(ep => ep.Ei.Count(ei => !ei.Equals(G.Neutral()))).ThenBy(ep => ep)
            .Zip(z1.Concat(zi)).ToDictionary(e => e.First, e => e.Second);
        return new(map);
    }

    public static (CrMap<Tn, Tg> cr, CrMap<Tn, Tg> cnext) LRCochain<Tn, Tg>(ConcreteGroup<Tn> N, ConcreteGroup<Tg> G, MapElt<Tg,
        Automorphism<Tn>> L, int r)
        where Tg : struct, IElt<Tg>
        where Tn : struct, IElt<Tn>
    {
        var cr = Cr(N, G, L, r);
        var cnext = Dr(cr);
        return (cr, cnext);
    }

    public static (SysReduction sred, CrMap<Tn, Tg> mapNext) SolveSystem<Tn, Tg>(CrMap<Tn, Tg> map, Queue<Xi> qis)
        where Tg : struct, IElt<Tg>
        where Tn : struct, IElt<Tn>
    {
        if (map.Count == 0 || map.Keys.Select(e => e.Ei.Length).Distinct().Count() != 1)
            throw new();

        var zero = map.First().Value.ZNZero;
        var pZero = zero.Zero;
        var Nab = zero.Nab;
        var gens = Nab.DecompElementary.ToDictionary(e => e.g, e => e);
        var sys0 = map.Values.SelectMany(v => gens.Select(e => (n: v[e.Key], e.Value.g, e.Value.o))).Where(e => !e.n.IsZero())
            .OrderBy(e => e.n.NbIndeterminates).ThenByDescending(e => e.n).ToArray();
        if (sys0.Length == 0)
            return (new(pZero, pZero.X(pZero.Indeterminates.First()), pZero, 0), map);

        var invs = Nab.DecompElementary.Select(e => e.o).Distinct().ToDictionary(e => e, e => IntExt.UnInvertible(e));

        var sys1 = sys0.GroupBy(e => e.n).ToDictionary(e => e.Key, e => e.Select(e0 => (e0.g, e0.o)).ToArray());
        var sys2 = sys1.Select(eq => (eq, invs: invs[eq.Value.Max(f => f.o)]))
            .Select(e => (e.invs, p: e.eq.Key, go: e.eq.Value, e.eq.Key.ExtractAllIndeterminates
                .Select(xi => (xi, new Monom<Xi>(e.eq.Key.Indeterminates, xi, 1)))
                .Select(m => (m.xi, m.Item2, e.eq.Key[m.Item2]))
                .OrderBy(m => e.invs.ContainsKey(m.Item3.K) ? 0 : 1)
                .ThenBy(m => m.Item3.K)
                .ThenBy(m => m.xi)
                .ToArray()))
            .OrderBy(e => e.Item4.Any(f => e.invs.ContainsKey(f.Item3.K)) ? 0 : 1)
            .ThenBy(e => e.Item4.Length)
            .ThenBy(e => e.Item4.Min(f => f.xi))
            .ThenBy(e => e.p)
            .ToArray();

        var (inv, eq0, go, unknowns) = sys2.First();
        var ord = go.Max(e => e.o);
        var m0 = unknowns[0];

        if (qis.Count == 0 || inv.ContainsKey(m0.Item3.K))
        {
            var p0 = eq0.X(m0.xi) * (m0.Item3.K);
            var p1 = (p0 - eq0).Mod(ord);

            var map0 = map.ToDictionary(e => e.Key, e => e.Value.SubstituteMod(p0, p1).Simplify());
            return (new(eq0, p0, p1, ord), new(map0));
        }
        else
        {
            var gcd = IntExt.Gcd(eq0.Coefs.Values.Where(e => !e.IsZero()).Select(e => e.K).ToArray());
            var eq1 = new Polynomial<ZnInt, Xi>(eq0.Indeterminates, eq0.KZero,
                new(eq0.Coefs.ToDictionary(e => e.Key, e => new ZnInt(e.Value.Mod, e.Value.K / gcd))));
            var xq = qis.Dequeue();
            var k0 = eq0[m0.Item2].K / gcd;
            var k1 = ord / gcd;
            var eq2 = eq1 - k1 * eq1.X(xq);
            var p0 = eq2.X(m0.xi) * k0;
            var p1 = (p0 - eq2).Mod(ord);
            // Console.WriteLine(new { eq0, eq1, eq2, p0, p1 });

            var map0 = map.ToDictionary(e => e.Key, e => e.Value.SubstituteMod(p0, p1).Simplify());
            return (new(eq0, p0, p1, ord), new(map0));
        }
    }

    public static (int nbPos, CrMap<Tn, Tg>) Reduce2Coboundaries<Tn, Tg>(ConcreteGroup<Tn> N, ConcreteGroup<Tg> G,
        MapElt<Tg, Automorphism<Tn>> L)
        where Tg : struct, IElt<Tg>
        where Tn : struct, IElt<Tn>
    {
        Console.WriteLine($"2Coboundaries {N.ShortName} and {G.ShortName}");
        var (c1, c2) = LRCochain(N, G, L, r: 1);

        var mapC2 = new CrMap<Tn, Tg>(c2.ToDictionary(e => e.Key, e => e.Value));
        var nbPos = 1;
        var listReds = new List<SysReduction>();
        var qis = new Queue<Xi>();
        while (mapC2.Values.Any(z => !z.IsZero()))
        {
            var (sred, mapNext) = SolveSystem(mapC2, qis);
            listReds.Add(sred);
            var p0 = sred.monomExpr;
            var zn = Group.Generate(new Zn(sred.mod));
            var nb = zn.Select(z => p0.LeadingDetails.lc.K * z).Distinct().Count();
            nbPos *= nb;
            Console.WriteLine($"Variable:{p0} Nb Possibilities:{nb} Total:{nbPos}");
            mapC2 = mapNext;
        }

        Console.WriteLine();
        listReds.Println("Step by step");
        Console.WriteLine($"## Coset Size : {nbPos}");
        if (listReds.Count == 0)
            return (1, new(c2.ToDictionary(e => e.Key, e => e.Value.ZNZero)));

        var zero = listReds.First().monomExpr.Zero;
        var ind = zero.Indeterminates;
        var subs = ind.Except(listReds.Select(s => s.monomExpr.ExtractIndeterminate)).Select(xi => (zero, xi)).ToArray();
        var map0 = c2.ToDictionary(e => e.Key, e => e.Value.Substitute(subs));
        return (nbPos, new CrMap<Tn, Tg>(map0));
    }

    public static void Reduce2Cocycles<Tn, Tg>(ConcreteGroup<Tn> N, ConcreteGroup<Tg> G,
        MapElt<Tg, Automorphism<Tn>> L)
        where Tg : struct, IElt<Tg>
        where Tn : struct, IElt<Tn>
    {
        Console.WriteLine($"2Cocycles {N.ShortName} and {G.ShortName}");
        var (c2, c3) = LRCochain(N, G, L, r: 2);

        var mapC3 = new CrMap<Tn, Tg>(c3.ToDictionary(e => e.Key, e => e.Value));
        var listReds = new List<SysReduction>();
        var qis = new Queue<Xi>(mapC3.Values.First().Indeterminates.Where(xi => xi.xi.Contains('q')));

        while (mapC3.Values.Any(z => !z.IsZero()))
        {
            var (sred, mapNext) = SolveSystem(mapC3, qis);
            listReds.Add(sred);
            mapC3 = mapNext;
        }

        Console.WriteLine();
        listReds.Println("Step by step");
        if (listReds.Count == 0)
            return;

        var subs = listReds.Select(s => (s.monomExpr, s.substitutionExpr)).ToArray();
        var map0 = c2.ToDictionary(e => e.Key, e => e.Value.SubstituteMod(subs));
        map0.OrderKeys(G).Println("2Cocycles");
        var cm = Dr(new CrMap<Tn, Tg>(map0));
        if (!cm.IsZero())
            throw new("@@@@@@@@@@@@@@@@@@@@@@");
    }
}
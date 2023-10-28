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
        ind.ExtendAppend(ind.Length.Range().Select(i => new Xi($"q{i}")).ToArray());
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

    public static (SysReduction sred, CrMap<Tn, Tg> mapNext) SysSolveStep<Tn, Tg>(CrMap<Tn, Tg> map, Queue<Xi> qis)
        where Tg : struct, IElt<Tg>
        where Tn : struct, IElt<Tn>
    {
        if (map.Count == 0 || map.Keys.Select(e => e.Ei.Length).Distinct().Count() != 1)
            throw new("Unexpected situation 1");

        var mapElt0 = map.First();
        var zero = mapElt0.Value.ZNZero;
        var pZero = zero.Zero;

        var sys1 = SysEquations(map);
        if (sys1.Length == 0)
        {
            var a0 = pZero.Indeterminates.First();
            return (new(pZero, pZero.X(a0), 1, a0, pZero, 0), map);
        }

        var (eq, mod, invs, ords, lt, lc) =  sys1.First();
        var xis = lt.ContentIndeterminates.ToArray();
        if (xis.Length != 1)
            throw new("Unexpected situation 2");

        var xq = qis.Dequeue();
        if (invs.TryGetValue(lc.K, out var ki))
        {
            var xi = xis[0];
            var X = eq.X(xi);
            var expr = lc * X;
            var sub = (X - ki * eq).Mod(mod) + mod * eq.X(xq);
            var map0 = map.ToDictionary(e => e.Key, e => e.Value.Substitute(sub, xi).Simplify());
            return (new(eq, expr, ords[lc.K], xi, sub, mod), new(map0));
        }
        else
        {
            var dec = IntExt.PrimesDec(mod);
            var (p, q) = dec.First();
            var pi = (q - 1).Range(1).Select(i => p.Pow(i)).ToArray();
            var pq = pi.Last(pq => eq.Coefs.Values.All(z => z.K % pq == 0));
            var eq1 = new Polynomial<ZnInt, Xi>(eq.Indeterminates, eq.KZero,
                new(eq.Coefs.ToDictionary(e => e.Key, e => new ZnInt(e.Value.Mod, e.Value.K / pq))));
            var (m1, k0) = eq1.Coefs.OrderByDescending(c => c.Key).First(c => invs.ContainsKey(c.Value.K));
            var xi = m1.ContentIndeterminates.First();
            var lc0 = eq[m1];
            var X = eq.X(xi);
            var expr = lc0 * X;
            var k1 = mod / pq;
            var eq2 = eq1 - k1 * eq1.X(xq);
            ki = invs[k0.K];
            var sub = (X - ki * eq2).Mod(mod);
            // Console.WriteLine(new { eq, eq1, eq2, expr, xi, sub });

            var map0 = map.ToDictionary(e => e.Key, e => e.Value.Substitute(sub, xi).Simplify());
            return (new(eq, expr, ords[lc0.K], xi, sub, mod), new(map0));
        }
    }

    public static EquationInfos[] SysEquations<Tn, Tg>(CrMap<Tn, Tg> map)
        where Tg : struct, IElt<Tg>
        where Tn : struct, IElt<Tn>
    {
        if (map.Count == 0 || map.Keys.Select(e => e.Ei.Length).Distinct().Count() != 1)
            throw new("Unexpected situation 1");

        var mapElt0 = map.First();
        var r = mapElt0.Key.Ei.Length;
        var zero = mapElt0.Value.ZNZero;
        var Nab = zero.Nab;
        var sys0 = map.Values.SelectMany(v => Nab.DecompElementary.Select(e => (n: v[e.g], e.o))).Where(e => !e.n.IsZero()).ToArray();
        if (sys0.Length == 0)
            return Array.Empty<EquationInfos>();

        var allOrders = Nab.ElemOrders;
        var invs = Nab.ElemInvertible;

        (Monom<Xi> x, ZnInt z) LT(Polynomial<ZnInt, Xi> eq, int o)
        {
            if (r == 2)
                return eq.Coefs.Select(e => (e.Key, e.Value)).OrderByDescending(e0 => e0.Key).MaxBy(e0 => allOrders[o][e0.Value.K]);

            return eq.Coefs.Select(e => (e.Key, e.Value)).OrderBy(e0 => e0.Key).MaxBy(e0 => allOrders[o][e0.Value.K]);
        }

        return sys0.GroupBy(e => e.n)
            .Select(e => (eq: e.Key, o: e.Max(a0 => a0.o)))
            .Select(e => (e.eq, e.o, mn: LT(e.eq, e.o)))
            .Select(e => new EquationInfos(e.eq, e.o, invs[e.o], allOrders[e.o], e.mn.x, e.mn.z))
            .OrderByDescending(ei => ei.orders[ei.lc.K])
            .ThenBy(ei => ei.equation.ExtractAllIndeterminates.Length)
            .ThenByDescending(e => e.lm)
            .ToArray();
    }

    public static CrMap<Tn, Tg> SysRewrite<Tn, Tg>(CrMap<Tn, Tg> map)
        where Tg : struct, IElt<Tg>
        where Tn : struct, IElt<Tn>
    {
        var map0 = map.ToDictionary(e => e.Key, e => e.Value);
        var mapElt0 = map.First();
        var zero = mapElt0.Value.Zero;
        var xis = zero.Indeterminates.ToArray();
        var subsMns = new Queue<Monom<Xi>>(xis.Select(xi => new Monom<Xi>(zero.Indeterminates, xi)).OrderDescending());
        var eqXis = map0.SelectMany(e => e.Value.Coefs.SelectMany(c => c.Value.ExtractAllIndeterminates)).Distinct().ToArray();
        var eqMns = new Queue<Monom<Xi>>(eqXis.Select(xi => new Monom<Xi>(zero.Indeterminates, xi)).OrderDescending());
        while (eqMns.Count != 0)
        {
            var em = eqMns.Dequeue();
            var sm = subsMns.Peek();
            if (em.Equals(sm))
                subsMns.Dequeue();
            else if (em.CompareTo(sm) == -1)
            {
                subsMns.Dequeue();
                var sxi = sm.ContentIndeterminates.First();
                var sX = zero.X(sxi);
                var exi = em.ContentIndeterminates.First();
                map0 = map0.ToDictionary(e => e.Key, e => e.Value.Substitute(sX, exi).Simplify());
            }
        }

        return new(map0);
    }
    
    public static (int nbPos, CrMap<Tn, Tg>) SysReduction<Tn, Tg>(CrMap<Tn, Tg> c2)
        where Tg : struct, IElt<Tg>
        where Tn : struct, IElt<Tn>
    {
        var map0 = SysRewrite(c2);
        var mapC2 = new CrMap<Tn, Tg>(map0.ToDictionary(e => e.Key, e => e.Value));
        var nbPos = 1;
        var listReds = new List<SysReduction>();
        var qis = new Queue<Xi>(mapC2.Values.First().Indeterminates.Where(xi => xi.xi.Contains('q')));
        while (mapC2.Values.Any(z => !z.IsZero()))
        {
            var (sred, mapNext) = SysSolveStep(mapC2, qis);
            listReds.Add(sred);
            var p0 = sred.expr;
            nbPos *= sred.order;
            Console.WriteLine($"Variable:{p0} Nb Possibilities:{sred.order} Total:{nbPos}");
            mapC2 = mapNext;
        }

        var mapElt0 = map0.First();
        var zero = mapElt0.Value.Zero;
        Console.WriteLine();
        listReds.Println("Step by step");
        Console.WriteLine($"## Solutions Size : {nbPos}");
        var ind = zero.Indeterminates;
        var subs = ind.Except(listReds.Select(s => s.expr.ExtractIndeterminate)).Select(xi => (zero, xi)).ToArray();
        var map1 = map0.ToDictionary(e => e.Key, e => e.Value.Substitute(subs));
        var map2 = SysRewrite(new CrMap<Tn, Tg>(map1));
        return (nbPos, map2);
    }

    public static (int nbPos, CrMap<Tn, Tg>) Reduce2Coboundaries<Tn, Tg>(ConcreteGroup<Tn> N, ConcreteGroup<Tg> G,
        MapElt<Tg, Automorphism<Tn>> L)
        where Tg : struct, IElt<Tg>
        where Tn : struct, IElt<Tn>
    {
        Console.WriteLine($"2Coboundaries {N.ShortName} and {G.ShortName}");
        var (c1, c2) = LRCochain(N, G, L, r: 1);
        var Nab = c2.Values.First().Nab;
        var (nb, map) = SysReduction(c2);
        map.OrderKeys(G).Println($"All 2Coboundaries with Z(N)~Z({Nab.AbElementaries})");
        return (nb, map);
    }

    public static (int nbPos, CrMap<Tn, Tg> map1) Reduce2Cocycles<Tn, Tg>(ConcreteGroup<Tn> N, ConcreteGroup<Tg> G,
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
            var (sred, mapNext) = SysSolveStep(mapC3, qis);
            listReds.Add(sred);
            mapC3 = mapNext;
        }

        var Nab = mapC3.Values.First().Nab;
        Console.WriteLine();
        listReds.Println("Step by step");

        var subs = listReds.Select(s => (s.substitutionExpr, s.xi)).ToArray();
        var map0 = c2.ToDictionary(e => e.Key, e => e.Value.Substitute(subs).Simplify());
        var (nb, map1) = SysReduction(new CrMap<Tn, Tg>(map0));
        map1.OrderKeys(G).Println($"All 2Cocycles with Z(N)~Z({Nab.AbElementaries})");
        var cm = Dr(map1);
        if (!cm.IsZero())
            throw new("@@@@@@@@@@@@@@@@@@@@@@");
        return (nb, map1);
    }
}
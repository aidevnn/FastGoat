using System.Numerics;
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
    public static void DisplayCrMap<Tn, Tg>(string title, params CrMap<Tn, Tg>[] crMaps)
        where Tg : struct, IElt<Tg>
        where Tn : struct, IElt<Tn>
    {
        var c0 = crMaps[0];
        var G = c0.G;
        var table = c0.OrderKeys(G).SelectMany(e => crMaps.Select(c => c[e.Key]).Cast<object>().Prepend(e.Key)).ToArray();
        var mat = Ring.Matrix(c0.Count, table);
        Console.WriteLine(title);
        Ring.DisplayMatrix(mat);
        Console.WriteLine();
    }

    public static void DisplayCrMap<Tn, Tg>(params CrMap<Tn, Tg>[] crMaps)
        where Tg : struct, IElt<Tg>
        where Tn : struct, IElt<Tn>
    {
        DisplayCrMap("", crMaps);
    }

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

    public static CrMap<Tn, Tg> TDr<Tn, Tg>(CrMap<Tn, Tg> map)
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
        foreach (var ep in Gr_next)
        {
            var other = r.Range().Select(j => map[ep.SkipAt(j)].Act((-1).Pow(j))).Aggregate(v.ZNZero, (acc, z0) => acc + z0);
            mapNext[ep] = other.Simplify();
        }

        return new(mapNext);
    }

    public static CrMap<Tn, Tg> Cr<Tn, Tg>(ConcreteGroup<Tn> N, ConcreteGroup<Tg> G,
        MapElt<Tg, Automorphism<Tn>> L, int r)
        where Tg : struct, IElt<Tg>
        where Tn : struct, IElt<Tn>
    {
        var Nab = new AbelianDirectSum<Tn>(N);
        var od = Nab.DecompElementary.Count;
        var Gr = Product.GpGenerate($"G{r}", Enumerable.Repeat(G, r).Cast<IGroup<Tg>>().ToArray());
        var eps = Gr.Where(ep => ep.Ei.All(e => !e.Equals(G.Neutral()))).OrderBy(ep => ep).ToArray();
        var xis = (od * eps.Length).Range().Select(i => $"a{i}").ToArray();
        var ind = new Indeterminates<Xi>(xis.Select(ei => new Xi(ei)).ToArray());
        var nbInds = ind.Length;
        ind.ExtendAppend(nbInds.Range().Select(i => new Xi($"q{i}")).ToArray());
        ind.SetOrder(MonomOrder.GrLex);
        var z0 = new ZNElt<Tn, Tg>(ind, Nab, L);
        var z1 = Enumerable.Repeat(z0.ZNZero, Gr.Count() - eps.Length).ToArray();
        var zi = ind.Chunk(od).Select(xi => new ZNElt<Tn, Tg>(ind, Nab, L, xi)).ToArray();
        var map = Gr.OrderBy(ep => ep.Ei.Count(ei => !ei.Equals(G.Neutral()))).ThenBy(ep => ep)
            .Zip(z1.Concat(zi)).ToDictionary(e => e.First, e => e.Second);

        return new(map);
    }

    public static CrMap<Tn, Tg> TCr<Tn, Tg>(ConcreteGroup<Tn> N, ConcreteGroup<Tg> G,
        MapElt<Tg, Automorphism<Tn>> L, int r)
        where Tg : struct, IElt<Tg>
        where Tn : struct, IElt<Tn>
    {
        var Nab = new AbelianDirectSum<Tn>(N);
        var od = Nab.DecompElementary.Count;
        var Gr = Product.GpGenerate($"G{r}", Enumerable.Repeat(G, r).Cast<IGroup<Tg>>().ToArray());
        var Grp = (Gp<Tg>)Gr.BaseGroup;
        var Gr0 = Gr.Where(ep => ep.Ei[0].Equals(G.Neutral())).ToArray();
        var eps = Gr0.Where(ep => ep.Ei.Skip(1).Any(e => !e.Equals(G.Neutral()))).OrderBy(ep => ep).ToArray();
        var xis = (od * eps.Length).Range().Select(i => $"a{i}").ToArray();
        var ind = new Indeterminates<Xi>(xis.Select(ei => new Xi(ei)).ToArray());
        var nbInds = ind.Length;
        ind.ExtendAppend(nbInds.Range().Select(i => new Xi($"q{i}")).ToArray());
        ind.SetOrder(MonomOrder.GrLex);
        var z0 = new ZNElt<Tn, Tg>(ind, Nab, L);
        var zi = ind.Chunk(od).Select(xi => new ZNElt<Tn, Tg>(ind, Nab, L, xi)).ToArray();
        var map0 = Gr0.OrderBy(ep => ep.Ei.Count(ei => !ei.Equals(G.Neutral()))).ThenBy(ep => ep)
            .Zip(zi.Prepend(z0)).ToDictionary(e => e.First, e => e.Second);
        var map = G.SelectMany(g =>
            map0.ToDictionary(e => Grp.Act(g, e.Key), e => e.Value.Act(g).Simplify())).ToDictionary(e => e.Key, e => e.Value);
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

    public static (CrMap<Tn, Tg> cr, CrMap<Tn, Tg> cnext) LRTCochain<Tn, Tg>(ConcreteGroup<Tn> N, ConcreteGroup<Tg> G, MapElt<Tg,
        Automorphism<Tn>> L, int r)
        where Tg : struct, IElt<Tg>
        where Tn : struct, IElt<Tn>
    {
        var cr = TCr(N, G, L, r + 1);
        var cnext = TDr(cr);
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
        var Nab = zero.Nab;
        var allOrders = Nab.ElemOrders;
        var invs = Nab.ElemInvertible;

        var sys1 = SysEquations(map);
        if (sys1.Length == 0)
        {
            var a0 = pZero.Indeterminates.First();
            return (new(pZero, pZero.X(a0), 1, a0, pZero, 0), map);
        }

        var (eq, mod, lt, lc) = sys1.First();
        var ords = allOrders[mod];
        var xis = lt.ContentIndeterminates.ToArray();
        if (xis.Length != 1)
            throw new("Unexpected situation 2");

        var xq = qis.Dequeue();
        if (invs[mod].TryGetValue(lc.K, out var ki))
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
            var pq = pi.Last(pq => eq.Coefs.All(c => c.Value.K % pq == 0));
            var eq1 = new Polynomial<ZnInt, Xi>(eq.Indeterminates, eq.KZero,
                new(eq.Coefs.ToDictionary(e => e.Key, e => new ZnInt(e.Value.Mod, e.Value.K / pq))));
            var (m1, k0) = eq1.Coefs.OrderByDescending(c => c.Key).First(c => invs[mod].ContainsKey(c.Value.K));
            var xi = m1.ContentIndeterminates.First();
            var lc0 = eq[m1];
            var X = eq.X(xi);
            var expr = lc0 * X;
            var k1 = mod / pq;
            var eq2 = eq1 - k1 * eq1.X(xq);
            ki = invs[mod][k0.K];
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

        var r = map.R;
        var mapElt0 = map.First();
        var zero = mapElt0.Value.ZNZero;
        var Nab = zero.Nab;
        var allOrders = Nab.ElemOrders;
        var invs = Nab.ElemInvertible;
        var sys0 = map.Values.SelectMany(v => Nab.DecompElementary.Select(e => (n: v[e.g], e.o))).Where(e => !e.n.IsZero())
            .ToArray();
        if (sys0.Length == 0)
            return Array.Empty<EquationInfos>();

        (Monom<Xi> x, ZnInt z) LT(Polynomial<ZnInt, Xi> eq, int o)
        {
            if (eq.Degree == 0)
            {
                var (lc, lm, _) = eq.LeadingDetails;
                return (lm, lc);
            }

            if (r == 2)
                return eq.Coefs.Where(e => !e.Key.IsOne).Select(e => (e.Key, e.Value)).OrderByDescending(e0 => e0.Key)
                    .MaxBy(e0 => allOrders[o][e0.Value.K]);

            return eq.Coefs.Where(e => !e.Key.IsOne).Select(e => (e.Key, e.Value)).OrderBy(e0 => e0.Key)
                .MaxBy(e0 => allOrders[o][e0.Value.K]);
        }

        return sys0.GroupBy(e => e.n)
            .Select(e => (eq: e.Key, o: e.Min(a0 => a0.o)))
            .Select(e => (e.eq, e.o, mn: LT(e.eq, e.o)))
            .Select(e => new EquationInfos(e.eq, e.o, e.mn.x, e.mn.z))
            .OrderByDescending(ei => invs[ei.mod].ContainsKey(ei.lca.K) ? ei.mod : (!ei.equation.ConstTerm.IsZero() ? 1 : 0))
            .ThenByDescending(ei => allOrders[ei.mod][ei.lca.K])
            .ThenBy(ei => ei.equation.ExtractAllIndeterminates.Length)
            .ThenByDescending(e => e.lma)
            // .ThenBy(e => e.equation.ConstTerm)
            .ThenBy(e => e.lca.K)
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

    public static (List<SysReduction> sreds, CrMap<Tn, Tg>) SysReduction<Tn, Tg>(CrMap<Tn, Tg> c2, bool details = false)
        where Tg : struct, IElt<Tg>
        where Tn : struct, IElt<Tn>
    {
        var map0 = SysRewrite(c2);
        var mapC2 = new CrMap<Tn, Tg>(map0.ToDictionary(e => e.Key, e => e.Value));
        var listReds = new List<SysReduction>();
        var qis = new Queue<Xi>(mapC2.Values.First().Indeterminates.Where(xi => xi.xi.Contains('q')));
        while (!mapC2.IsZero())
        {
            var (sred, mapNext) = SysSolveStep(mapC2, qis);
            listReds.Add(sred);
            mapC2 = mapNext;
        }

        var nbPos = listReds.Aggregate(1, (acc, sred) => acc * sred.order);
        var mapElt0 = map0.First();
        var zero = mapElt0.Value.Zero;
        var ind = zero.Indeterminates;
        var subs = ind.Except(listReds.Select(s => s.expr.ExtractIndeterminate)).Select(xi => (zero, xi)).ToArray();
        var map1 = map0.ToDictionary(e => e.Key, e => e.Value.Substitute(subs));
        var map2 = SysRewrite(new CrMap<Tn, Tg>(map1));
        if (details)
        {
            Console.WriteLine();
            listReds.Println("Step by step SysReduction");
            Console.WriteLine($"## Solutions Size : {nbPos}");
        }

        return (listReds, map2.Add(mapC2));
    }

    public static (SysReduction, CrMap<Tn, Tg>) SysFullRewriteStep<Tn, Tg>(CrMap<Tn, Tg> map, Queue<Xi> qis)
        where Tg : struct, IElt<Tg>
        where Tn : struct, IElt<Tn>
    {
        var mapElt0 = map.First();
        var zero = mapElt0.Value.ZNZero;
        var pZero = zero.Zero;
        var Nab = zero.Nab;
        var allOrders = Nab.ElemOrders;
        var invs = Nab.ElemInvertible;

        var sys1 = SysEquations(map);
        // sys1.Println("Sys1");
        var a0 = pZero.Indeterminates.Last();
        if (sys1.All(e => e.equation.Degree == 0))
        {
            if (sys1.Length != 0)
            {
                var (eq0, mod0, lm0, lc0) = sys1.First();
                return (new(eq0, eq0.X(a0), 1, a0, pZero, mod0), map);
            }

            return (new(pZero, pZero.X(a0), 1, a0, pZero, 0), map);
        }

        var eqInfos = sys1.First(e => e.equation.Degree == 1);
        var (eq, mod, lm, lc) = eqInfos;
        var ords = allOrders[mod];
        var xis = lm.ContentIndeterminates.ToArray();
        if (xis.Length != 1)
            throw new("Unexpected situation 2");

        var c = eq.ConstTerm;
        var xi = xis[0];
        if (!invs[mod].TryGetValue(lc.K, out var ki))
        {
            if (!c.IsZero())
            {
                var (k, d) = mod.Range().Select(k => (k, IntExt.AmodP(lc.K * k + c.K, mod))).OrderBy(e => e.Item2).First();
                var sub = (pZero + k).Mod(mod);
                var sred = new SysReduction(eq, eq.X(xi), ords[lc.K], xi, sub, mod);
                var map2 = map.ToDictionary(e => e.Key, e => e.Value.Substitute(sub, xi).Simplify());
                return (sred, new(map2));
            }
            else
            {
                var map2 = map.ToDictionary(e => e.Key, e => e.Value.Substitute(pZero, xi).Simplify());
                return (new(eq, eq.X(xi), ords[lc.K], xi, pZero, mod), new(map2));
            }
        }
        else
        {
            if (eq.ExtractAllIndeterminates.Length == 1)
            {
                if (lc.K == 1)
                {
                    var xq = qis.Dequeue();
                    var map2 = map.ToDictionary(e => e.Key, e => e.Value.Substitute(mod * eq.X(xq), xi).Simplify());
                    return (new(eq, eq, ords[lc.K], xi, mod * eq.X(xq), mod), new(map2));
                }
                else
                {
                    var xq = qis.Dequeue();
                    var map2 = map.ToDictionary(e => e.Key, e => e.Value.Substitute(mod * eq.X(xq), xi).Simplify());
                    return (new(eq, eq.X(xi), ords[lc.K], xi, ki * eq.X(xi), mod), new(map2));
                }
            }
            else
            {
                var xq = qis.Dequeue();
                var sub = (eq.X(xi) - ki * eq).Mod(mod) + mod * eq.X(xq);
                var map2 = map.ToDictionary(e => e.Key, e => e.Value.Substitute(sub, xi).Simplify());
                return (new(eq, eq.X(xi), ords[lc.K], xi, sub, mod), new(map2));
            }
        }
    }

    public static CrMap<Tn, Tg> SysFullRewrite<Tn, Tg>(CrMap<Tn, Tg> map, bool details = false)
        where Tg : struct, IElt<Tg>
        where Tn : struct, IElt<Tn>
    {
        var map0 = new CrMap<Tn, Tg>(map.ToDictionary(e => e.Key, e => e.Value));
        var map1 = new CrMap<Tn, Tg>(map.ToDictionary(e => e.Key, e => e.Value.ZNZero));
        var listReds = new List<SysReduction>();
        var qis = new Queue<Xi>(map0.Values.First().Indeterminates.Where(xi => xi.xi.Contains('q')));
        var ltMap = new List<CrMap<Tn, Tg>>() { map0.Clone };
        while (!map0.Equals(map1))
        {
            map1 = map0.Clone;
            (var sred, map0) = SysFullRewriteStep(map0, qis);
            ltMap.Add(map0);

            if (sred.mod != 0)
                listReds.Add(sred);
        }

        var subs = listReds.Where(s => !s.substitutionExpr.IsZero()).Select(s => (s.substitutionExpr, s.xi)).ToArray();
        var map2 = map.ToDictionary(e => e.Key, e => e.Value.Substitute(subs).Simplify());
        var map3 = SysRewrite(new CrMap<Tn, Tg>(map2));
        if (details)
        {
            ltMap.Add(map3);
            listReds.Println("Steps");
            DisplayCrMap(ltMap.ToArray());
        }

        return map3;
    }

    public static SysSolution<Tn, Tg> Reduce2Coboundaries<Tn, Tg>(ConcreteGroup<Tn> N, ConcreteGroup<Tg> G,
        MapElt<Tg, Automorphism<Tn>> L, bool details = false)
        where Tg : struct, IElt<Tg>
        where Tn : struct, IElt<Tn>
    {
        var (c1, c2) = LRCochain(N, G, L, r: 1);
        var Nab = c2.Values.First().Nab;
        var (sreds, map) = SysReduction(c2);
        var nb = sreds.Aggregate(1, (acc, sred) => acc * sred.order);

        if (details)
        {
            Console.WriteLine($"2Coboundaries N:{N.ShortName} and G:{G.ShortName} Total:{nb}");
            map.OrderKeys(G).Println($"All 2Coboundaries N:{N.ShortName} and G:{G.ShortName} with Z(N)~Z({Nab.AbElementaries})");
        }

        return new(nb, sreds, map);
    }

    public static SysSolution<Tn, Tg> Reduce2Cocycles<Tn, Tg>(ConcreteGroup<Tn> N, ConcreteGroup<Tg> G,
        MapElt<Tg, Automorphism<Tn>> L, bool details = false)
        where Tg : struct, IElt<Tg>
        where Tn : struct, IElt<Tn>
    {
        var (c2, c3) = LRCochain(N, G, L, r: 2);

        var mapC3 = new CrMap<Tn, Tg>(c3.ToDictionary(e => e.Key, e => e.Value));
        var listReds = new List<SysReduction>();
        var qis = new Queue<Xi>(mapC3.Values.First().Indeterminates.Where(xi => xi.xi.Contains('q')));

        while (!mapC3.IsZero())
        {
            var (sred, mapNext) = SysSolveStep(mapC3, qis);
            listReds.Add(sred);
            mapC3 = mapNext;
        }

        var Nab = mapC3.Values.First().Nab;
        var subs = listReds.Select(s => (s.substitutionExpr, s.xi)).ToArray();
        var map0 = c2.Substitute(subs);
        var (sreds, map1) = SysReduction(map0, details);

        var nb = sreds.Aggregate(1, (acc, sred) => acc * sred.order);

        if (details)
        {
            Console.WriteLine($"2Cocycles N:{N.ShortName} and G:{G.ShortName} Total:{nb}");
            Console.WriteLine();
            listReds.Println("Step by step Reduce2Cocycles");
            map1.OrderKeys(G).Println($"All 2Cocycles N:{N.ShortName} and G:{G.ShortName} with Z(N)~Z({Nab.AbElementaries})");
        }

        var cm = Dr(map1);
        if (!cm.IsZero())
            throw new("@@@@@@@@@@@@@@@@@@@@@@");

        return new(nb, sreds, map1);
    }

    public static void Reduce2Cohomologies<Tn, Tg>(ConcreteGroup<Tn> N, ConcreteGroup<Tg> G,
        MapElt<Tg, Automorphism<Tn>> L, string lbl = "test", bool details = false)
        where Tg : struct, IElt<Tg>
        where Tn : struct, IElt<Tn>
    {
        Console.WriteLine($"############# {lbl,-12} #############");
        Console.WriteLine($"H2(G, N) with N:{N.ShortName} and G:{G.ShortName}");
        var sols2Cobs = Reduce2Coboundaries(N, G, L);
        var sols2Cocs = Reduce2Cocycles(N, G, L);
        var (nb2Cobs, sred2Cobs, map2Cobs) = sols2Cobs;
        var (nb2Cocs, sred2Cocs, map2Cocs) = sols2Cocs;
        var nb2Cohs = nb2Cocs / nb2Cobs;
        var nab = map2Cobs.Nab;
        // if (nab.DecompElementary.Any(e => !IntExt.Primes10000.Contains(e.o)))
        //     throw new("Elementary decomposition of N must contains only simple pGroups");

        Console.WriteLine($"#### {lbl} |B2|:{nb2Cobs} |Z2|:{nb2Cocs} |H2|:{nb2Cohs}");
        Console.WriteLine();

        if (nb2Cobs == 1 || nb2Cohs == 1 || nb2Cocs == 1)
            return;

        var ind = sols2Cocs.allMaps.Values.First().Indeterminates;
        var zero = map2Cocs.PZero;

        var map2Cobs0 = map2Cobs.Recreate(ind);
        var gens2Cobs = map2Cobs0.Generators();
        var gens2Cocs = map2Cocs.Generators();

        Console.WriteLine($"B2(G,N):{nb2Cobs}");
        if (details)
            DisplayCrMap($"2Coboundaries gens ords:[{gens2Cobs.Select(g => g.ord).Glue(" ")}]",
                gens2Cobs.Select(g => g.map).ToArray());

        Console.WriteLine($"Z2(G,N):{nb2Cocs}");
        if (details)
            DisplayCrMap($"2Cocycles ords:[{gens2Cocs.Select(g => g.ord).Glue(" ")}]", gens2Cocs.Select(g => g.map).ToArray());

        var (sreds, _) = SysReduction(map2Cobs0);
        var nb0 = sreds.Aggregate(1, (acc, s) => acc * s.order);
        var lt = new List<(int ord, CrMap<Tn, Tg> map)>();
        var map02 = map2Cobs0.Clone;

        foreach (var (ord, map) in gens2Cocs.OrderByDescending(e => e.ord))
        {
            var ai = zero.Indeterminates.First(xi => xi.xi.Contains('q'));

            var mapi = map.Mul(zero.X(ai));
            var mapai = map02.Add(mapi);
            var (sr, mr) = SysReduction(mapai);
            var nb1 = sr.Aggregate(1, (acc, s) => acc * s.order);
            if (nb0 < nb1)
            {
                lt.Add((ord, map));
                map02 = mr.Clone;
                nb0 = nb1;
            }
        }

        var ais = zero.Indeterminates.Where(xi => xi.xi.Contains('a')).ToArray();
        var mapf = ais.Zip(lt).Select(e => e.Second.map.Mul(zero.X(e.First))).Aggregate((a0, a1) => a0.Add(a1));
        var (sf, mf) = SysReduction(mapf);
        var genf = mf.Generators().Select(e => e.map).Prepend(mf.Zero).ToArray();
        var nb2Cohs0 = genf.Length;

        if (details)
        {
            DisplayCrMap($"2Cohomologies Reprs Tmp:{nb2Cohs0}", genf);
        }

        if (nb2Cohs == nb2Cohs0)
        {
            Console.WriteLine($"@H2(G,N):{nb2Cohs0} Expected:{nb2Cohs}");
            if (details)
                DisplayCrMap("2Cohomologies Reprs", genf);
        }
        else
        {
            var set = new HashSet<CrMap<Tn, Tg>>() { mf.Zero };
            var listReprs = new List<(int ord, CrMap<Tn, Tg> map, HashSet<CrMap<Tn, Tg>> reprs)>();
            listReprs.Add((1, mf.Zero, new() { mf.Zero }));
            foreach (var (ord, map) in lt)
            {
                var set0 = new HashSet<CrMap<Tn, Tg>>();
                var arr = set.ToArray();
                var set1 = new HashSet<CrMap<Tn, Tg>>(set);
                for (int i = 1; i < ord; i++)
                {
                    var m1 = map.Mul(i);
                    var m2 = SysFullRewrite(map2Cobs0.Add(m1), details).ConsTerm();

                    if (!m2.IsZero())
                    {
                        set0.Add(m2);
                        foreach (var m3 in arr)
                        {
                            var m4 = SysFullRewrite(map2Cobs0.Add(m2).Add(m3)).ConsTerm();
                            set1.Add(m4);
                        }
                    }
                }

                set.Clear();
                set.UnionWith(set1);
                listReprs.Add((ord, map, set0));
            }

            if (listReprs.Sum(e => e.reprs.Count) == nb2Cohs)
            {
                Console.WriteLine($"%H2(G,N):{set.Count} Expected:{nb2Cohs}");
                if (details)
                    DisplayCrMap("2Cohomologies Reprs", listReprs.SelectMany(e => e.reprs).ToArray());
            }
            else if (set.Count == nb2Cohs)
            {
                Console.WriteLine($"*H2(G,N):{set.Count} Expected:{nb2Cohs}");
                if (details)
                    DisplayCrMap("2Cohomologies Reprs", set.ToArray());
            }
            else
            {
                DisplayCrMap($"ERROR genf:{genf.Length} Expected:{nb2Cohs}", genf);
                DisplayCrMap($"ERROR listReprs:{listReprs.Sum(e => e.reprs.Count)} Expected:{nb2Cohs}", listReprs.SelectMany(e => e.reprs).ToArray());
                DisplayCrMap($"ERROR set:{set.Count} Expected:{nb2Cohs}", set.ToArray());

                throw new("########");
            }
        }
    }

    public static (SysSolution<Tn, Tg> sols2Cobs, SysSolution<Tn, Tg> sols2Cocs)
        TwoCohomologyOrder<Tn, Tg>(ConcreteGroup<Tn> N, ConcreteGroup<Tg> G, MapElt<Tg, Automorphism<Tn>> L, string lbl = "test")
        where Tg : struct, IElt<Tg>
        where Tn : struct, IElt<Tn>
    {
        Console.WriteLine($"############# {lbl,-12} #############");
        Console.WriteLine($"H2(G, N) with N:{N.ShortName} and G:{G.ShortName}");
        var sols2Cobs = Reduce2Coboundaries(N, G, L);
        var sols2Cocs = Reduce2Cocycles(N, G, L);
        var (nb2Cobs, sred2Cobs, map2Cobs) = sols2Cobs;
        var (nb2Cocs, sred2Cocs, map2Cocs) = sols2Cocs;
        var nb2Cohs = nb2Cocs / nb2Cobs;
        Console.WriteLine($"#### {lbl} |B2|:{nb2Cobs} |Z2|:{nb2Cocs} |H2|:{nb2Cohs}");
        Console.WriteLine();

        var max = BigInteger.Pow(N.Count(), G.Count() - 1);
        if (max < 7000)
        {
            var all = CocyclesDFS.TwoCocycles(N, G, L, lbl);
            all.AllTwoCocycles();
            if (all.AllCoboundaries.Count != nb2Cobs || all.AllCosets.Count != nb2Cohs)
                throw new("############### Error in order H2(G, N)");
        }

        return (sols2Cobs, sols2Cocs);
    }
}
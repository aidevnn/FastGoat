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
        var G = map.G;
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
            var first = e0 * map[ep0];
            var last = (-1).Pow(r + 1) * map[ep.SkipAt(r)];
            var other = r.Range().Select(j => (-1).Pow(j + 1) * map[Chg(ep, j)]).Aggregate(v.ZNZero, (acc, z0) => acc + z0);
            mapNext[ep] = (first + other + last).Simplify();
        }

        return new(G, mapNext);
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
        var G = map.G;
        var Gr_next = Product.GpGenerate($"G{r + 1}", Enumerable.Repeat(G, r + 1).Cast<IGroup<Tg>>().ToArray());

        var mapNext = new Dictionary<Ep<Tg>, ZNElt<Tn, Tg>>();
        foreach (var ep in Gr_next)
        {
            var other = r.Range().Select(j => map[ep.SkipAt(j)].Act((-1).Pow(j))).Aggregate(v.ZNZero, (acc, z0) => acc + z0);
            mapNext[ep] = other.Simplify();
        }

        return new(G, mapNext);
    }

    public static CrMap<Tn, Tg> Cr<Tn, Tg>(ConcreteGroup<Tn> N, ConcreteGroup<Tg> G, MapElt<Tg, Automorphism<Tn>> L, int r)
        where Tg : struct, IElt<Tg>
        where Tn : struct, IElt<Tn>
    {
        if (r < 0)
            throw new();

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

        return new(G, map);
    }

    public static CrMap<Tn, Tg> TCr<Tn, Tg>(ConcreteGroup<Tn> N, ConcreteGroup<Tg> G, MapElt<Tg, Automorphism<Tn>> L, int r)
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
        var map = G.SelectMany(g => map0.ToDictionary(e => Grp.Act(g, e.Key), e => e.Value.Act(g).Simplify()))
            .ToDictionary(e => e.Key, e => e.Value);
        return new(G, map);
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

    static (Monom<Xi> x, ZnInt z) LT(Polynomial<ZnInt, Xi> eq, Dictionary<int, int> orders)
    {
        if (eq.Degree == 0)
        {
            var (lc, lm, _) = eq.LeadingDetails;
            return (lm, lc);
        }

        return eq.Coefs.Where(e => !e.Key.IsOne).Select(e => (e.Key, e.Value)).OrderBy(e0 => e0.Key)
            .MaxBy(e0 => orders[e0.Value.K]);
    }

    public static (SysReduction sred, CrMap<Tn, Tg> mapNext) SysSolveStep<Tn, Tg>(CrMap<Tn, Tg> cr, Queue<Xi> qis)
        where Tg : struct, IElt<Tg>
        where Tn : struct, IElt<Tn>
    {
        if (cr.Count == 0 || cr.Keys.Select(e => e.Ei.Length).Distinct().Count() != 1)
            throw new("Unexpected situation 1");

        var zero = cr.ZNZero;
        var pZero = zero.Zero;
        var Nab = zero.Nab;
        var G = cr.G;
        var allOrders = Nab.ElemOrders;
        var invs = Nab.ElemInvertible;
        var modSolve = Nab.ElemSolve;
        var r = cr.R;
        var decomp = Nab.DecompElementary;

        var sys0 = cr.OrderKeys(G).SelectMany(v => decomp.Select(e => (n: v.Value[e.g], e.o))).Where(e => e.n.Degree != 0).ToArray();
        var sys1 = sys0.Select(e => (eq: e.n, eql: e.n.NbIndeterminates, e.o, mn: LT(e.n, allOrders[e.o])))
            .Select(e => (e.eq, e.eql, mod: e.o, lm: e.mn.x, lc: e.mn.z))
            .OrderByDescending(ei => invs[ei.mod].ContainsKey(ei.lc.K) ? 1 : 0)
            .ThenByDescending(ei => allOrders[ei.mod][ei.lc.K])
            .ThenByDescending(ei => ei.eq.ConstTerm.IsZero())
            .ThenBy(ei => ei.eql)
            .ThenByDescending(e => e.lm)
            .ThenBy(e => e.lc.K)
            .ToArray();

        if (sys1.Length == 0)
        {
            var a0 = pZero.Indeterminates.First();
            return (new(pZero, 1, a0, pZero, 0), cr);
        }

        var (eq, eql, mod, lm, lc) = sys1.First();
        var ords = allOrders[mod];
        var xis = lm.ContentIndeterminates.ToArray();
        if (xis.Length != 1)
        {
            var a0 = pZero.Indeterminates.First();
            return (new(pZero, 1, a0, pZero, 0), cr);
        }

        var xq = qis.Dequeue();
        if (invs[mod].TryGetValue(lc.K, out var ki))
        {
            var xi = xis[0];
            var X = eq.X(xi);
            var sub = (X - ki * eq).Mod(mod) + mod * eq.X(xq);
            var map0 = cr.Substitute(sub, xi);
            return (new(eq, ords[lc.K], xi, sub, mod), map0);
        }
        else
        {
            var gcd = IntExt.Gcd(eq.Coefs.Where(c => !c.Value.IsZero() && !c.Key.IsOne).Select(e => e.Value.K).ToArray());
            var eq1 = new Polynomial<ZnInt, Xi>(eq.Indeterminates, eq.KZero,
                new(eq.Coefs.Where(c => !c.Value.IsZero() && !c.Key.IsOne)
                    .ToDictionary(e => e.Key, e => new ZnInt(e.Value.Mod, e.Value.K / gcd))));
            var (m1, k0) = eq1.Coefs.OrderBy(c => c.Value.K).ThenByDescending(c => c.Key).First(c => invs[mod].ContainsKey(c.Value.K));
            var xi = m1.ContentIndeterminates.First();
            var lc0 = eq[m1];

            var k = modSolve[mod][(lc0.K, eq.ConstTerm.K)];
            var X = eq.X(xi);
            var eq2 = eq1 - ords[gcd] * eq1.X(xq);
            ki = invs[mod][k0.K];
            var sub = (X - ki * eq2 + k).Mod(mod);

            var map0 = cr.Substitute(sub, xi);
            return (new(eq, ords[lc.K], xi, sub, mod), map0);
        }
    }

    public static CrMap<Tn, Tg> SysRewrite<Tn, Tg>(CrMap<Tn, Tg> cr)
        where Tg : struct, IElt<Tg>
        where Tn : struct, IElt<Tn>
    {
        var map = cr.Clone();
        var zero = map.ZZero;
        var xis = zero.Indeterminates.ToArray();
        var eqXis = map.ExtractAllIndeterminates;
        var eqMns = new Queue<Monom<Xi>>(eqXis.Select(xi => new Monom<Xi>(zero.Indeterminates, xi)).OrderDescending());
        var subsMns = new Queue<Monom<Xi>>(xis.Select(xi => new Monom<Xi>(zero.Indeterminates, xi)).OrderDescending());
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
                map = map.Substitute(sX, exi);
            }
        }

        return new(map);
    }

    public static (List<SysReduction> sreds, CrMap<Tn, Tg>) SysReduction<Tn, Tg>(CrMap<Tn, Tg> c2, bool details = false)
        where Tg : struct, IElt<Tg>
        where Tn : struct, IElt<Tn>
    {
        var mapI = SysRewrite(c2);
        var mapC2 = new CrMap<Tn, Tg>(c2.G, mapI.ToDictionary(e => e.Key, e => e.Value));
        var listReds = new List<SysReduction>();
        var qis = new Queue<Xi>(mapC2.ZZero.Indeterminates.Where(xi => xi.xi.Contains('q')));
        while (!mapC2.IsZero())
        {
            (var sred, mapC2) = SysSolveStep(mapC2, qis);
            listReds.Add(sred);
        }

        var zero = mapI.ZZero;
        var ind = zero.Indeterminates;
        var subs = ind.Except(listReds.Select(s => s.xi)).Select(xi => (zero, xi)).ToArray();
        var mapF = mapI.Substitute(subs);
        if (details)
        {
            Console.WriteLine();
            listReds.Println("Step by step SysReduction");
            var dim = listReds.Aggregate(1, (acc, sred) => acc * sred.order);
            Console.WriteLine($"## Solutions Size : {dim}");
        }

        return (listReds, mapF);
    }

    public static CrMap<Tn, Tg> SysFullReduction<Tn, Tg>(CrMap<Tn, Tg> cr)
        where Tg : struct, IElt<Tg>
        where Tn : struct, IElt<Tn>
    {
        var zero = cr.ZZero;
        var ind = zero.Indeterminates;
        var gens = cr.Generators();
        var lt = new List<(int ord, CrMap<Tn, Tg> map)>();
        foreach (var (ord, g) in gens.OrderByDescending(e => e.ord))
        {
            var tmp = lt.Zip(ind).Aggregate(cr.Zero, (acc, e) => acc.Add(e.First.map.Mul(zero.X(e.Second))));
            var g0 = SysRepresentative(tmp.Add(g));
            var ord0 = g0.Cycle().Length;
            if (ord0 != 1)
                lt.Add((ord0, g0));
        }

        return lt.OrderByDescending(e => e.ord).Zip(ind)
            .Aggregate(cr.Zero, (acc, e) => acc.Add(e.First.map.Mul(zero.X(e.Second))));
    }

    public static CrMap<Tn, Tg> SysRepresentative<Tn, Tg>(CrMap<Tn, Tg> cr, bool details = false)
        where Tg : struct, IElt<Tg>
        where Tn : struct, IElt<Tn>
    {
        var map0 = cr.Clone();
        var listReds = new List<SysReduction>();
        var qis = new Queue<Xi>(map0.Values.First().Indeterminates.Where(xi => xi.xi.Contains('q')));
        var ltMap = new List<CrMap<Tn, Tg>>() { cr.ConsTerm(), cr.Clone() };
        while (!map0.Equals(map0.ConsTerm()))
        {
            (var sred, map0) = SysSolveStep(map0, qis);
            ltMap.Add(map0);
            if (sred.mod != 0)
                listReds.Add(sred);
        }

        var subs = listReds.Select(s => (s.substitutionExpr, s.xi)).ToArray();
        var map2 = cr.Substitute(subs);
        var map3 = SysRewrite(new CrMap<Tn, Tg>(map2));
        if (details)
        {
            ltMap.Add(map3);
            listReds.Println("Steps");
            DisplayCrMap(ltMap.ToArray());
        }

        return map3;
    }

    public static SysSolution<Tn, Tg> ReduceCoboundaries<Tn, Tg>(ConcreteGroup<Tn> N, ConcreteGroup<Tg> G,
        MapElt<Tg, Automorphism<Tn>> L, int r = 1, bool details = false)
        where Tg : struct, IElt<Tg>
        where Tn : struct, IElt<Tn>
    {
        var (cr, crNext) = LRCochain(N, G, L, r);
        return ReduceCoboundaries(crNext, details);
    }

    public static SysSolution<Tn, Tg> ReduceCoboundaries<Tn, Tg>(CrMap<Tn, Tg> cr, bool details = false)
        where Tg : struct, IElt<Tg>
        where Tn : struct, IElt<Tn>
    {
        var Nab = cr.Values.First().Nab;
        var (sreds, map) = SysReduction(SysFullReduction(cr));
        var nb = sreds.Aggregate(1, (acc, sred) => acc * sred.order);

        if (details)
        {
            var (N, G) = (cr.N, cr.G);
            Console.WriteLine($"{cr.R}Coboundaries N:{N.ShortName} and G:{G.ShortName} Total:{nb}");
            map.OrderKeys(G).Println($"All {cr.R}Coboundaries N:{N.ShortName} and G:{G.ShortName} with Z(N)~Z({Nab.AbElementaries})");
        }

        return new(nb, sreds, map);
    }

    public static SysSolution<Tn, Tg> ReduceCocycles<Tn, Tg>(ConcreteGroup<Tn> N, ConcreteGroup<Tg> G,
        MapElt<Tg, Automorphism<Tn>> L, int r = 2, bool details = false)
        where Tg : struct, IElt<Tg>
        where Tn : struct, IElt<Tn>
    {
        var (cr, crNext) = LRCochain(N, G, L, r);
        return ReduceCocycles(cr, crNext, details);
    }

    public static SysSolution<Tn, Tg> ReduceCocycles<Tn, Tg>(CrMap<Tn, Tg> cr, CrMap<Tn, Tg> crNext, bool details = false)
        where Tg : struct, IElt<Tg>
        where Tn : struct, IElt<Tn>
    {
        var mapTmp = new CrMap<Tn, Tg>(crNext.G, crNext.ToDictionary(e => e.Key, e => e.Value));
        var listReds = new List<SysReduction>();
        var qis = new Queue<Xi>(mapTmp.Values.First().Indeterminates.Where(xi => xi.xi.Contains('q')));

        while (!mapTmp.IsZero())
        {
            var (sred, mapNext) = SysSolveStep(mapTmp, qis);
            listReds.Add(sred);
            mapTmp = mapNext;
        }

        var Nab = mapTmp.Values.First().Nab;
        var subs = listReds.Select(s => (s.substitutionExpr, s.xi)).ToArray();
        var map0 = cr.Substitute(subs);
        var (sreds, map1) = SysReduction(SysFullReduction(map0), details);

        var nb = sreds.Aggregate(1, (acc, sred) => acc * sred.order);

        if (details)
        {
            var (N, G) = (cr.N, cr.G);
            Console.WriteLine($"{cr.R}Cocycles N:{N.ShortName} and G:{G.ShortName} Total:{nb}");
            Console.WriteLine();
            listReds.Println($"Step by step Reduce{cr.R}Cocycles");
            map1.OrderKeys(G).Println($"All {cr.R}Cocycles N:{N.ShortName} and G:{G.ShortName} with Z(N)~Z({Nab.AbElementaries})");
        }

        var cm = Dr(map1);
        if (!cm.IsZero())
            throw new("@@@@@@@@@@@@@@@@@@@@@@");

        return new(nb, sreds, map1);
    }

    static (HashSet<CrMap<Tn, Tg>>[] reprs, SysSolution<Tn, Tg> solsCobs, SysSolution<Tn, Tg> solsCocs) AllRepresentatives<Tn, Tg>
        (SysSolution<Tn, Tg> solsCobs, SysSolution<Tn, Tg> solsCocs)
        where Tg : struct, IElt<Tg>
        where Tn : struct, IElt<Tn>
    {
        var (nbCobs, sredCobs, mapCobs) = solsCobs;
        var (nbCocs, sredCocs, mapCocs) = solsCocs;
        var nbCohs = nbCocs / nbCobs;
        var r = mapCocs.R;
        var crZero = mapCocs.Zero;

        Console.WriteLine($"#### |B{r}|:{nbCobs} |Z{r}|:{nbCocs} |H{r}|:{nbCohs}");
        Console.WriteLine();

        var zero = mapCocs.ZZero;
        var ind = zero.Indeterminates;

        if (nbCohs == 1)
            return (Array.Empty<HashSet<CrMap<Tn, Tg>>>(), solsCobs, solsCocs);

        var mapCobs0 = mapCobs.Recreate(ind);
        var gensCocs = mapCocs.Generators();

        var nb0 = nbCobs;
        var lt = new List<(int ord, int nb, CrMap<Tn, Tg> map)>();
        var map02 = mapCobs0.Clone();

        foreach (var (k0, ord, map) in gensCocs.OrderBy(e => e.ord).Select((e, i) => (i, e.ord, e.map)))
        {
            if (nb0 == nbCocs)
                break;

            Console.CursorLeft = 0;
            var ai = ind.Except(map02.ExtractAllIndeterminates).Last();

            var mapi = map.Mul(zero.X(ai));
            var mapai = map02.Add(mapi);
            var (sr, mr) = SysReduction(mapai);
            var nb1 = sr.Aggregate(1, (acc, s) => acc * s.order);

            if (nb0 < nb1)
            {
                lt.Add((ord, nb1 / nb0, map));
                map02 = mr.Clone();
                nb0 = nb1;
            }

            Console.Write($"Step:{k0} Gens:{lt.Count}/{gensCocs.Length} Dim:{nb0}/{nbCocs}");
        }

        Console.CursorLeft = 40;
        Console.WriteLine();

        var listReprs = new List<HashSet<CrMap<Tn, Tg>>>();
        listReprs.Add(new() { crZero });
        var gens = mapCobs0.Generators().Select(e => e.map).ToList();
        foreach (var (ord, nb, map) in lt)
        {
            var set0 = new HashSet<CrMap<Tn, Tg>>() { crZero };
            var mapTmp = gens.Zip(ind).Aggregate(mapCobs0.Zero, (acc, e) => acc.Add(e.First.Mul(zero.X(e.Second))));
            for (int i = 1; i < ord; i++)
            {
                var m1 = map.Mul(i);
                var m2 = SysRepresentative(mapTmp.Add(m1)).ConsTerm();
                set0.Add(m2);
                if (set0.Count == nb)
                    break;
            }

            listReprs.Add(set0);
            gens.Add(map);
        }

        var (sreds, mapCobs1) = SysReduction(mapCobs0);
        var solsCobs0 = new SysSolution<Tn, Tg>(nbCobs, sreds, mapCobs1);
        return (listReprs.ToArray(), solsCobs0, solsCocs);
    }

    public static (SysSolution<Tn, Tg> solsCobs, SysSolution<Tn, Tg> solsCocs)
        RCohomologyOrder<Tn, Tg>(ConcreteGroup<Tn> N, ConcreteGroup<Tg> G, MapElt<Tg, Automorphism<Tn>> L, int r = 2,
            string lbl = "test")
        where Tg : struct, IElt<Tg>
        where Tn : struct, IElt<Tn>
    {
        Console.WriteLine($"############# {lbl,-12} #############");
        Console.WriteLine($"H2(G, N) with N:{N.ShortName} and G:{G.ShortName}");
        var solsCobs = ReduceCoboundaries(N, G, L, r - 1);
        var solsCocs = ReduceCocycles(N, G, L, r);
        var (nbCobs, sredCobs, mapCobs) = solsCobs;
        var (nbCocs, sredCocs, mapCocs) = solsCocs;
        var nb2Cohs = nbCocs / nbCobs;
        Console.WriteLine($"#### {lbl} |B{r}|:{nbCobs} |Z{r}|:{nbCocs} |H{r}|:{nb2Cohs}");
        Console.WriteLine();

        var max = BigInteger.Pow(N.Count(), G.Count() - 1);
        if (max < 7000)
        {
            var all = CocyclesDFS.TwoCocycles(N, G, L, lbl);
            all.AllTwoCocycles();
            if (all.AllCoboundaries.Count != nbCobs || all.AllCosets.Count != nb2Cohs)
                throw new($"############### Error in order H{r}(G, N)");
        }

        return (solsCobs, solsCocs);
    }

    public static (CrMap<Tn, Tg>[] solsCohs, SysSolution<Tn, Tg> solsCobs, SysSolution<Tn, Tg> solsCocs)
        ReduceCohomologies<Tn, Tg>(ConcreteGroup<Tn> N, ConcreteGroup<Tg> G, MapElt<Tg, Automorphism<Tn>> L, int r = 2,
            string lbl = "test", bool details = false)
        where Tg : struct, IElt<Tg>
        where Tn : struct, IElt<Tn>
    {
        Console.WriteLine($"############# {lbl,-12} #############");
        Console.WriteLine($"H{r}(G, N) with N:{N.ShortName} and G:{G.ShortName}");
        var solsCobs = ReduceCoboundaries(N, G, L, r - 1);
        var solsCocs = ReduceCocycles(N, G, L, r);

        return ReduceCohomologies(solsCobs, solsCocs, lbl, details);
    }

    public static (CrMap<Tn, Tg>[] solsCohs, SysSolution<Tn, Tg> solsCobs, SysSolution<Tn, Tg> solsCocs)
        ReduceCohomologies<Tn, Tg>(SysSolution<Tn, Tg> solsCobs, SysSolution<Tn, Tg> solsCocs, string lbl = "test",
            bool details = false)
        where Tg : struct, IElt<Tg>
        where Tn : struct, IElt<Tn>
    {
        var (listReprs, sys0, sys1) = AllRepresentatives(solsCobs, solsCocs);
        var (nbCobs, sredCobs, mapCobs) = sys0;
        var (nbCocs, sredCocs, mapCocs) = sys1;
        var nbCohs = nbCocs / nbCobs;
        var allCosets = new HashSet<CrMap<Tn, Tg>>();
        var r = mapCocs.R;
        var crZero = mapCocs.Zero;

        if (nbCohs == 1)
            return (new[] { crZero }, sys0, sys1);

        Console.WriteLine($"B{r}(G,N):{nbCobs}={sredCobs.Select(s => s.order).Glue("x")}");
        Console.WriteLine($"Z{r}(G,N):{nbCocs}={sredCocs.Select(s => s.order).Glue("x")}");

        var k0 = 1;
        foreach (var comb in listReprs.MultiLoop())
        {
            Console.Write($"Cosets:{k0++}/{nbCohs}");
            Console.CursorLeft = 0;
            var map0 = comb.Aggregate(crZero, (acc, l) => acc.Add(l));
            allCosets.Add(map0);
        }

        Console.CursorLeft = 20;
        Console.WriteLine();
        Console.WriteLine($"H{r}(G,N):{allCosets.Count}={listReprs.Select(e => e.Count).Glue("x")} Expected:{nbCohs}");
        if (details)
            DisplayCrMap($"{r}Cohomologies Reprs", allCosets.ToArray());

        if (allCosets.Count != nbCohs)
            Console.WriteLine("?????????????????????????????????");

        return (allCosets.ToArray(), sys0, sys1);

        // {
        //     Console.WriteLine($"ERROR H{r}(G,N):Expected:{nbCohs}");
        //     Console.WriteLine($"    listReprs:{listReprs.Select(e => e.Count).Glue(" x ")}={nbCohs2}");
        //
        //     if (listReprs.Sum(l => l.Count) < 20)
        //         DisplayCrMap($"listReprs:{nbCohs2} Expected:{nbCohs}", listReprs.SelectMany(e => e).ToArray());
        //
        //     DisplayCrMap($"{r}Coboundaries", mapCobs.Generators().Select(g => g.map).ToArray());
        //     DisplayCrMap($"{r}Cocycles", mapCocs.Generators().Select(g => g.map).ToArray());
        //     throw new("########");
        // }
    }
}
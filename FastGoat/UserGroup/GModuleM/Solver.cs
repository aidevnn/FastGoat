using System.Numerics;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Integers;

namespace FastGoat.UserGroup.GModuleM;

public static class Solver
{
    public static Dictionary<Tg, GZNElt<Tn, Tg>> MapsH1<Tn, Tg>(ConcreteGroup<Tn> N, ConcreteGroup<Tg> G)
        where Tn : struct, IElt<Tn>
        where Tg : struct, IElt<Tg>
    {
        if (N.GroupType != GroupType.AbelianGroup)
            throw new GroupException(GroupExceptionType.OnlyAbelianGroups);

        var og = G.Count();
        var nb = og - 1;
        var digits = Enumerable.Repeat("0", $"{nb}".Length).Glue();
        var fmt = $"x{{0:{digits}}}";
        var xis = nb.Range(1).Select(k => new Xi(string.Format(fmt, k))).ToArray();

        var GZX = new GZNElt<Tn, Tg>(N, G, xis);
        var z0 = GZX.ZnEltZero;
        var epz = new Queue<ZNElt<Tn>>(xis.Select(xi => z0.GetUnknown(xi)));

        var mapElt = new Dictionary<Tg, GZNElt<Tn, Tg>>(og);
        foreach (var g in G)
        {
            if (g.Equals(G.Neutral()))
                mapElt[g] = GZX.Zero;
            else
            {
                var map0 = GZX.Zero.Coefs;
                map0[G.Neutral()] = epz.Dequeue();
                mapElt[g] = new(N, G, map0);
            }
        }

        return mapElt;
    }

    public static Dictionary<Ep<Tg>, GZNElt<Tn, Tg>> MapCoboundaries<Tn, Tg>(ConcreteGroup<Tn> N, ConcreteGroup<Tg> G)
        where Tn : struct, IElt<Tn>
        where Tg : struct, IElt<Tg>
    {
        var c = MapsH1(N, G);
        var GxG = Product.GpGenerate(G, G);

        var og = G.Count();
        var mapElt = new Dictionary<Ep<Tg>, GZNElt<Tn, Tg>>(og * og);
        foreach (var ep in GxG)
        {
            mapElt[ep] = ep[0] * c[ep[1]] - c[G.Op(ep[0], ep[1])] + c[ep[0]];
        }

        mapElt.Println($"Map2Coboundaries N:{N.ShortName} by G:{G.ShortName} ");
        return mapElt;
    }

    public static Dictionary<Ep<Tg>, GZNElt<Tn, Tg>> MapCocycles<Tn, Tg>(ConcreteGroup<Tn> N, ConcreteGroup<Tg> G)
        where Tn : struct, IElt<Tn>
        where Tg : struct, IElt<Tg>
    {
        if (N.GroupType != GroupType.AbelianGroup)
            throw new GroupException(GroupExceptionType.OnlyAbelianGroups);

        var og = G.Count();
        var GxG = Product.GpGenerate(G, G);
        var nb = og * (og - 2) + 1;
        var digits = Enumerable.Repeat("0", $"{nb}".Length).Glue();
        var fmt = $"y{{0:{digits}}}";
        var xis = nb.Range(1).Select(k => new Xi(string.Format(fmt, k))).ToArray();

        var GZX = new GZNElt<Tn, Tg>(N, G, xis);
        var z0 = GZX.ZnEltZero;
        var epz = new Queue<ZNElt<Tn>>(xis.Select(xi => z0.GetUnknown(xi)));

        var mapElt = new Dictionary<Ep<Tg>, GZNElt<Tn, Tg>>(og * og);
        foreach (var ep in GxG)
        {
            if (ep[0].Equals(G.Neutral()) || ep[1].Equals(G.Neutral()))
                mapElt[ep] = GZX.Zero;
            else
            {
                var map0 = GZX.Zero.Coefs;
                map0[G.Neutral()] = epz.Dequeue();
                mapElt[ep] = new(N, G, map0);
            }
        }

        mapElt.Println($"Map2Cocycles N:{N.ShortName} by G:{G.ShortName} ");
        return mapElt;
    }

    public static HashSet<GZNElt<Tn, Tg>> TwoCocycleCondition<Tn, Tg>(Dictionary<Ep<Tg>, GZNElt<Tn, Tg>> map)
        where Tg : struct, IElt<Tg>
        where Tn : struct, IElt<Tn>
    {
        var m0 = map.First().Value.Zero;
        var (N, G) = (m0.N, m0.G);
        var sys = new HashSet<GZNElt<Tn, Tg>>();
        foreach (var (r, s, t) in G.Grid3D(G, G))
        {
            var eq = r * map[new(s, t)] + map[new(r, G.Op(s, t))] - map[new(r, s)] - map[new(G.Op(r, s), t)];
            if (!eq.IsZero())
                sys.Add(eq);
        }

        sys.Order().Println($"N:{N.ShortName} by G:{G.ShortName} Sys:{sys.Count()}");
        return sys;
    }

    public static (Xi, Tn)[][] SolveEq2Cocycles<Tn, Tg>(GZNElt<Tn, Tg> gz, MapElt<Tg, Automorphism<Tn>> L)
        where Tg : struct, IElt<Tg>
        where Tn : struct, IElt<Tn>
    {
        if (gz.IsKnown())
            return Array.Empty<(Xi, Tn)[]>();

        var (N, G) = (gz.N, gz.G);
        var xis = gz.GetUnknowns();

        return N.MultiLoop(xis.Length).Select(l => xis.Zip(l).ToArray()).Where(b => gz.Substitute(b).Act(L).IsZero()).ToArray();
    }

    public static IEnumerable<List<(Xi, Tn)>> SolveEq2Cocycles<Tn, Tg>
        (HashSet<GZNElt<Tn, Tg>> sys, List<(Xi, Tn)> sols, MapElt<Tg, Automorphism<Tn>> L)
        where Tg : struct, IElt<Tg>
        where Tn : struct, IElt<Tn>
    {
        var eq = sys.Order().First();
        var newSols = SolveEq2Cocycles(eq, L);
        foreach (var sol in newSols)
        {
            var sols0 = sols.Concat(sol).ToList();
            var newSys = sys.Select(gz0 => gz0.Substitute(sol).Act(L)).Where(gz0 => !gz0.IsZero()).Distinct().ToHashSet();
            if (newSys.Count == 0)
            {
                yield return sols0;
                continue;
            }
            
            foreach (var sol1 in SolveEq2Cocycles(newSys, sols0, L))
                yield return sol1;
        }
    }

    public static IEnumerable<List<(Xi, Tn)>> SolveEq2Cocycles<Tn, Tg>(ConcreteGroup<Tn> N, ConcreteGroup<Tg> G, MapElt<Tg, Automorphism<Tn>> L)
        where Tg : struct, IElt<Tg>
        where Tn : struct, IElt<Tn>
    {
        var map = MapCocycles(N, G);
        var sys = TwoCocycleCondition(map);
        foreach (var sol in SolveEq2Cocycles(sys, new(), L))
            yield return sol;
    }
}
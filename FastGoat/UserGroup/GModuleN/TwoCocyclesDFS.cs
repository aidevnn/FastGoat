using System.Numerics;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;

namespace FastGoat.UserGroup.GModuleN;

public class TwoCocyclesDFS<Tn, Tg> where Tn : struct, IElt<Tn> where Tg : struct, IElt<Tg>
{
    public ConcreteGroup<Tg> G { get; }
    public ConcreteGroup<Ep<Tg>> GxG { get; }
    public Gp<Tg> Gr { get; }
    public ConcreteGroup<Ep2<Tg, Tg>> GxG2 { get; }
    public ConcreteGroup<Ep3<Tg, Tg, Tg>> GxGxG { get; }
    public ConcreteGroup<Tn> ZentrumN { get; }
    public ConcreteGroup<Tn> N { get; }
    public MapElt<Tg, Automorphism<Tn>> L { get; }
    public Dictionary<MapElt<Ep<Tg>, Tn>, HashSet<MapElt<Ep<Tg>, Tn>>> AllCosets { get; }
    public HashSet<MapElt<Tg, Tn>> AllLambda { get; private set; }
    public HashSet<MapElt<Ep<Tg>, Tn>> AllCoboundaries { get; private set; }
    public Tn[] StartN { get; set; }
    public MapGroupBase<Ep<Tg>, Tn> Group_GxG_N { get; }
    public string Lbl { get; }

    public TwoCocyclesDFS(ConcreteGroup<Tn> N0, ConcreteGroup<Tg> G0, MapElt<Tg, Automorphism<Tn>> L0, string lbl)
    {
        (N, G, L) = (N0, G0, L0);
        GxG = Product.GpGenerate(G, G);
        GxG2 = Product.Generate(G, G);
        GxGxG = Product.Generate(G, G, G);
        Gr = Product.Gp(G, G);
        Group_GxG_N = new MapGroupBase<Ep<Tg>, Tn>(GxG, N);
        Lbl = lbl;

        AllCoboundaries = new();
        AllCosets = new();
        StartN = Group.AllSubGroups(N).Keys.SelectMany(sg => sg.GetGenerators()).Prepend(N.Neutral()).Distinct().ToArray();
        AllLambda = new();

        ZentrumN = Group.Zentrum(N);
    }

    public HashSet<MapElt<Ep<Tg>, Tn>> AllTwoCocycles()
    {
        Console.WriteLine($"MaxCoBoundaries:{N.Count()}^{G.Count() - 1} = {BigInteger.Pow(N.Count(), G.Count() - 1)}");
        AllCoboundaries.Clear();
        AllCosets.Clear();
        AllCoboundaries.Add(Group_GxG_N.Neutral());
        step = 0;
        All2Coboundaries(new() { [G.Neutral()] = N.Neutral() });

        var nGxG_N = Group_GxG_N.Neutral();
        AllCosets[nGxG_N] = AllCoboundaries.ToHashSet();
        var nG = G.Neutral();
        var map = G.SelectMany(e => new[] { new Ep<Tg>(nG, e), new Ep<Tg>(e, nG) }).Distinct()
            .ToDictionary(e => e, _ => N.Neutral());

        step = 0;
        All2Cocycles(map);
        Console.WriteLine();
        Console.WriteLine(
            $"## COMPLETE ## Coset Size:{AllCoboundaries.Count,4} Nb Cosets:{AllCosets.Count,4} Total Sols:{AllCoboundaries.Count * AllCosets.Count}");

        return AllCosets.Keys.ToHashSet();
    }

    private int step = 0;

    private void All2Coboundaries(Dictionary<Tg, Tn> c)
    {
        // Console.WriteLine($"{Lbl,-6} step:{++step,-5} MapLength:{c.Count,4} NbSols:{AllCoboundariesGens.Count,3}/{AllCoboundaries.Count}");

        var ng = G.Neutral();
        var g = G.Where(e => !c.ContainsKey(e)).OrderBy(e => G.ElementsOrders[e]).ThenAscending().FirstOrDefault(ng);

        if (g.Equals(G.Neutral()))
            return;

        foreach (var n in ZentrumN)
        {
            ++step;
            var c1 = c.Append(new(g, n)).ToDictionary(e => e.Key, e => e.Value);
            var c1xc1 = c1.Keys.Grid2D(c1.Keys).ToHashSet();
            if (c1xc1.Any(e => !c1.ContainsKey(G.Op(e.t1, e.t2))))
            {
                All2Coboundaries(c1);
                continue;
            }

            var map = c1xc1.Select(e => (e, N.Op(L[e.t1][c1[e.t2]], N.Op(N.Invert(c1[G.Op(e.t1, e.t2)]), c1[e.t1]))))
                .ToDictionary(e => new Ep<Tg>(e.e.t1, e.e.t2), e => e.Item2);

            if (c1.Keys.Grid3D(c1.Keys, c1.Keys).Any(e => !TwoCocycleCondition(map, e.t1, e.t2, e.t3)))
            {
                continue;
            }

            if (c1.Count < G.Count())
            {
                All2Coboundaries(c1);
            }
            else
            {
                var sol = new MapElt<Ep<Tg>, Tn>(GxG, N, map);
                AllCoboundaries.Add(sol);
            }
        }
    }

    public void All2Coboundaries()
    {
        var arrG = G.Where(g => !g.Equals(G.Neutral())).ToArray();
        AllLambda = ZentrumN.MultiLoop(arrG.Length)
            .Select(l => arrG.Zip(l).Append((G.Neutral(), N.Neutral())).ToDictionary(e => e.Item1, e => e.Item2))
            .Select(m => new MapElt<Tg, Tn>(G, N, m)).ToHashSet();

        var maps = AllLambda.Select(c => GxG.Select(e =>
                    (e, N.Op(L[e.Ei[0]][c[e.Ei[1]]], N.Op(N.Invert(c[G.Op(e.Ei[0], e.Ei[1])]), c[e.Ei[0]]))))
                .ToDictionary(e => e.e, e => e.Item2))
            .ToArray();

        AllCoboundaries = maps.Select(map => new MapElt<Ep<Tg>, Tn>(GxG, N, map)).Distinct().Where(f => ValidMap(f)).ToHashSet();
    }

    private void All2Cocycles(Dictionary<Ep<Tg>, Tn> map)
    {
        Console.Write($"{Lbl,-8} step:{++step,-5} NbSols:{AllCosets.Count,3}/{AllCoboundaries.Count * AllCosets.Count}");
        Console.CursorLeft = 0;

        var lG = G.GetGenerators().Where(e => !e.Equals(G.Neutral())).Distinct().ToArray();
        var startG = lG.Grid2D(G).Select(e => new Ep<Tg>(e.t1, G.Op(e.t1, e.t2)))
            .Where(e => !map.ContainsKey(e) && !e.Ei[1].Equals(G.Neutral()))
            .OrderBy(e => (e[0], G.Op(e[0], e[1]))).FirstOrDefault(GxG.Neutral());

        if (startG.Equals(GxG.Neutral()))
            return;

        var (r, s) = startG.Ei.Deconstruct();
        foreach (var n in ZentrumN)
        {
            var next = TwoCocyclesUpdate(map, new() { (r, s, n) });
            if (next.Count == map.Count)
                continue;

            if (next.Count < GxG.Count())
            {
                All2Cocycles(next);
            }
            else
            {
                var sol = new MapElt<Ep<Tg>, Tn>(GxG, N, next);
                if (AllCosets.Values.All(cos => !cos.Contains(sol)) && ValidMap(sol))
                {
                    var prevMaps = AllCosets.Keys.ToArray();
                    var allSols = AllCoboundaries.Select(e => Group_GxG_N.Op(e, sol)).ToHashSet();
                    var sol0 = allSols.OrderMaps(G).First();
                    if (AllCosets.ContainsKey(sol0))
                        Console.WriteLine($"#??????????????? Zentrum bug fix");


                    AllCosets[sol0] = allSols;
                    Console.Write(
                        $"{Lbl,-8} step:{++step,-5} NbSols:{AllCosets.Count,3}/{AllCoboundaries.Count * AllCosets.Count}");
                    Console.CursorLeft = 0;

                    // Z2 group of 2cocycles structure
                    var newSols = Group.GenerateElements(Group_GxG_N, new[] { sol }).ToArray();
                    foreach (var (e0, e1) in newSols.Grid2D(prevMaps))
                    {
                        var sol1 = Group_GxG_N.Op(e0, e1);
                        if (AllCosets.Values.All(cos => !cos.Contains(sol1)) && ValidMap(sol1))
                        {
                            var allSols1 = AllCoboundaries.Select(e => Group_GxG_N.Op(e, sol1)).ToHashSet();
                            var sol2 = allSols1.OrderMaps(G).First();
                            if (AllCosets.ContainsKey(sol2))
                                Console.WriteLine($"##?????????????? Zentrum bug fix");

                            AllCosets[sol2] = allSols1;

                            Console.Write(
                                $"{Lbl,-8} step:{++step,-5} NbSols:{AllCosets.Count,3}/{AllCoboundaries.Count * AllCosets.Count}");
                            Console.CursorLeft = 0;
                        }
                    }
                }
            }
        }
    }

    public bool TwoCocycleCondition(Dictionary<Ep<Tg>, Tn> map, Tg r, Tg s, Tg t)
    {
        // ω(r, s) ω(rs, t) = L(r)[ω(s, t)] ω(r, st).

        var r_s = map[new(r, s)];
        var rs_t = map[new(G.Op(r, s), t)];
        var s_t = map[new(s, t)];
        var r_st = map[new(r, G.Op(s, t))];
        var lhs = N.Op(r_s, rs_t);
        var rhs = N.Op(L[r][s_t], r_st);
        return lhs.Equals(rhs);
    }

    public bool TwoCocycleCondition(Dictionary<Ep<Tg>, Tn> map, Ep<Tg> e0, Ep<Tg> e1)
    {
        // ω(r, s)ω(rs, t) = L(r)[ω(s, t)] ω(r, st).
        // ω(rs, t) ω(r, st)-1 = ω(r, s)-1 L(r)[ω(s, t)].

        if (e0[1].Equals(e1[0]))
        {
            var (r, s, t) = (e0[0], e0[1], e1[1]);

            var r_s = new Ep<Tg>(r, s);
            var s_t = new Ep<Tg>(s, t);
            var rs_t = new Ep<Tg>(G.Op(r, s), t);
            var r_st = new Ep<Tg>(r, G.Op(s, t));
            if (!map.ContainsKey(r_s) || !map.ContainsKey(rs_t) || !map.ContainsKey(s_t) || !map.ContainsKey(r_st))
                return true;

            var v0 = N.Op(map[r_s], map[rs_t]);
            var v1 = N.Op(L[r][map[s_t]], map[r_st]);
            return v0.Equals(v1);
        }
        else if (G.Op(e0[0], e0[1]).Equals(e1[0]))
        {
            var (r, s, rs, t) = (e0[0], e0[1], e1[0], e1[1]);

            var r_s = new Ep<Tg>(r, s);
            var s_t = new Ep<Tg>(s, t);
            var rs_t = new Ep<Tg>(rs, t);
            var r_st = new Ep<Tg>(r, G.Op(s, t));
            if (!map.ContainsKey(r_s) || !map.ContainsKey(rs_t) || !map.ContainsKey(s_t) || !map.ContainsKey(r_st))
                return true;

            var v0 = N.Op(map[r_s], map[rs_t]);
            var v1 = N.Op(L[r][map[s_t]], map[r_st]);
            return v0.Equals(v1);
        }

        return true;
    }

    public bool ValidMap(MapElt<Ep<Tg>, Tn> f)
    {
        return G.Grid3D(G, G).All(e => TwoCocycleCondition(f.map, e.t1, e.t2, e.t3));
    }

    public Dictionary<Ep<Tg>, Tn> TwoCocyclesUpdate(Dictionary<Ep<Tg>, Tn> prev, List<(Tg r, Tg s, Tn n)> set)
    {
        if (set.Any(e => prev.ContainsKey(new(e.r, e.s))))
        {
            Console.WriteLine("never reached code");
            return prev;
        }

        var next = prev.ToDictionary(e => e.Key, e => e.Value);
        foreach (var (r, s, n) in set)
            next[new(r, s)] = n;

        var pmap1 = next.Select(e => (e.Key, e.Value)).ToArray();
        if (pmap1.Grid2D(pmap1).Any(e => !TwoCocycleCondition(next, e.t1.Item1, e.t2.Item1)))
            return prev;

        foreach (var (g1, g2, g3) in GxGxG)
        {
            var g1g2 = G.Op(g1, g2);
            var g2g3 = G.Op(g2, g3);
            var (g1_g2, g1g2_g3, g2_g3, g1_g2g3) = (new Ep<Tg>(g1, g2), new Ep<Tg>(g1g2, g3), new Ep<Tg>(g2, g3),
                new Ep<Tg>(g1, g2g3));

            var (bg1_g2, bg1g2_g3, bg2_g3, bg1_g2g3) = (next.ContainsKey(g1_g2), next.ContainsKey(g1g2_g3),
                next.ContainsKey(g2_g3),
                next.ContainsKey(g1_g2g3));

            // template for copy-paste
            // ω(r, s)ω(rs, t) = L(r)[ω(s, t)] ω(r, st).
            // if (bg1_g2 && bg1g2_g3 && bg2_g3 && bg1_g2g3)
            // {
            //      var (ng1_g2, ng1g2_g3, ng2_g3, ng1_g2g3) = (next[g1_g2], next[g1g2_g3], next[g2_g3], next[g1_g2g3]);
            //      
            // }

            if (bg1_g2 && bg1g2_g3 && bg2_g3 && bg1_g2g3)
            {
                var (ng1_g2, ng1g2_g3, ng2_g3, ng1_g2g3) = (next[g1_g2], next[g1g2_g3], next[g2_g3], next[g1_g2g3]);
                if (!N.Op(ng1_g2, ng1g2_g3).Equals(N.Op(L[g1][ng2_g3], ng1_g2g3)))
                    return prev;

                continue;
            }

            // ω(r, s)ω(rs, t) = L(r)[ω(s, t)] ω(r, st).
            if (!bg1_g2 && bg1g2_g3 && bg2_g3 && bg1_g2g3)
            {
                var (ng1g2_g3, ng2_g3, ng1_g2g3) = (next[g1g2_g3], next[g2_g3], next[g1_g2g3]);
                var ng1_g2 = N.Op(N.Op(L[g1][ng2_g3], ng1_g2g3), N.Invert(ng1g2_g3));
                return TwoCocyclesUpdate(prev, set.Append((g1, g2, ng1_g2)).ToList());
            }

            // ω(r, s)ω(rs, t) = L(r)[ω(s, t)] ω(r, st).
            if (bg1_g2 && !bg1g2_g3 && bg2_g3 && bg1_g2g3)
            {
                var (ng1_g2, ng2_g3, ng1_g2g3) = (next[g1_g2], next[g2_g3], next[g1_g2g3]);
                var ng1g2_g3 = N.Op(N.Invert(ng1_g2), N.Op(L[g1][ng2_g3], ng1_g2g3));
                return TwoCocyclesUpdate(prev, set.Append((g1g2, g3, ng1g2_g3)).ToList());
            }

            // ω(r, s)ω(rs, t) = L(r)[ω(s, t)] ω(r, st).
            if (bg1_g2 && bg1g2_g3 && !bg2_g3 && bg1_g2g3)
            {
                var (ng1_g2, ng1g2_g3, ng1_g2g3) = (next[g1_g2], next[g1g2_g3], next[g1_g2g3]);
                var ng2_g3 = L[G.Invert(g1)][N.Op(N.Op(ng1_g2, ng1g2_g3), N.Invert(ng1_g2g3))];
                return TwoCocyclesUpdate(prev, set.Append((g2, g3, ng2_g3)).ToList());
            }

            // ω(r, s)ω(rs, t) = L(r)[ω(s, t)] ω(r, st).
            if (bg1_g2 && bg1g2_g3 && bg2_g3 && !bg1_g2g3)
            {
                var (ng1_g2, ng1g2_g3, ng2_g3) = (next[g1_g2], next[g1g2_g3], next[g2_g3]);
                var ng1_g2g3 = N.Op(N.Invert(L[g1][ng2_g3]), N.Op(ng1_g2, ng1g2_g3));
                return TwoCocyclesUpdate(prev, set.Append((g1, g2g3, ng1_g2g3)).ToList());
            }
        }

        return next;
    }
}
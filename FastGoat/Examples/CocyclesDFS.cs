using System.Numerics;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Matrix;

namespace FastGoat.Examples;

public static class CocyclesDFS
{
    public class OneCocyclesDFS<Tn, Tg> where Tn : struct, IElt<Tn> where Tg : struct, IElt<Tg>
    {
        public ConcreteGroup<Tg> G { get; }
        public ConcreteGroup<Tn> N { get; }
        public MapElt<Tg, Automorphism<Tn>> L { get; }

        public OneCocyclesDFS(ConcreteGroup<Tg> G0, ConcreteGroup<Automorphism<Tn>> AutN, Homomorphism<Tg, Automorphism<Tn>> L0)
        {
            (N, G) = (AutN.Neutral().AutGroup.G, G0);
            L = new MapElt<Tg, Automorphism<Tn>>(G, AutN, L0.HomMap.ToDictionary(e => e.Key, e => e.Value));
        }


        public OneCocyclesDFS(ConcreteGroup<Tn> N0, ConcreteGroup<Tg> G0, MapElt<Tg, Automorphism<Tn>> L0)
        {
            (N, G, L) = (N0, G0, L0);
        }

        public (Tg, Tn) OneCocycleExpression((Tg, Tn) e0, (Tg, Tn) e1)
        {
            // ω(rs) = L(r)[ω(s)] ω(r)

            var ((g0, n0), (g1, n1)) = (e0, e1);
            var g2 = G.Op(g0, g1);
            var n2 = N.Op(L[g0][n1], n0);
            return (g2, n2);
        }

        public List<MapElt<Tg, Tn>> All1Cocycles()
        {
            return SearchOneCocycles(new() { [G.Neutral()] = N.Neutral() }).Distinct().ToList();
        }

        private IEnumerable<MapElt<Tg, Tn>> SearchOneCocycles(Dictionary<Tg, Tn> current)
        {
            if (current.Count == G.Count())
            {
                var sol = new MapElt<Tg, Tn>(G, N, current);
                yield return sol;
            }

            foreach (var next in NextChildsOneCocycles(current))
            {
                foreach (var sol in SearchOneCocycles(next))
                    yield return sol;
            }
        }

        private Dictionary<Tg, Tn> OneCocycleUpdate(Dictionary<Tg, Tn> prev, Tg g0, Tn n0)
        {
            if (prev.ContainsKey(g0))
                return prev;

            var next = prev.ToDictionary(e => e.Key, e => e.Value);
            next[g0] = n0;
            var g1g2 = next.Grid2D(next).ToArray();
            foreach (var ((g1, n1), (g2, n2)) in g1g2)
            {
                var (g3, n3) = OneCocycleExpression((g1, n1), (g2, n2));
                if (next.TryGetValue(g3, out Tn n4))
                {
                    if (!n4.Equals(n3))
                        return prev;
                }
                else
                {
                    next[g3] = n3;
                }
            }

            return next;
        }

        private IEnumerable<Dictionary<Tg, Tn>> NextChildsOneCocycles(Dictionary<Tg, Tn> current)
        {
            var rems = G.Except(current.Keys).OrderByDescending(e => G.ElementsOrders[e]).ThenAscending().ToArray();
            var nb = current.Count;
            foreach (var (g0, n0) in rems.Grid2D(N))
            {
                var next = OneCocycleUpdate(current, g0, n0);
                if (next.Count == nb)
                    continue;
                else if (next.Count > nb)
                {
                    yield return next;
                }
                else
                    throw new();
            }
        }
    }

    static void All_1_Cocycles<Tn, Tg>(ConcreteGroup<Tn> N, ConcreteGroup<Tg> G)
        where Tg : struct, IElt<Tg>
        where Tn : struct, IElt<Tn>
    {
        var autN = Group.AutomorphismGroup(N);
        var allOps = Group.AllHomomorphisms(G, autN);
        Console.WriteLine($"N:{N.ShortName} and G:{G.ShortName}");
        var all = new Dictionary<MapElt<Tg, Automorphism<Tn>>, HashSet<MapElt<Tg, Tn>>>();
        foreach (var L in allOps)
        {
            var L0 = new MapElt<Tg, Automorphism<Tn>>(G, autN, new(L.HomMap));
            var homol = new OneCocyclesDFS<Tn, Tg>(N, G, L0);
            var allOneCocycles = homol.All1Cocycles();
            Console.WriteLine($"L:{L0}");
            Console.WriteLine($"Count:{allOneCocycles.Count}");
            all[L0] = allOneCocycles.ToHashSet();
            foreach (var co in allOneCocycles)
                Console.WriteLine($"  1co:{co}");

            Console.WriteLine();
        }

        Console.WriteLine($"Total:{all.Values.Sum(v => v.Count)}");
        Console.WriteLine();
    }

    public static void ExampleAll1Cocycle()
    {
        var (c2, c4, c2c2, c2c2c2) = (new Cn(2), new Cn(4), FG.Abelian(2, 2), FG.Abelian(2, 2, 2));
        All_1_Cocycles(c4, c4);
        // All_1_Cocycles(c2c2, c2c2c2);
    }

    public static MapElt<Tg, Tn> ToMapElt<Tg, Tn>(this IMap<Tg, Tn> imap, ConcreteGroup<Tn> N)
        where Tg : struct, IElt<Tg>
        where Tn : struct, IElt<Tn>
    {
        return new(imap.Domain, N, new(imap.Domain.ToDictionary(e => e, e => imap[e])));
    }

    private static Tn[] ArrImages<Tn, Tg>(MapElt<Ep<Tg>, Tn> mapElt, ConcreteGroup<Tg> G)
        where Tg : struct, IElt<Tg>
        where Tn : struct, IElt<Tn>
    {
        return mapElt.map.OrderKeys(G).Select(m => m.Value).ToArray();
    }

    public static IOrderedEnumerable<MapElt<Ep<Tg>, Tn>> OrderMaps<Tn, Tg>(this IEnumerable<MapElt<Ep<Tg>, Tn>> maps,
        ConcreteGroup<Tg> G)
        where Tg : struct, IElt<Tg>
        where Tn : struct, IElt<Tn>
    {
        return maps.OrderBy(m => m,
            Comparer<MapElt<Ep<Tg>, Tn>>.Create((m0, m1) => ArrImages(m0, G).SequenceCompareTo(ArrImages(m1, G))));
    }

    public static IOrderedEnumerable<KeyValuePair<Ep<Tg>, Tn>> OrderKeys<Tn, Tg>(this IEnumerable<KeyValuePair<Ep<Tg>, Tn>> kvs,
        ConcreteGroup<Tg> G)
        where Tg : struct, IElt<Tg>
        where Tn : struct, IElt<Tn>
    {
        return kvs.OrderByDescending(kv => kv.Key.Ei.Count(e => e.Equals(G.Neutral()))).ThenBy(kv => kv.Key);
    }

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

    public static Dictionary<MapElt<Tg, Automorphism<Tn>>, TwoCocyclesDFS<Tn, Tg>> TwoCocycles<Tn, Tg>(ConcreteGroup<Tn> N,
        ConcreteGroup<Tg> G, bool trivialActionOnly = true)
        where Tg : struct, IElt<Tg>
        where Tn : struct, IElt<Tn>
    {
        var autN = Group.AutomorphismGroup(N);
        var allOps = Group.AllHomomorphisms(G, autN);
        var trivL = new Homomorphism<Tg, Automorphism<Tn>>(G, G.ToDictionary(e => e, _ => autN.Neutral()));
        if (trivialActionOnly)
            allOps = new() { trivL };
        return allOps.Select(L => new MapElt<Tg, Automorphism<Tn>>(G, autN, new(L.HomMap)))
            .Select((L, lbl) => new TwoCocyclesDFS<Tn, Tg>(N, G, L, $"Lbl{lbl + 1}/{allOps.Count}"))
            .ToDictionary(e => e.L, e => e);
    }

    public static TwoCocyclesDFS<Tn, Tg> TwoCocycles<Tn, Tg>(ConcreteGroup<Tn> N, ConcreteGroup<Tg> G, MapElt<Tg, Automorphism<Tn>> L,
        string lbl = "test")
        where Tg : struct, IElt<Tg>
        where Tn : struct, IElt<Tn>
    {
        return new TwoCocyclesDFS<Tn, Tg>(N, G, L, lbl);
    }

    public static Dictionary<MapElt<Tg, Automorphism<Tn>>, HashSet<MapElt<Ep<Tg>, Tn>>> All_2_Cocycles_N_G<Tn, Tg>(ConcreteGroup<Tn> N,
        ConcreteGroup<Tg> G, bool trivialActionOnly = true)
        where Tg : struct, IElt<Tg>
        where Tn : struct, IElt<Tn>
    {
        var autN = Group.AutomorphismGroup(N);
        var allOps = Group.AllHomomorphisms(G, autN);
        Console.WriteLine("Start Search All ext N -> E -> G");
        Console.WriteLine($"    with N:{N.ShortName} and G:{G.ShortName}");
        Ring.MatrixDisplayForm = Ring.MatrixDisplay.SquareBracket;
        var all = new Dictionary<MapElt<Tg, Automorphism<Tn>>, HashSet<MapElt<Ep<Tg>, Tn>>>();
        var lbl = 0;
        var trivL = new Homomorphism<Tg, Automorphism<Tn>>(G, G.ToDictionary(e => e, _ => autN.Neutral()));
        if (trivialActionOnly)
            allOps = new() { trivL };

        foreach (var L in allOps)
        {
            var L0 = new MapElt<Tg, Automorphism<Tn>>(G, autN, new(L.HomMap));
            var homol = new TwoCocyclesDFS<Tn, Tg>(N, G, L0, $"Lbl{++lbl}/{allOps.Count}");
            var all2Cocycles = homol.AllTwoCocycles();
            all[L0] = all2Cocycles.ToHashSet();
        }

        Console.WriteLine($"N:{N.ShortName} and G:{G.ShortName}");
        Console.WriteLine($"Total:{all.Values.Sum(v => v.Count)}");
        Console.WriteLine();

        return all;
    }

    public static void TwoCocyclesExamples()
    {
        var (c2, c4, c2c2, c8, c2c2c2, c4c4, c2c4, d8, q8, c16) = (new Cn(2), new Cn(4), FG.Abelian(2, 2), new Cn(8),
            FG.Abelian(2, 2, 2), FG.Abelian(4, 4), FG.Abelian(4, 2), FG.Dihedral(4), FG.Quaternion(8), new Cn(16));

        // Twisted Actions for trivial action of G by N
        // it takes longuer time than before but it is more accure
        All_2_Cocycles_N_G(c4, c4);
        All_2_Cocycles_N_G(c2c2, c4);
        All_2_Cocycles_N_G(c2c4, c4);

        All_2_Cocycles_N_G(c4, c2c2);
        All_2_Cocycles_N_G(c2c2, c2c2);
        All_2_Cocycles_N_G(c2c4, c2c2);

        All_2_Cocycles_N_G(c8, c2);
        All_2_Cocycles_N_G(d8, c2);
        All_2_Cocycles_N_G(c16, c2);
    }

    static HashSet<ExtensionGroup<Tn, Tg>> AllExt<Tn, Tg>(ConcreteGroup<Tn> N, ConcreteGroup<Tg> G)
        where Tg : struct, IElt<Tg>
        where Tn : struct, IElt<Tn>
    {
        var all = All_2_Cocycles_N_G(N, G, trivialActionOnly: false);
        var allExt = new HashSet<ExtensionGroup<Tn, Tg>>();
        foreach (var (L, tws) in all)
        {
            foreach (var w in tws)
            {
                var ext = Group.ExtensionGroup(N, L, w, G);
                if (!ext.ExtBase.IsGroup)
                    continue;

                allExt.Add(ext);
            }
        }

        return allExt;
    }

    public static IEnumerable<(ConcreteGroup<Ep2<Tn, Ep<ZnInt>>>, (int, int, int))> BuildExtensions<Tn>(HashSet<ConcreteGroup<Tn>> gr)
        where Tn : struct, IElt<Tn>
    {
        foreach (var tuple in BuildExtensions(new HashSet<ConcreteGroup<Ep<ZnInt>>>() { FG.Abelian(2) }, gr))
            yield return tuple;
    }

    public static IEnumerable<(ConcreteGroup<Ep2<Tn, Tg>>, (int, int, int))> BuildExtensions<Tn, Tg>(HashSet<ConcreteGroup<Tg>> Gs,
        HashSet<ConcreteGroup<Tn>> Ns)
        where Tn : struct, IElt<Tn>
        where Tg : struct, IElt<Tg>
    {
        var allNg = Ns.Grid2D(Gs).ToArray();
        foreach (var tuple in BuildExtensions(allNg))
            yield return tuple;
    }

    public static (int, int, int) SubGroupsDetails<T>(Dictionary<ConcreteGroup<T>, List<ConcreteGroup<T>>> sg0)
        where T : struct, IElt<T>
    {
        return (sg0.Values.Sum(s => s.Count), sg0.Count, sg0.Count(s => s.Value.Count == 1));
    }

    public static (int, int, int) SubGroupsDetails<T>(ConcreteGroup<T> g) where T : struct, IElt<T>
    {
        return SubGroupsDetails(Group.AllSubGroups(g));
    }

    public static IEnumerable<(ConcreteGroup<Ep2<Tn, Tg>>, (int, int, int))> BuildExtensions<Tn, Tg>(
        (ConcreteGroup<Tn> N, ConcreteGroup<Tg> G)[] AllNG)
        where Tn : struct, IElt<Tn>
        where Tg : struct, IElt<Tg>
    {
        var allExts = AllNG.SelectMany(e => AllExt(e.N, e.G)).ToList();
        var infosSubGroups =
            new Dictionary<ConcreteGroup<Ep2<Tn, Tg>>,
                Dictionary<ConcreteGroup<Ep2<Tn, Tg>>, List<ConcreteGroup<Ep2<Tn, Tg>>>>>();

        var abTypes = new HashSet<int[]>(new SequenceEquality<int>());
        var listIsos = new HashSet<ConcreteGroup<Ep2<Tn, Tg>>>(new IsomorphEquality<Ep2<Tn, Tg>>());

        foreach (var (ext, i) in allExts.Select((ext, i) => (ext, i + 1)))
        {
            Console.Write($"Progress:{i,5} / {allExts.Count} Found:{listIsos.Count}");
            Console.CursorLeft = 0;

            if (ext.GroupType == GroupType.AbelianGroup)
            {
                try
                {
                    var abType = AbelianInvariantsFactors.Reduce(ext).Order().ToArray();
                    if (!abTypes.Add(abType))
                        continue;
                }
                catch (Exception e)
                {
                    Console.WriteLine("@???????????????????????????????????????????????????????????");
                    // TODO why some twisted actions wont give a group
                }
            }

            if (!Group.IsGroup(ext))
            {
                Console.WriteLine("@@??????????????????????????????????????????????????????????");
                continue;
            }

            // var cg = (ConcreteGroup<Ep2<Tn, Tg>>)ext;
            // var res = Parallel.ForEach(listIsos, (g0, state) =>
            // {
            //     // unexpected behaviour
            //     if (cg.IsIsomorphicTo(g0))
            //         state.Stop();
            // });
            //     
            // if (!res.IsCompleted)
            //     continue;
            // else
            //     listIsos.Add(ext);

            if (!listIsos.Add(ext))
                continue;

            var allSubs = Group.AllSubGroups(ext);
            var info0 = SubGroupsDetails(allSubs);

            yield return (ext, info0);

            infosSubGroups[ext] = allSubs;
        }

        Console.WriteLine();
        Console.WriteLine($"AllExts Found : {infosSubGroups.Count}");
    }

    public static void DisplayInfosGroups<Tg>((ConcreteGroup<Tg>, (int, int, int))[] elts, int countStart = 0, bool naming = true)
        where Tg : struct, IElt<Tg>
    {
        foreach (var (g, k, infos) in elts.Select((e, k) => (e.Item1, countStart + k + 1, e.Item2)))
        {
            var og = g.Count();
            var name = naming ? $"Sm{og}[{k}]" : g.Name;
            Console.WriteLine("##########################################################");
            Console.WriteLine($"#################  {name,10} found   ####################");
            Console.WriteLine("##########################################################");
            g.SetName(name);
            DisplayGroup.HeadOrders(g);
            Console.WriteLine($"AllSubGr:{infos.Item1} AllConjsCl:{infos.Item2} AllNorms:{infos.Item3}");
            Console.WriteLine();
        }
    }

    public static void DisplayInfosGroups<Tg>(IEnumerable<ConcreteGroup<Tg>> gs, bool naming = false) where Tg : struct, IElt<Tg>
    {
        DisplayInfosGroups(
            gs.Select(e => (e, SubGroupsDetails(e)))
                .OrderBy(e => e.Item1.GroupType)
                .ThenByDescending(e => e.Item1.ElementsOrders.Values.Max())
                .ThenBy(e => e.Item2).ToArray()
            , naming: naming);
    }

    public static void ExampleAll16Orders()
    {
        var allOrder4 = BuildExtensions(new HashSet<ConcreteGroup<ZnInt>>() { new Cn(2) })
            .OrderBy(e => e.Item1.GroupType)
            .ThenByDescending(e => e.Item1.ElementsOrders.Values.Max())
            .ThenBy(e => e.Item2).ToArray();

        var allOrder8 = BuildExtensions(allOrder4.Select(e => e.Item1)
                .ToHashSet(new IsomorphEquality<Ep2<ZnInt, Ep<ZnInt>>>()))
            .OrderBy(e => e.Item1.GroupType)
            .ThenByDescending(e => e.Item1.ElementsOrders.Values.Max())
            .ThenBy(e => e.Item2).ToArray();

        var allOrder16 = BuildExtensions(allOrder8.Select(e => e.Item1)
                .ToHashSet(new IsomorphEquality<Ep2<Ep2<ZnInt, Ep<ZnInt>>, Ep<ZnInt>>>()))
            .OrderBy(e => e.Item1.GroupType)
            .ThenByDescending(e => e.Item1.ElementsOrders.Values.Max())
            .ThenBy(e => e.Item2).ToArray();

        Console.Clear();

        DisplayInfosGroups(allOrder4);
        DisplayInfosGroups(allOrder8);
        DisplayInfosGroups(allOrder16);

        Console.Beep();
    }

    public static void ExampleAll32Orders()
    {
        var allAb16 = IntExt.Partitions32[4].Select(l => FG.Abelian(l.Select(k => 2.Pow(k)).ToArray()))
            .ToHashSet(); // 5 groups

        var allSdp16 = Group.AllSemiDirectProd("C8x:C2", FG.Abelian(8), new Cn(2))
            .Concat(Group.AllSemiDirectProd("(C4xC2)x:C2", FG.Abelian(4, 2), new Cn(2)))
            .Concat(Group.AllSemiDirectProd("C4x:C4", FG.Abelian(4), new Cn(4)))
            .ToHashSet(new IsomorphEquality<Ep2<Ep<ZnInt>, ZnInt>>()); // 7 groups

        var c2 = new Cn(2);
        var c2q8 = Product.Generate(new Cn(2), FG.Quaternion(8));
        var q16 = Product.Generate(Group.Generate("C1", c2, c2[0]), FG.Quaternion(16));
        q16.SetName("Q16");

        var l0 = new List<(ConcreteGroup<Ep2<Ep2<ZnInt, Mat>, Ep<ZnInt>>>, (int, int, int))>();
        var l1 = new List<(ConcreteGroup<Ep2<Ep2<Ep<ZnInt>, ZnInt>, Ep<ZnInt>>>, (int, int, int))>();
        var l2 = new List<(ConcreteGroup<Ep2<Ep<ZnInt>, Ep<ZnInt>>>, (int, int, int))>();
        foreach (var (g, infos) in BuildExtensions(new HashSet<ConcreteGroup<Ep2<ZnInt, Mat>>>() { q16, c2q8 }))
            l0.Add((g, infos));

        foreach (var (g, infos) in BuildExtensions(allSdp16))
        {
            if (l0.Count + l1.Count == 51)
                break;

            if (l0.Any(g0 => g.IsIsomorphicTo(g0.Item1)))
                continue;

            l1.Add((g, infos));
        }

        foreach (var (g, infos) in BuildExtensions(allAb16))
        {
            if (l0.Any(g0 => g.IsIsomorphicTo(g0.Item1)))
                continue;

            if (l1.Any(g0 => g.IsIsomorphicTo(g0.Item1)))
                continue;

            l2.Add((g, infos));
            if (l0.Count + l1.Count + l2.Count == 51)
                break;
        }

        DisplayInfosGroups(l0.ToArray(), countStart: 0);
        DisplayInfosGroups(l1.ToArray(), countStart: l0.Count);
        DisplayInfosGroups(l2.ToArray(), countStart: l0.Count + l1.Count);

        Console.Beep(); // ~50min
    }

    public static void AllGroupsOrder_12_24()
    {
        var c2 = FG.Abelian(2);
        var c3 = FG.Abelian(3);
        var ab4 = new[] { FG.Abelian(2, 2), FG.Abelian(4) };

        var allOrder8 = BuildExtensions(ab4.ToHashSet())
            .OrderBy(e => e.Item1.GroupType)
            .ThenByDescending(e => e.Item1.ElementsOrders.Values.Max())
            .ThenBy(e => e.Item2).ToArray();

        var tuple12 = ab4.SelectMany(e => new[] { (c3, e), (e, c3) }).ToArray();
        var allOrder12 = BuildExtensions(tuple12)
            .OrderBy(e => e.Item1.GroupType)
            .ThenByDescending(e => e.Item1.ElementsOrders.Values.Max())
            .ThenBy(e => e.Item2).ToArray();

        var tuple24 = allOrder12.Select(e => (e.Item1, c2))
            .Concat(allOrder8.Select(e => (e.Item1, c3)))
            .ToArray();
        var allOrder24 = BuildExtensions(tuple24)
            .OrderBy(e => e.Item1.GroupType)
            .ThenByDescending(e => e.Item1.ElementsOrders.Values.Max())
            .ThenBy(e => e.Item2).ToArray();

        Console.Clear();

        DisplayInfosGroups(allOrder8);
        DisplayInfosGroups(allOrder12);
        DisplayInfosGroups(allOrder24);

        Console.Beep();
    }

    public static void AllGroupsOrder_10_20_40()
    {
        var c5 = FG.Abelian(5);
        var c10 = FG.Abelian(2, 5);
        var ab4 = new[] { FG.Abelian(2, 2), FG.Abelian(4) };
        var tuples20 = ab4.Select(g => (c5, g)).ToArray();
        var tuples40 = ab4.Select(g => (c10, g)).ToArray();

        var allOrder10 = BuildExtensions(new[] { (c5, FG.Abelian(2)) })
            .OrderBy(e => e.Item1.GroupType)
            .ThenByDescending(e => e.Item1.ElementsOrders.Values.Max())
            .ThenBy(e => e.Item2).ToArray();

        var allOrder20 = BuildExtensions(tuples20)
            .OrderBy(e => e.Item1.GroupType)
            .ThenByDescending(e => e.Item1.ElementsOrders.Values.Max())
            .ThenBy(e => e.Item2).ToArray();

        var allOrder40 = BuildExtensions(tuples40)
            .OrderBy(e => e.Item1.GroupType)
            .ThenByDescending(e => e.Item1.ElementsOrders.Values.Max())
            .ThenBy(e => e.Item2).ToArray();

        Console.Clear();

        DisplayInfosGroups(allOrder10);
        DisplayInfosGroups(allOrder20);
        DisplayInfosGroups(allOrder40);

        Console.Beep();
    }

    public static void ExampleAll81Orders()
    {
        var c3 = (ConcreteGroup<ZnInt>)(new Cn(3));
        var allOrder9 = BuildExtensions(new[] { (c3, c3) })
            .OrderBy(e => e.Item1.GroupType)
            .ThenByDescending(e => e.Item1.ElementsOrders.Values.Max())
            .ThenBy(e => e.Item2).ToArray();

        var allOrder27 = BuildExtensions(allOrder9.Select(e => (e.Item1, c3)).ToArray())
            .OrderBy(e => e.Item1.GroupType)
            .ThenByDescending(e => e.Item1.ElementsOrders.Values.Max())
            .ThenBy(e => e.Item2).ToArray();

        (allOrder27[2], allOrder27[3]) = (allOrder27[3], allOrder27[2]);
        var allOrder81 = BuildExtensions(allOrder27.Select(e => (e.Item1, c3)).ToArray())
            .Take(15)
            .OrderBy(e => e.Item1.GroupType)
            .ThenByDescending(e => e.Item1.ElementsOrders.Values.Max())
            .ThenBy(e => e.Item2).ToArray();

        Console.Clear();

        DisplayInfosGroups(allOrder9);
        DisplayInfosGroups(allOrder27);
        DisplayInfosGroups(allOrder81);

        Console.Beep(); // ~30min
    }

    public static void ExampleAll36Order()
    {
        var c2 = FG.Abelian(2);
        var c3 = FG.Abelian(3);
        var ab4 = new[] { FG.Abelian(2, 2), FG.Abelian(4) };
        var ab9 = new[] { FG.Abelian(3, 3), FG.Abelian(9) };

        var tuple12 = ab4.SelectMany(e => new[] { (c3, e), (e, c3) }).ToArray();
        var allOrder12 = BuildExtensions(tuple12)
            .OrderBy(e => e.Item1.GroupType)
            .ThenByDescending(e => e.Item1.ElementsOrders.Values.Max())
            .ThenBy(e => e.Item2).ToArray();
    
        var allOrder18 = BuildExtensions(ab9.ToHashSet())
            .OrderBy(e => e.Item1.GroupType)
            .ThenByDescending(e => e.Item1.ElementsOrders.Values.Max())
            .ThenBy(e => e.Item2).ToArray();

        var ab9a = ab9.Select(e9 => Product.Generate(FG.Abelian(1), e9)).ToList();
        ab9a.ForEach(e => e.SetName(((Gp2<Ep<ZnInt>, Ep<ZnInt>>)e.BaseGroup).G2.Name));
        var tuple36 = ab9a.Select(e9 => (e9, FG.Abelian(4)))
            .Concat(allOrder12.Select(e => (e.Item1, c3)))
            .Concat(allOrder18.Select(e => (e.Item1, c2)))
            .ToArray();
    
        var allOrder36 = BuildExtensions(tuple36)
            // .Take(14) // stop 
            .OrderBy(e => e.Item1.GroupType)
            .ThenByDescending(e => e.Item1.ElementsOrders.Values.Max())
            .ThenBy(e => e.Item2).ToArray();
    
        Console.Clear();

        DisplayInfosGroups(allOrder12);
        DisplayInfosGroups(allOrder18);
        DisplayInfosGroups(allOrder36);

        Console.Beep();
    }
    
    public static void ExampleAll42Order()
    {
        var (c7, c3, c2) = (FG.Abelian(7), FG.Abelian(3), FG.Abelian(2));
        var allOrder14 = BuildExtensions(new[] { (c7, c2) })
            .OrderBy(e => e.Item1.GroupType)
            .ThenByDescending(e => e.Item1.ElementsOrders.Values.Max())
            .ThenBy(e => e.Item2).ToArray();

        var allOrder21 = BuildExtensions(new[] { (c7, c3) })
            .OrderBy(e => e.Item1.GroupType)
            .ThenByDescending(e => e.Item1.ElementsOrders.Values.Max())
            .ThenBy(e => e.Item2).ToArray();

        var tuple42 = allOrder14.Select(e => (e.Item1, c3)).Concat(allOrder21.Select(e => (e.Item1, c2))).ToArray();
        var allOrder42 = BuildExtensions(tuple42)
            .OrderBy(e => e.Item1.GroupType)
            .ThenByDescending(e => e.Item1.ElementsOrders.Values.Max())
            .ThenBy(e => e.Item2).ToArray();

        DisplayInfosGroups(allOrder14);
        DisplayInfosGroups(allOrder21);
        DisplayInfosGroups(allOrder42);
    }

    public static void ExampleAll56Order()
    {
        var allOrder14 = new[] { Group.SemiDirectProd(FG.Abelian(7), FG.Abelian(2)), Product.Generate(FG.Abelian(2), FG.Abelian(7)) }
            .Select(e => (e, SubGroupsDetails(Group.AllSubGroups(e))))
            .OrderBy(e => e.Item1.GroupType)
            .ThenByDescending(e => e.Item1.ElementsOrders.Values.Max())
            .ThenBy(e => e.Item2).ToArray();

        var allOrder28 = BuildExtensions(allOrder14.Select(e => e.Item1).ToHashSet())
            .OrderBy(e => e.Item1.GroupType)
            .ThenByDescending(e => e.Item1.ElementsOrders.Values.Max())
            .ThenBy(e => e.Item2).ToArray();

        var allOrder56 = BuildExtensions(allOrder28.Select(e => e.Item1).ToHashSet())
            .OrderBy(e => e.Item1.GroupType)
            .ThenByDescending(e => e.Item1.ElementsOrders.Values.Max())
            .ThenBy(e => e.Item2).ToArray();

        DisplayInfosGroups(allOrder14);
        DisplayInfosGroups(allOrder28);
        DisplayInfosGroups(allOrder56);

        Console.Beep();
    }

    public static void ExampleAll60Order()
    {
        var ab4 = new[] { Product.Generate(new Cn(2), new Cn(2)), Product.Generate("C4", new Cn(1), new Cn(4)) };
        var c15 = Product.Generate("C15", FG.Abelian(1), Product.Generate(new Cn(3), new Cn(5)));
        var allOrder20 = BuildExtensions(
                ab4.Select(e => (FG.Abelian(5), e))
                    .Append((FG.Abelian(10), Product.Generate("C2", new Cn(1), new Cn(2)))).ToArray())
            .OrderBy(e => e.Item1.GroupType)
            .ThenByDescending(e => e.Item1.ElementsOrders.Values.Max())
            .ThenBy(e => e.Item2).ToArray();
    
        var tuple6 = new[] { FG.DihedralSdp(3), Product.Generate(new Cn(2), new Cn(3)) };
        var allOrder30 = BuildExtensions(tuple6.Select(e => (FG.Abelian(5), e)).ToArray())
            .OrderBy(e => e.Item1.GroupType)
            .ThenByDescending(e => e.Item1.ElementsOrders.Values.Max())
            .ThenBy(e => e.Item2).ToArray();

        var tuple60 = allOrder20.Select(e => (e.Item1, FG.Abelian(3)))
            .Concat(allOrder30.Select(e => (e.Item1, FG.Abelian(2))))
            .Append((c15, FG.Abelian(4))).ToArray();
        var allOrder60 = BuildExtensions(tuple60)
            .OrderBy(e => e.Item1.GroupType)
            .ThenByDescending(e => e.Item1.ElementsOrders.Values.Max())
            .ThenBy(e => e.Item2).ToArray();

        DisplayInfosGroups(allOrder20);
        DisplayInfosGroups(allOrder30);
        DisplayInfosGroups(allOrder60); // Except A5 which is Simple

        Console.Beep();
    }

}
using System.Numerics;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Perms;

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

    public class TwoCocyclesDFS<Tn, Tg> where Tn : struct, IElt<Tn> where Tg : struct, IElt<Tg>
    {
        public ConcreteGroup<Tg> G { get; }
        public ConcreteGroup<Ep<Tg>> GxG { get; }
        public ConcreteGroup<Tn> N { get; }
        public MapElt<Tg, Automorphism<Tn>> L { get; }
        public HashSet<MapElt<Ep<Tg>, Tn>> AllMapsTwo { get; private set; }
        public MapGroupBase<Ep<Tg>, Tn> Group_GxG_N { get; }
        public string Lbl { get; }

        public TwoCocyclesDFS(ConcreteGroup<Tn> N0, ConcreteGroup<Tg> G0, MapElt<Tg, Automorphism<Tn>> L0, string lbl)
        {
            (N, G, L) = (N0, G0, L0);
            GxG = Product.GpGenerate(G, G);
            Group_GxG_N = new MapGroupBase<Ep<Tg>, Tn>(GxG, N);
            AllMapsTwo = new() { Group_GxG_N.Neutral() };
            Lbl = lbl;
        }

        public HashSet<MapElt<Ep<Tg>, Tn>> AllTwoCocycles()
        {
            AllMapsTwo.Clear();
            var nG = G.Neutral();
            var map = G.SelectMany(e => new[] { new Ep<Tg>(nG, e), new Ep<Tg>(e, nG) }).Distinct()
                .ToDictionary(e => e, _ => N.Neutral());
            var contents = All2Cocycles(map).ToHashSet();
            Console.WriteLine("############# COMPLETE ##############");
            return contents;
        }

        private IEnumerable<MapElt<Ep<Tg>, Tn>> All2Cocycles(Dictionary<Ep<Tg>, Tn> map)
        {
            var nGxG_N = Group_GxG_N.Neutral();
            AllMapsTwo.Add(nGxG_N);
            yield return nGxG_N;

            var lG = G.GetGenerators().Where(e => !e.Equals(G.Neutral())).Distinct().ToArray();
            var startG = G.GetGenerators().Grid2D(lG).Select(e => new Ep<Tg>(e.t1, e.t2))
                .OrderBy(e => (e[0], G.Op(e[0], e[1]))).ToArray();
            var nbStartG = startG.Length;
            var startN = Group.AllSubGroups(N).Keys.SelectMany(sg => sg.GetGenerators()).Prepend(N.Neutral()).Distinct().ToArray();
            startN.Println("startN");
            startG.Println("startG");
            var max = BigInteger.Pow(startN.Length, nbStartG);
            Console.WriteLine($"WARNING !!!!!!!!!!!!! Max Iterations : {max}; StartN:{startN.Length} StartG:{nbStartG}");
            var ct = 0;
            foreach (var l in startN.MultiLoop(nbStartG))
            {
                Console.WriteLine($"{Lbl} search : {++ct} Max Iterations : {max}; StartN:{startN.Length} StartG:{nbStartG}");
                var pmap1 = startG.Zip(l).ToArray();

                var map2 = map.Select(e => (e.Key, e.Value)).Concat(pmap1).Distinct().ToDictionary(e => e.Item1, e => e.Item2);
                if (pmap1.Grid2D(pmap1).Any(e => !TwoCocycleCondition(map2, e.t1.Item1, e.t2.Item1)))
                {
                    var first = pmap1.Grid2D(pmap1).First(e => !TwoCocycleCondition(map2, e.t1.Item1, e.t2.Item1));
                    Console.WriteLine("Reject : {0}", pmap1.Glue(";"));
                    continue;
                }

                var setKeys = map2.Keys.ToHashSet();
                if (AllMapsTwo.Any(map3 => SubMap(map3, setKeys, map2)))
                    continue;

                Console.WriteLine("PASS   : {0}", pmap1.Glue(";"));
                var nbSolsBefore = AllMapsTwo.Count;
                foreach (var sol in SearchTwoCocycles(map2, startN))
                {
                    Console.Write('$');
                    yield return sol;
                }

                var nbSolsAfter = AllMapsTwo.Count;
                if (nbSolsBefore != nbSolsAfter)
                {
                    Console.WriteLine($" AllSols:{AllMapsTwo.Count}");
                }
            }
        }

        public int ct = 0;

        private IEnumerable<MapElt<Ep<Tg>, Tn>> SearchTwoCocycles(Dictionary<Ep<Tg>, Tn> current, Tn[] startN)
        {
            if (current.Count == GxG.Count())
            {
                var sol = new MapElt<Ep<Tg>, Tn>(GxG, N, current);
                if (!AllMapsTwo.Contains(sol) && ValidMap(sol))
                {
                    var prevMaps = AllMapsTwo.ToHashSet();
                    AllMapsTwo.Add(sol);
                    yield return sol;
                    
                    // Z2 group of 2cocycles structure
                    var newSols = Group.GenerateElements(Group_GxG_N, new[] { sol }); 
                    foreach (var (e0, e1) in newSols.Grid2D(prevMaps))
                    {
                        var sol1 = Group_GxG_N.Op(e0, e1);
                        if (!AllMapsTwo.Contains(sol1) && ValidMap(sol1))
                        {
                            AllMapsTwo.Add(sol1);
                            yield return sol1;
                        }
                    }
                }
            }
            else
            {
                foreach (var next in NextChildsTwoCocycles(current, startN))
                {
                    var allSols = AllMapsTwo.Count;
                    foreach (var sol in SearchTwoCocycles(next, startN))
                    {
                        // Console.WriteLine($"Current: {current.Count} then Next:{next.Count} then Next2:{sol.map.Count}");
                        yield return sol;
                    }

                    if (allSols != AllMapsTwo.Count)
                        break;
                }
            }
        }

        public IEnumerable<Dictionary<Ep<Tg>, Tn>> NextChildsTwoCocycles(Dictionary<Ep<Tg>, Tn> current, Tn[] startN)
        {
            var nG = G.Neutral();
            var candidates = current.Keys.Where(e => !e[0].Equals(nG) && !e[1].Equals(nG)).SelectMany(e => new[] { e[0], e[1] })
                .ToHashSet();
            var rems = GxG.Except(current.Keys)
                .Where(e => !G.Op(e[0], e[1]).Equals(nG))
                .OrderByDescending(e =>
                    (e[0].Equals(e[1]) ? 1 : 0) + (candidates.Contains(e[0]) ? 1 : 0) + (candidates.Contains(e[1]) ? 1 : 0))
                .ThenByDescending(e => GxG.ElementsOrders[e])
                .ThenAscending().ToArray();
            var nb = current.Count;
            var pos = rems.Grid2D(startN.OrderByDescending(e => N.ElementsOrders[e]))
                .OrderBy(e => (e.t1[0], G.Op(e.t1[0], e.t1[1]))).ToArray();

            // Console.WriteLine($"WARNING2 !!!!!!!!!!!!! nb Iterations : {pos.Length} RemG:{rems.Length}");
            foreach (var (g0, n0) in pos)
            {
                var next = TwoCocyclesUpdate(current, new[] { (g0[0], g0[1], n0) }.ToList());
                if (next.Count == nb)
                    continue;
                else if (next.Count > nb)
                    yield return next;
                else
                    throw new();
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

        private static bool SubMap(MapElt<Ep<Tg>, Tn> map, HashSet<Ep<Tg>> keySubmap, Dictionary<Ep<Tg>, Tn> subMap)
        {
            if (map.map.Count < subMap.Count)
                return false;

            if (!keySubmap.IsSubsetOf(map.map.Keys))
                return false;

            return keySubmap.All(k => subMap[k].Equals(map[k]));
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

            var g1g2g3 = G.Grid3D(G, G).ToArray();
            var setKeys = next.Keys.ToHashSet();
            if (AllMapsTwo.Any(map => SubMap(map, setKeys, next)))
                return prev;

            foreach (var (g1, g2, g3) in g1g2g3)
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
                    var ng2_g3 = L[g1].Invert()[N.Op(N.Invert(ng1_g2g3), N.Op(ng1_g2, ng1g2_g3))];
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

    static Dictionary<MapElt<Tg, Automorphism<Tn>>, HashSet<MapElt<Ep<Tg>, Tn>>> All_2_Cocycles_N_G<Tn, Tg>(ConcreteGroup<Tn> N,
        ConcreteGroup<Tg> G, bool trivialActionOnly = true)
        where Tg : struct, IElt<Tg>
        where Tn : struct, IElt<Tn>
    {
        var autN = Group.AutomorphismGroup(N);
        var allOps = Group.AllHomomorphisms(G, autN);
        Console.WriteLine($"N:{N.ShortName} and G:{G.ShortName}");
        Ring.MatrixDisplayForm = Ring.MatrixDisplay.SquareBracket;
        var all = new Dictionary<MapElt<Tg, Automorphism<Tn>>, HashSet<MapElt<Ep<Tg>, Tn>>>();
        var lbl = 0;
        if (trivialActionOnly)
            allOps = new() { new(G, G.ToDictionary(e => e, _ => autN.Neutral())) };

        foreach (var L in allOps)
        {
            var L0 = new MapElt<Tg, Automorphism<Tn>>(G, autN, new(L.HomMap));
            var homol = new TwoCocyclesDFS<Tn, Tg>(N, G, L0, $"Lbl{++lbl}/{allOps.Count}");
            var all2Cocycles = homol.AllTwoCocycles();
            Console.WriteLine($"L:{L0}");
            Console.WriteLine($"Count:{all2Cocycles.Count}");
            all[L0] = all2Cocycles.ToHashSet();
            // foreach (var co in all2Cocycles)
            //     Console.WriteLine($"  2co:{co}");

            Console.WriteLine();
        }

        Console.WriteLine($"Total:{all.Values.Sum(v => v.Count)}");
        Console.WriteLine($"N:{N.ShortName} and G:{G.ShortName}");
        Console.WriteLine("Press key to continue...");
        Console.ReadLine();
        return all;
    }

    public static void TwoCocyclesExamples()
    {
        var (c2, c4, c2c2, c8, c2c2c2, c4c4, c4c2, d8, q8, c16) = (new Cn(2), new Cn(4), FG.Abelian(2, 2), new Cn(8),
            FG.Abelian(2, 2, 2), FG.Abelian(4, 4), FG.Abelian(4, 2), FG.Dihedral(4), FG.Quaternion(8), new Cn(16));

        // Twisted Actions for trivial action of G by N
        All_2_Cocycles_N_G(c4, c8);
        All_2_Cocycles_N_G(c2c2, c8);
        All_2_Cocycles_N_G(c4, c4c2);
        All_2_Cocycles_N_G(c4, d8);
        All_2_Cocycles_N_G(c4, q8);
        All_2_Cocycles_N_G(c2c2, c4c2);
        All_2_Cocycles_N_G(c2c2, d8);
        // C2 x C2 x C2 is problematic
    }

    public static void NonSplitExample_C4_D8()
    {
        var (c4, d8) = (new Cn(4), FG.Dihedral(4));
        var all = All_2_Cocycles_N_G(c4, d8, trivialActionOnly: false);
        var allNonSplit = new HashSet<ExtensionGroup<ZnInt, Perm>>(new IsomorphEquality<Ep2<ZnInt, Perm>>());
        foreach (var (L, tws) in all)
        {
            foreach (var w in tws)
            {
                var ext = Group.ExtensionGroup(c4, L, w, d8);
                var sp = NonSplitExtension.AllSplittingGroups(c4, ext, d8, details: false).First();
                if (sp.s.IsNull)
                    allNonSplit.Add(ext);
            }
        }

        foreach (var ext in allNonSplit)
        {
            DisplayGroup.HeadOrders(ext);
        }

        Console.WriteLine($"Total Non Split Ext : {allNonSplit.Count}");
        // Results of 51 non isometrics groups of order 32, all are NonAbelianGroup, it is invalid. 
    }
}
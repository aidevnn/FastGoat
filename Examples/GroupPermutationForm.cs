using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.Subgroups;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Perms;

namespace Examples;

public static class GroupPermutationForm
{
    static GroupPermutationForm()
    {
        Perm.Style = DisplayPerm.CyclesComma;
    }

    static (int m, int n, int r)[] MetaCyclicSdp(int order)
    {
        return IntExt.Dividors(order).Where(d => d > 1)
            .SelectMany(m => FG.MetaCyclicSdpGetR(m, order / m).Select(r => (m, n: order / m, r)))
            .Where(e => e.r > 1)
            .ToArray();
    }

    static Perm ConcatPerm(params Perm[] perms)
    {
        var Ns = perms.Select(p => p.Sn.N).ToArray();
        var sn = new Sn(Ns.Sum());
        var cum = Ns.Aggregate(new[] { 0 }, (acc, e) => acc.Append(acc.Last() + e).ToArray());
        var table = perms.Zip(cum).SelectMany(e => e.First.Table.Select(i => i + e.Second)).ToArray();
        return sn.CreateElementTable(table);
    }

    static Perm PaddingRight(Perm perm, int pad) => pad == 0 ? perm : ConcatPerm(perm, (new Sn(pad)).Neutral());
    static Perm PaddingLeft(Perm perm, int pad) => pad == 0 ? perm : ConcatPerm((new Sn(pad)).Neutral(), perm);
    static Perm Padding(int padLeft, Perm perm, int padRight) => PaddingLeft(PaddingRight(perm, padRight), padLeft);

    static Perm[] CyclesSplit(int m)
    {
        var dec = IntExt.PrimesDec(m).Select(e => e.Key.Pow(e.Value)).Order().ToArray();
        return dec.Select(e => new Sn(e).Cycle(e.Range(1))).ToArray();
    }

    static Perm Cycles(int m)
    {
        if (m == 1)
            return new Sn(1).Neutral();

        return ConcatPerm(CyclesSplit(m));
    }

    static (string name, Sn sn, Perm[] gens) AbelianPermGens(params int[] seq)
    {
        var cycles = seq.Select(c => Cycles(c)).ToArray();
        var dims = cycles.Select(a => a.Sn.N).ToArray();
        var name = cycles.Select(c => c.Order).ToAbString();
        var gens = cycles.Select((a, i) => Padding(dims.Take(i).Sum(), a, dims.Skip(i + 1).Sum())).ToArray();
        return (name, gens[0].Sn, gens);
    }

    static ConcreteGroup<Perm> AbelianPerm(params int[] seq)
    {
        var (name, sn, gens) = AbelianPermGens(seq);
        return Group.Generate(name, sn, gens.ToArray());
    }

    static (Perm a, Perm b) rUmToPerm(int m, int r)
    {
        var (seqA, seqB) = (new List<Perm>(), new List<Perm>());
        foreach (var c in CyclesSplit(m))
        {
            var N = c.Sn.N;
            seqA.Add(c);
            if (IntExt.Gcd(N, r) != 1)
            {
                seqB.Add(c.Sn.Neutral());
                continue;
            }

            var ri = new ZnInt(N, r).Inv();
            var b = c.Sn.CreateElementTable(N.Range().Select(i => (ri * i).K).ToArray());
            seqB.Add(b);
        }

        return (ConcatPerm(seqA.ToArray()), ConcatPerm(seqB.ToArray()));
    }

    static IEnumerable<(Perm a, Perm b)> OpsCnByUm(int n, Perm a, Perm b)
    {
        foreach (var perms in CyclesSplit(n).AllCombinations().Where(e => e.Length > 0))
        {
            var table = IntExt.PermAndCyclesFromType(perms.SelectMany(e => e.PermType).Order().ToArray()).perm;
            var b0 = new Sn(table.Length).CreateElementTable(table);
            var a1 = PaddingRight(a, b0.Sn.N);
            var b1 = ConcatPerm(b, b0);
            yield return (a1, b1);
        }
    }

    static (string name, Sn sn, Perm[] gens) MetaCyclicGens(int m, int n, int r)
    {
        var (a, b) = rUmToPerm(m, r);
        var sn = a.Sn;
        var ar = a ^ r;

        var ord = Group.ElementIsOrder(sn, b, n);
        var c1 = (b ^ n).Equals(sn.Neutral());
        var c2 = (b * a * (b ^ -1)).Equals(ar);
        if (!ord && c1 && c2)
        {
            foreach (var (a0, b0) in OpsCnByUm(n, a, b))
            {
                var H1 = Group.GenerateElements(a0.Sn, a0, b0);
                if (H1.Count == m * n)
                {
                    (sn, a, ar, b) = (a0.Sn, a0, a0 ^ r, b0);
                    break;
                }
            }

            ord = Group.ElementIsOrder(sn, b, n);
            c1 = (b ^ n).Equals(sn.Neutral());
            c2 = (b * a * (b ^ -1)).Equals(ar);
        }

        if (ord && c1 && c2)
        {
            var name = FG.MetaCyclicName(m, n, r);
            return (name, sn, [a, b]);
        }

        return ("C1", sn, [sn.Neutral(), sn.Neutral()]);
    }

    static ConcreteGroup<Perm> SearchPermGroupMetaCyclic(int m, int n, int r)
    {
        var (name, sn, gens) = MetaCyclicGens(m, n, r);
        return Group.Generate(name, sn, gens);
    }

    static (string name, Sn sn, Perm[] gens) DihedralGens(int n) => MetaCyclicGens(n, 2, n - 1);
    static ConcreteGroup<Perm> Dihedral(int n) => SearchPermGroupMetaCyclic(n, 2, n - 1);
    static ConcreteGroup<Perm> SemiDihedral(int n) => SearchPermGroupMetaCyclic(1 << (n - 1), 2, (1 << (n - 2)) - 1);
    static ConcreteGroup<Perm> ModularMax(int n) => SearchPermGroupMetaCyclic(1 << (n - 1), 2, (1 << (n - 2)) + 1);

    static (string name, Sn sn, Perm[] gens) QuaternionGens(int n)
    {
        if (!int.IsPow2(n))
            throw new();

        var k1 = n / 2;
        var k2 = n / 4;
        var sn = new Sn(n);
        var a = sn.Cycle(k1.Range(1)) * sn.Cycle(k1.Range(k1 + 1));
        var b = Group.OpSeq(sn, k2.SeqLazy().Select(i => sn.Cycle(i + 1, n - i, i + 1 + k2, n - i - k2)));
        return ($"Q{n}", sn, [a, b]);
    }

    static (string name, Sn sn, Perm[] gens) DicyclicGens(int n)
    {
        var k = IntExt.PrimesDecomposition(n).Count(i => i == 2);
        if (k != 0)
        {
            var m = n / (1 << k);
            var (Qname, sn0, gens) = QuaternionGens(1 << (k + 2));
            if (m == 1)
                return (Qname, sn0, gens);

            var (b, c) = gens.Deconstruct();
            var (a0, b0) = rUmToPerm(m, m - 1);
            var a1 = PaddingRight(a0, b.Sn.N);
            var b1 = PaddingLeft(b, b0.Sn.N);
            var c1 = ConcatPerm(b0, c);

            var H1 = Group.GenerateElements(a1.Sn, a1, b1, c1);
            if (H1.Count == 4 * n)
                return ($"C{m} x: {Qname}", a1.Sn, [a1, b1, c1]);

            throw new();
        }

        var (b2, c2) = rUmToPerm(n, n - 1);
        foreach (var (b3, c3) in OpsCnByUm(n: 4, b2, c2))
        {
            var H1 = Group.GenerateElements(b3.Sn, b3, c3);
            if (H1.Count == 4 * n)
                return ($"Dic{n}", b3.Sn, [b3, c3]);
        }

        throw new();
    }

    static HashSet<T> GenerateElementsLimited<T>(IGroup<T> bg, T[] generators, int limits)
        where T : struct, IElt<T>
    {
        var elements = new HashSet<T>([bg.Neutral()]);
        var q = new Queue<T>(elements);
        HashSet<T> generatedElements = new HashSet<T>(elements);
        while (q.Count != 0 && limits >= generatedElements.Count)
        {
            var e1 = q.Dequeue();
            foreach (var e2 in generators)
            {
                var e3 = bg.Op(e2, e1);
                if (generatedElements.Add(e3))
                    q.Enqueue(e3);

                if (limits < generatedElements.Count)
                    return [];
            }
        }

        if (limits < generatedElements.Count)
            return [];

        return generatedElements;
    }

    static Dictionary<ConcreteGroup<T>, List<GroupSubset<T>>> RestrictionConjugates<T>(ConcreteGroup<T> g,
        List<GroupSubset<T>> conjugates) where T : struct, IElt<T>
    {
        var lt = conjugates.Where(h => h.SubSetOf(g)).ToHashSet();
        var all = new Dictionary<ConcreteGroup<T>, List<GroupSubset<T>>>(conjugates.Count);
        var gens = g.GetGenerators().ToHashSet();
        var act = Group.ByConjugateSet(g);
        while (lt.Count != 0)
        {
            var sg = lt.First();
            var conjs = Group.Orbits(gens, act, sg);
            var subConjs = lt.Where(c0 => conjs.Any(c1 => c0.SetEquals(c1))).ToList();
            var g1 = Group.Generate("H", g, sg.Generators.ToArray());
            all[g1] = subConjs;
            lt.ExceptWith(subConjs);
        }

        return all;
    }

    static Dictionary<ConcreteGroup<T>, List<GroupSubset<T>>> Restriction<T>(ConcreteGroup<T> G,
        Dictionary<ConcreteGroup<T>, List<GroupSubset<T>>> allSubgroups) where T : struct, IElt<T>
    {
        var setG = G.ToSet();
        return allSubgroups.Where(e => e.Value.Any(f => f.SubSetOf(setG)))
            .SelectMany(e => RestrictionConjugates(G, e.Value)).ToDictionary();
    }

    static ConcreteGroup<Perm> PermutationForm<T>(Dictionary<ConcreteGroup<T>, List<GroupSubset<T>>> subGroups)
        where T : struct, IElt<T>
    {
        var G = subGroups.MaxBy(g => g.Key.Order).Key;
        if (G.GroupType == GroupType.AbelianGroup)
            return AbelianPerm(Group.AbelianGroupType(G));

        var ordG = G.Order;
        var byOrder = subGroups.Keys.GroupBy(e => e.Order).ToDictionary(e => e.Key, e => e.ToHashSet());
        var (H, Facts) = subGroups.Keys
            .Where(e => e.Order != 1 && e.Order != ordG && subGroups[e].Count == 1).OrderBy(g => g.Order).ToArray()
            .ToDictionary(d => d, d =>
            {
                if (byOrder.ContainsKey(ordG / d.Order))
                    return byOrder[ordG / d.Order].Where(e => e.Intersect(d).Count() == 1).Take(1).ToList();

                Console.WriteLine($"############ Warning {G.ShortName} D:{d.ShortName}");
                return [];
            })
            .Where(e => e.Value.Count > 0)
            .OrderByDescending(f => f.Value.Take(1).Count(g => subGroups[g].Count == 1))
            .ThenByDescending(f => f.Value.Take(1).Sum(g => g.Order))
            .ThenBy(f => f.Key.Order)
            .FirstOrDefault(_ => true, new(G, new()));
        if (Facts.Count == 0)
        {
            if (G.Order == 60)
            {
                var A5 = FG.Alternate(5);
                if (G.IsIsomorphicTo(A5))
                    return A5;
            }

            return G.ToPermGroup().Item1;
        }

        var K = Facts.First();
        if (Group.IsNormalSubgroup(G, K))
        {
            var h = PermutationForm(Restriction(H, subGroups));
            var k = PermutationForm(Restriction(K, subGroups));
            var nh = h.Neutral().Sn.N;
            var nk = k.Neutral().Sn.N;
            var hGens = h.GetGenerators().Select(e => PaddingRight(e, nk));
            var kGens = k.GetGenerators().Select(e => PaddingLeft(e, nh));
            var sn = new Sn(nh + nk);
            var HK = Group.Generate(G.Name, sn, hGens.Concat(kGens).ToArray());
            if (!HK.IsIsomorphicTo(G))
                throw new($"###1 {G.ShortName}");

            return HK;
        }
        else
        {
            var k = PermutationForm(Restriction(K, subGroups));
            var h = PermutationForm(Restriction(H, subGroups)).ToPermGroup().Item1;
            var (snH, snK) = (h.Neutral().Sn, k.Neutral().Sn);
            var (nH, nK) = (snH.N, snK.N);

            var autHpg = FG.RegPermAutGroup(Group.AutomorphismGroup(h)).Item1;
            var homKAutH = Group.AllHomomorphisms(k, autHpg);
            foreach (var hom in homKAutH.OrderBy(e => e.Kernel().Count()))
            {
                var img = hom.Image().ToArray();
                if (img.Length == k.Count())
                {
                    var HK0 = Group.Generate(G.Name, snH, h.GetGenerators().Concat(img).ToArray());
                    if (HK0.IsIsomorphicTo(G))
                        return HK0;
                }

                var hGens = h.GetGenerators().Select(f => PaddingRight(f, nK)).ToArray();
                var kGens = k.GetGenerators().Select(e => ConcatPerm(hom[e], e)).ToArray();
                var sn = new Sn(nH + nK);
                var hkGens = hGens.Concat(kGens).ToArray();
                if (GenerateElementsLimited(sn, hkGens, ordG).Count != ordG)
                    continue;

                var HK = Group.Generate(G.Name, sn, hkGens);
                if (HK.IsIsomorphicTo(G))
                    return HK;
            }

            Console.WriteLine($"NOT FOUND  H:{H.ShortName} K:{K.ShortName} {autHpg.ShortName}");
            Console.WriteLine($"H in {h.Neutral().Sn} K in {k.Neutral().Sn}");
            throw new($"###2 {G.ShortName}");
        }
    }

    public static void ExampleAbelian()
    {
        for (int ord = 1; ord <= 128; ord++)
        {
            foreach (var abType in AbelianExt.AllAbTypes(ord))
            {
                var ab = AbelianPerm(abType.can);
                DisplayGroup.HeadGenerators(ab);
                var abType2 = Group.AbelianGroupType(ab);
                if (!abType.can.SequenceEqual(abType2))
                    throw new();
            }
        }
    }

    public static void ExampleDihedral()
    {
        for (int n = 3; n < 33; n++)
        {
            var D2pg = Dihedral(n);
            var D2n = FG.Dihedral(n);
            DisplayGroup.HeadGenerators(D2pg);
            DisplayGroup.AreIsomorphics(D2pg, D2n);
            Console.WriteLine();
            if (!D2n.IsIsomorphicTo(D2pg))
                throw new();
        }
    }

    public static void ExampleDicyclic()
    {
        for (int m = 3; m < 33; m++)
        {
            var Dic_m = FG.DicyclicGL2p(m);
            var (name, sn, gens) = DicyclicGens(m);
            var Dic_m_pg = Group.Generate(name + "pg", sn, gens);
            DisplayGroup.HeadGenerators(Dic_m_pg);
            DisplayGroup.AreIsomorphics(Dic_m_pg, Dic_m);
            Console.WriteLine();
            if (!Dic_m_pg.IsIsomorphicTo(Dic_m))
                throw new();
        }
    }

    public static void ExampleSemiDihedralAndModularMax()
    {
        for (int n = 4; n < 9; n++)
        {
            var qdpg = SemiDihedral(n);
            DisplayGroup.HeadGenerators(qdpg);
            if (!qdpg.IsIsomorphicTo(FG.SemiDihedralGL2p(n)))
                throw new();

            var mmpg = ModularMax(n);
            DisplayGroup.HeadGenerators(mmpg);
            if (!mmpg.IsIsomorphicTo(FG.ModularMaxGL2p(n)))
                throw new();

            Console.WriteLine();
        }
    }

    public static void ExampleAllMetaCyclicSemiDirectProducts()
    {
        GlobalStopWatch.Restart();
        for (int i = 2; i < 129; i++)
        {
            foreach (var (m, n, r) in MetaCyclicSdp(i))
            {
                var mtpg = SearchPermGroupMetaCyclic(m, n, r);
                var mt = FG.MetaCyclicSdp(m, n, r);
                DisplayGroup.HeadGenerators(mtpg);
                if (!mt.IsIsomorphicTo(mtpg))
                    throw new();
            }
        }

        GlobalStopWatch.Show();
    }

    static void MetaCyclicPermFormUpTo(int maxOrd)
    {
        var (errors, total) = (0, 0);
        GlobalStopWatch.Restart();

        for (int ord = 1; ord <= maxOrd; ord++)
        {
            foreach (var (m, n, r) in MetaCyclicSdp(ord))
            {
                ++total;
                var Gpg = PermutationForm(Group.AllSubGroups(FG.MetaCyclicSdpWg(m, n, r))
                    .ToDictionary(e => e.Key, e => e.Value.Select(f => f.ToSet()).ToList()));
                DisplayGroup.HeadOrders(Gpg);
                DisplayGroup.Generators(Gpg, showBaseGroup: true);
                var sn = Gpg.Neutral().Sn;
                var sn0 = MetaCyclicGens(m, n, r).sn;
                var test = sn.N == sn0.N;
                Console.WriteLine($"{Gpg.ShortName} in {sn} expected {sn0} Test:{test}");
                if (!test)
                    ++errors;

                Console.WriteLine();
            }
        }

        GlobalStopWatch.Show($"Errors:{errors}/{total} Max Order:{maxOrd}");
    }

    static void AllGroupPermFormUpTo(int maxOrd)
    {
        var total = 0;
        GlobalStopWatch.Restart();

        for (int ord = 1; ord <= maxOrd; ord++)
        {
            foreach (var G in FG.AllGroupsOfOrder(ord))
            {
                ++total;
                var Gpg = PermutationForm(Group.AllSubGroups(G)
                    .ToDictionary(e => e.Key, e => e.Value.Select(f => f.ToSet()).ToList()));
                DisplayGroup.HeadGenerators(Gpg);
                var sn = Gpg.Neutral().Sn;
                Console.WriteLine($"{G.ShortName} in {sn} Ratio:{G.Count() / (sn.N + 0.0):f3}");
                Console.WriteLine();
            }
        }

        GlobalStopWatch.Show($"Total:{total}  Max Order:{maxOrd}");
    }

    public static void ExampleAllMetaCyclicPermForm()
    {
        MetaCyclicPermFormUpTo(128);
        // TODO: fix errors
        // |F(7x:12)3| = 84 in S14 expected S11 Test:False
        // |M(9x:12)2| = 108 in S16 expected S13 Test:False
        // |F(7x:18)3| = 126 in S18 expected S16 Test:False
        // 
    }

    public static void ExampleAllGroupsPermFormUpTo63()
    {
        AllGroupPermFormUpTo(63);
    }
}
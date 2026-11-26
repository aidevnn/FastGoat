using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.Naming;
using FastGoat.Structures.Subgroups;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Perms;
using FastGoat.UserGroup.Words;

namespace Craft.Craft;

public struct PolyCyclicSequence<T> : IElt<PolyCyclicSequence<T>> where T : struct, IElt<T>
{
    public ConcreteGroup<T> PCGroup { get; }
    public List<ConcreteGroup<T>> Normals { get; }
    public List<(ConcreteGroup<Coset<T>> K, T e, bool Split)> Generators { get; }
    public int[] RX => Generators.Select(e => e.K.Order).ToArray();

    public PolyCyclicSequence(ConcreteGroup<T> pcGroup, List<(ConcreteGroup<Coset<T>> K, T e, bool Split)> gens)
    {
        Generators = gens;
        Normals = gens.Select(e => e.K.Neutral().H).Prepend(pcGroup).ToList();
        PCGroup = pcGroup;
    }

    public bool Equals(PolyCyclicSequence<T> other) =>
        Generators.Count == other.Generators.Count
        && Generators.Select(e => (e.K.Order, e.Split)).SequenceEqual(other.Generators.Select(e => (e.K.Order, e.Split)))
        && Normals.Last().SetEquals(other.Normals.Last());

    public int CompareTo(PolyCyclicSequence<T> other)
    {
        var compCount = Generators.Count.CompareTo(other.Generators.Count);
        if (compCount != 0)
            return -compCount;

        return Generators.Select(e => (e.K.Order, e.Split)).SequenceCompareTo(other.Generators.Select(e => (e.K.Order, e.Split)));
    }

    public int Hash => Generators.Count;

    public override int GetHashCode() => Hash;

    public PolyCyclicSequence<T> Add(ConcreteGroup<Coset<T>> K, T e, bool Split) =>
        new(PCGroup, Generators.Append((K, e, Split)).ToList());

    public void Display()
    {
        var G = PCGroup;
        Generators.Select(e => (e.e, G.ElementsOrders[e.e], e.K.Order, e.Split, e.K.GetGenerators().Select(c => c.X).ToXArray()))
            .Println($"{G.Name} = {Normals.Select((N, i) => $"G{i}[{N.Order}]").Glue(" > ")}");
        Console.WriteLine();
    }

    public Dictionary<T, (T e, int pow)[]> ExponentVectors()
    {
        var G = PCGroup;
        var expVec = Generators.Select(e => e.K.Order.SeqLazy().Select(pow => (e.e, pow, ep: G.Times(e.e, pow))).ToArray())
            .MultiLoop().Select(l => l.ToArray())
            .Select(l => (G.OpSeq(l.Select(e => e.ep)), l.Select(e => (e.e, e.pow)).ToArray()))
            .DistinctBy(e => e.Item1)
            .ToDictionary(e => e.Item1, e => e.Item2);
        if (!G.SetEquals(expVec.Keys))
        {
            expVec.Println(e => $"{e.Key} = {e.Value.ToXArray()}");
            throw new("ExpVec problem");
        }
        return expVec;
    }

    public PolyCyclicSequence<T> InducedPolyCyclicSequence(ConcreteGroup<T> H)
    {
        var G = PCGroup;
        var R = RX;
        var Y = InducedPolyCyclicSequenceStep(H);
        var expVec = ExponentVectors();
        var X = Y.Where(e => !e.Equals(G.Neutral())).Select(e => (e, R[Depth(expVec, e)])).ToArray();
        return Create(H, X);
    }

    public T[] InducedPolyCyclicSequenceStep(ConcreteGroup<T> H)
    {
        var expVec = ExponentVectors();
        var R = RX;
        var G = PCGroup;
        var Z = Generators.Select(_ => G.Neutral()).ToArray();
        var n = Z.Length;
        var gensH = new HashSet<T>(H.GetGenerators());
        GlobalStopWatch.InfiniteLoopBreakerReset();
        while (gensH.Any())
        {
            GlobalStopWatch.InfiniteLoopBreaker(50);
            var g = gensH.OrderBy(e => Depth(expVec, e)).ThenBy(e => Led(expVec, e)).First();
            gensH.Remove(g);
            var h = Sift(G, expVec, R, Z, g);
            var d = Depth(expVec, h);
            if (d < n)
            {
                gensH.Add(G.Times(h, R[d]));
                gensH.UnionWith(Z.Select(z => G.Commute(h, z)));
                Z[d] = h;
            }

            gensH.Remove(G.Neutral());
        }

        var Y = CanonicalInducedPolyCyclicSeq(G, expVec, R, Z.Where(z => !z.Equals(G.Neutral())).ToArray());
        return Y.ToArray();
    }

    public IEnumerable<(PolyCyclicSequence<T> N, PolyCyclicSequence<T> K, bool isDirectProduct)> Split()
    {
        var G = PCGroup;
        foreach (var N in Normals.Skip(1).SkipLast(1).Reverse())
        {
            if (!Group.IsNormalSubgroup(G, N))
                continue;

            var lhs = InducedPolyCyclicSequence(N);
            var K1 = Group.FactorGroup(G, N);
            if (K1.Order * N.Order == G.Order && K1.IntersectCount(N) == 1)
            {
                var isDirectProduct = Group.IsNormalSubgroup(G, K1);
                var rhs = InducedPolyCyclicSequence(K1);
                yield return (lhs, rhs, isDirectProduct);
            }
        }
    }

    // public (PolyCyclicSequence<Word>, Homomorphism<T, Word>, Homomorphism<Word, T>) ToWordGroup()
    // {
    //     var wg = FG.WordGroup(PCGroup.Name, Relators().ToLower());
    //     var gens = wg.WGbase.GetGenerators().ToArray();
    //     var pmap1 = Generators.Select(e => e.e).Zip(gens.Reverse()).ToDictionary(e => e.First, e => e.Second).ToDictionary();
    //     var iso1 = Group.Hom(PCGroup, Group.HomomorphismMap(PCGroup, wg, pmap1));
    //     var iso2 = Group.Hom(wg, Group.HomomorphismMap(wg, PCGroup, pmap1.ToDictionary(e => e.Value, e => e.Key)));
    //     var pcSeq = PolyCyclicSequence<Word>.Create(wg, gens.Reverse().ToArray());
    //     return (pcSeq, iso1, iso2);
    // }

    public List<List<(PolyCyclicSequence<T>, bool IsDirectProd)>> FullSplit()
    {
        var all = new List<List<(PolyCyclicSequence<T>, bool IsDirectProd)>>();
        var splits = Split().ToArray();
        foreach (var (lhs, rhs, isDirectProduct) in splits.OrderBy(e => e.N.PCGroup.Order))
        {
            foreach (var fullSplit in rhs.FullSplit())
                all.Add(fullSplit.Prepend((lhs, isDirectProduct)).ToList());
        }

        if (all.Count == 0)
            all.Add([(this, true)]);

        return all;
    }

    public string Relators()
    {
        var G = Normals.Last();
        var gens0 = Generators.SkipLast(1).Reverse().ToArray();
        var allNormalForms = gens0.Select(e => e.K.Order.SeqLazy().ToDictionary(i => i, _ => e.e).ToArray()).MultiLoop()
            .ToDictionary(l => G.OpSeq(l.Select(e => G.Times(e.Value, e.Key))), l => l.ToArray());
        var Xis = gens0.Length.SeqLazy().ToDictionary(i => i, i => $"{(char)('A' + i)}");

        var relators = new Dictionary<string, string>();
        foreach (var i in gens0.Length.SeqLazy())
        {
            var xi = gens0[i].e;
            var ord = G.ElementsOrders[xi];

            var Rii = ord.SeqLazy(1).Select(j => (j, xi_j: G.Times(xi, j)))
                .Select(e => (e.j, e.xi_j, seq: allNormalForms[e.xi_j]))
                .First(e => e.seq.Take(i + 1).All(f => f.Key == 0));
            relators[$"R{i}{i}"] = $"{Xis[i]}{Rii.j}{StrRels(Xis, Rii.seq)}";
        }

        foreach (var i in gens0.Length.SeqLazy())
        {
            var xi = gens0[i].e;
            foreach (var Rij in gens0.Length.SeqLazy().Where(j => j < i)
                         .Select(j => (j, xj1_xi_xj: G.OpSeq(G.Invert(gens0[j].e), xi, gens0[j].e)))
                         .Select(e => (e.j, e.xj1_xi_xj, seq: allNormalForms[e.xj1_xi_xj]))
                         .Where(e => e.seq.Take(e.j + 1).All(f => f.Key == 0)))
            {
                relators[$"R{i}{Rij.j}"] = $"{Xis[Rij.j]}-1{Xis[i]}{Xis[Rij.j]}{StrRels(Xis, Rij.seq)}";
            }
        }

        var strRelators = relators.Values.Distinct().Glue(", ");
        // Console.WriteLine();
        // Console.WriteLine("PCGroup Relators: {0}", strRelators);
        // Console.WriteLine(Xis.ToDictionary(e => gens0[e.Key].e, e => e.Value).GlueMap());
        // Console.WriteLine();
        return strRelators;
    }

    public static PolyCyclicSequence<T> Create(ConcreteGroup<T> G, (T e, int pow)[] gens)
    {
        var all = new HashSet<PolyCyclicSequence<T>>(10 * G.Order) { new(G, []) };
        while (all.Count != 0)
        {
            var tmp = all.OrderBy(pol => pol.Generators.Count(e => !e.Split)).ToArray();
            all.Clear();

            foreach (var partialDec in tmp)
            {
                foreach (var (K, e, split) in DecompStep(partialDec.Normals.Last(), gens[partialDec.Generators.Count]))
                {
                    var N = K.Neutral().H;
                    var dec = partialDec.Add(K, e, split);
                    if (N.Order == 1)
                        return dec;
                    else
                        all.Add(dec);
                }
            }
        }
        
        throw new();
    }

    private static string StrRels(Dictionary<int, string> Xis, KeyValuePair<int, T>[] seq)
    {
        var strRii = seq.Index()
            .Select(e => e.Item.Key == 0 ? "" : e.Item.Key == 1 ? Xis[e.Index] : $"{Xis[e.Index]}{e.Item.Key}")
            .Glue();
        return string.IsNullOrEmpty(strRii) ? "" : $" = {strRii}";
    }

    public static IEnumerable<XSet<T>> KerHom(ConcreteGroup<T> G, int k)
    {
        var Cn = new Cn(k);
        return GroupCraft.AllMorphismsWithPruning(G, Cn).Where(hom => !hom.IsNull)
            .Select(hom => (ker: hom.Kernel().ToXSet(), img: hom.Image().ToXSet()))
            .Where(e => e.img.Count == Cn.Order && e.ker.Count * e.img.Count == G.Order)
            .Select(e => e.ker)
            .Distinct()
            .OrderByDescending(e => e.Count);
    }

    public static IEnumerable<XSet<T>> KerHom(ConcreteGroup<T> G, T a, int relOrd)
    {
        var aPows = relOrd.SeqLazy().Select(i => G.Times(a, i)).ToArray();
        return KerHom(G, relOrd)
            .Where(e => G.SetEquals(e.Grid2D(aPows).Select(t => G.Op(t.t2, t.t1)).ToHashSet()))
            .Distinct();
    }

    public static IEnumerable<(ConcreteGroup<Coset<T>> K, T e, bool Split)> DecompStep(ConcreteGroup<T> G, (T a, int relOrd) f)
    {
        var(a, relOrd) = f;
        var seq = KerHom(G, a, relOrd)
            .Select(ker => Group.Generate("N", G, ker.ToArray()))
            .Select(N => (N, K: G.Over2(N), e: a))
            .Select(e => (e.K, e.e, inter: e.N.IntersectCount(Group.GenerateElements(G, e.e)) == 1));

        return seq;
    }

    public static IEnumerable<(ConcreteGroup<Coset<T>> K, T e, bool Split)> DecompStep(ConcreteGroup<T> G)
    {
        var seq = IntExt.PrimesDecomposition(G.Order).Distinct().Order()
            .SelectMany(k => KerHom(G, k).Select(ker => Group.Generate("N", G, ker.ToArray())))
            .Select(N => (N, K: G.Over2(N)))
            .Select(e => (e.N, e.K, e: e.K.GetGenerators().First().X))
            .Select(e => (e.K, e.e, inter: e.N.IntersectCount(Group.GenerateElements(G, e.e)) == 1))
            .OrderBy(e => !e.inter)
            .ThenBy(e => e.K.Order);

        if (G.GroupType == GroupType.AbelianGroup)
            return seq.Take(1);

        return seq;
    }

    public static IEnumerable<PolyCyclicSequence<T>> Decomp(ConcreteGroup<T> G)
    {
        var all = new HashSet<PolyCyclicSequence<T>>(10 * G.Order) { new(G, []) };
        do
        {
            var tmp = all.OrderBy(pol => pol.Generators.Count(e => !e.Split)).ToArray();
            all.Clear();

            foreach (var partialDec in tmp)
            {
                foreach (var (K, e, split) in DecompStep(partialDec.Normals.Last()))
                {
                    var N = K.Neutral().H;
                    var dec = partialDec.Add(K, e, split);
                    if (N.Order == 1)
                        yield return dec;
                    else
                        all.Add(dec);
                }
            }
        } while (all.Count != 0);
    }

    public static int Depth(Dictionary<T, (T e, int pow)[]> expVec, T e) => expVec[e].TakeWhile(a => a.pow == 0).Count();

    public static int Led(Dictionary<T, (T e, int pow)[]> expVec, T e) => expVec[e][Depth(expVec, e)].pow;

    public static Dictionary<T, T[]> AllZ(ConcreteGroup<T> G, Dictionary<T, (T e, int pow)[]> expVec, int[] R)
    {
        return expVec.ToDictionary(e => e.Key, e => e.Value.Select(r => G.Times(r.e, r.pow)).ToArray())
            .Where(Z => Z.Value.Index().Any(z => G.Neutral().Equals(z.Item) || Depth(expVec, z.Item) == z.Index))
            .ToDictionary();
    }

    public static T Sift(ConcreteGroup<T> G, Dictionary<T, (T e, int pow)[]> expVec, int[] R, T[] Z, T g)
    {
        var n = Z.Length;
        var h = g;
        var d = Depth(expVec, h);
        while (d < n && !Z[d].Equals(G.Neutral()))
        {
            var k = Z[d];
            var rd = R[d];
            var ledk = Led(expVec, k);
            var e = IntExt.AmodP(Led(expVec, h) * IntExt.InvModPbez(ledk, rd), rd);
            var h0 = G.Op(G.Times(k, -e), h);
            d = Depth(expVec, h0);
            h = h0;
        }

        return h;
    }

    private static T[] CanonicalInducedPolyCyclicSeq(ConcreteGroup<T> G, Dictionary<T, (T e, int pow)[]> expVec, int[] R, T[] Y)
    {
        var Y0 = Y.ToArray();
        var ds = Y0.Select(y => Depth(expVec, y)).ToArray();
        for (int i = Y0.Length - 1; i >= 0; i--)
        {
            var b0 = expVec[Y0[i]];
            var f = IntExt.InvModPbez(b0[ds[i]].pow, R[ds[i]]);
            if (IntExt.AmodP(f * b0[ds[i]].pow, R[ds[i]]) != 1)
                continue;
            Y0[i] = G.Times(Y0[i], f);
            for (int j = i + 1; j < Y0.Length; j++)
            {
                var b1 = expVec[Y0[i]];
                Y0[i] = G.Op(Y0[i], G.Times(Y0[j], -b1[ds[j]].pow));
            }
        }

        return Y0;
    }
}
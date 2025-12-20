using FastGoat.Commons;

public static class AbelianExt
{
    static AbelianExt()
    {
        seqEq = EqualityComparer<int[]>.Create(
            (a, b) => a is not null && b is not null && a.SequenceEqual(b),
            obj => obj.Length);
    }

    private static EqualityComparer<int[]> seqEq { get; }

    public static IEnumerable<int[]> AllAbTypes(int k)
    {
        if (k == 1)
            return [[1]];

        var dec = IntExt.PrimesDec(k);
        return dec.Select(e => IntExt.Partitions32[e.Value].Select(l => l.Select(i => e.Key.Pow(i)).ToArray()))
            .MultiLoop()
            .Select(l => l.SelectMany(i => i).OrderDescending().ToArray());
    }

    public static (int p, int r)[] AbToElemPows(int[] seq)
    {
        return seq.SelectMany(n =>
        {
            if (n < 1)
                throw new();

            if (n == 1)
                return new (int, int)[0];

            return IntExt.PrimesDec(n).Select(e => (e.Key, e.Value)).Order().ToArray();
        }).Order().ToArray();
    }

    public static int[] ElemsToCan(int[] seq)
    {
        if (seq.Length == 0 || (seq.Length == 1 && seq[0] == 1))
            return new[] { 1 };

        var set = seq.ToList();
        while (true)
        {
            var arr = set.OrderDescending().ToArray();
            var (e, l) = set.OrderDescending().Select(e => (e, arr.FirstOrDefault(c => IntExt.Gcd(e, c) == 1, -1)))
                .FirstOrDefault(e => e.Item2 != -1, (-1, -1));
            if (e == -1)
                break;

            set.Remove(e);
            set.Remove(l);
            set.Add(e * l);
        }

        return set.OrderDescending().ToArray();
    }

    public static int[] AbToElems(params int[] mods) => AbToElemPows(mods).Select(e => e.p.Pow(e.r)).ToArray();

    public static int[] AbToCan(params int[] mods) => ElemsToCan(AbToElems(mods));

    public static IEnumerable<int[]> Stair(int[] pows)
    {
        yield return [0];

        var powsDesc = pows.OrderDescending().ToArray();
        var temp1 = new List<int[]>();
        temp1.AddRange(powsDesc[0].SeqLazy(1).Select(i => new[] { i }).ToArray());
        for (int i = 0; i < powsDesc.Length; i++)
        {
            var temp2 = temp1.ToList();
            temp1.Clear();
            foreach (var current in temp2)
            {
                yield return current;
                if (current.Length == pows.Length)
                    continue;

                var max = int.Min(powsDesc[current.Length], current[i]);
                foreach (var j in max.Range(1))
                    temp1.Add(current.Append(j).ToArray());
            }
        }
    }

    public static IEnumerable<int[]> SubGroups(int[] g)
    {
        var g1 = AbToElemPows(g).GroupBy(e => e.p).ToDictionary(
            e => e.Key,
            e => e.Select(f => f.r).OrderDescending().ToArray());

        var gStairs = g1.ToDictionary(e => e.Key, e => Stair(e.Value).ToArray());
        return gStairs.Select(e => e.Value.Select(f => (e.Key, f)).ToArray()).MultiLoop()
            .Select(l =>
                l.Where(li => li.f.All(i => i != 0)).SelectMany(li => li.f.Select(i => li.Key.Pow(i))).ToArray())
            .Select(l => AbToCan(l));
    }

    public static IEnumerable<int[]> Morphs(int[] g, int[] n)
    {
        var g1Subs = SubGroups(g).ToHashSet(seqEq);
        g1Subs.IntersectWith(SubGroups(n));
        return g1Subs;
    }

    public static IEnumerable<int> EltOrds(int mod, int ord)
    {
        var gcd = IntExt.Gcd(mod, ord);
        var n = mod / gcd;
        return IntExt.Coprimes(gcd).Select(e => e * n % mod);
    }

    public static IEnumerable<int[]> AbEltsOrds(int[] mods, int ord)
    {
        return mods.Select(mod => EltOrds(mod, ord).ToArray()).MultiLoop().Select(l => l.ToArray());
    }

    public record DetailsQuotient(int p, int[] gr, int[] nr, int[] idx, int[] ker)
    {
        public override string ToString()
        {
            return $"p={p} g:[{gr.Glue(",")}] n:[{nr.Glue(",")}] ker:[{ker.Glue(",")}] idx:[{idx.Glue(",")}]";
        }

        public static int[] Kernel(DetailsQuotient[] quos) =>
            AbToCan(quos.SelectMany(q => q.ker.Select(i => q.p.Pow(i))).ToArray());
    }

    public record AbQuotient(
        int[] g,
        int[] gElts,
        Dictionary<int, int[]> gP,
        int[] n,
        int[] nElts,
        Dictionary<int, int[]> nP,
        DetailsQuotient[][] quos)
    {
        public void Display()
        {
            Console.WriteLine($"g:[{gElts.Glue(",")}]");
            Console.WriteLine($"n:[{nElts.Glue(",")}]");
            quos.Println(l => $"[{DetailsQuotient.Kernel(l).Glue(",")}] {l.Glue(", ")}", "Quotients");
        }
    }
    
    static void AbQuotientStep(int[] g, int[] n, List<int> s, List<(int[], int[])> all)
    {
        if (n.Length == 0)
        {
            all.Add((s.ToArray(), g.Where(e => e != 0).ToArray()));
            return;
        }

        var mx = n.Max();
        var ns = n.SkipLast(1).ToArray();
        foreach (var idx in g.Length.Range().Where(i => !s.Contains(i) && g[i] >= mx))
        {
            var gs = g.Select((e, i) => i == idx && !s.Contains(i) ? e - mx : e).ToArray();
            var s0 = s.Append(idx).ToList();
            AbQuotientStep(gs, ns, s0, all);
        }
    }

    public static AbQuotient Quotients(int[] g, int[] n)
    {
        var gElts = AbToElems(g);
        var nElts = AbToElems(n);

        var gP = AbToElemPows(g).GroupBy(e => e.p).ToDictionary(e => e.Key, e => e.Select(f => f.r).ToArray());
        var nP = AbToElemPows(n).GroupBy(e => e.p).ToDictionary(e => e.Key, e => e.Select(f => f.r).ToArray());
        var sols = new Dictionary<int, List<(int[], int[])>>();
        foreach (var k in gP.Keys.Where(k => !nP.ContainsKey(k)))
            nP[k] = [];

        if (nP.Any(e =>
                e.Value.OrderDescending().Zip(gP[e.Key].OrderDescending()).Any(f => f.First > f.Second)))
            return new(g, gElts, gP, n, nElts, nP, [[]]);

        foreach (var (p, ln) in nP)
        {
            var lg = gP[p];
            var seq = new List<(int[], int[])>();
            AbQuotientStep(lg, ln, [], seq);
            sols[p] = seq;
        }

        var quos = sols
            .Select(e => e.Value
                .Select(l => (p: e.Key, ni: l.Item1, ker: l.Item2, pr: l.Item2.Select(i => e.Key.Pow(i)).ToArray()))
                .ToArray())
            .MultiLoop()
            .Select(l => l.ToArray())
            .Select(k => k.Select(e => new AbelianExt.DetailsQuotient(e.p, gP[e.p], nP[e.p], e.ni, e.ker)).ToArray())
            .ToArray();

        return new(g, gElts, gP, n, nElts, nP, quos);
    }

}
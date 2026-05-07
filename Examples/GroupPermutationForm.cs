using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Perms;

namespace Examples;

public static class GroupPermutationForm
{
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

    static ConcreteGroup<Perm> AbelianPerm(params int[] seq)
    {
        var cycles = seq.Select(c => Cycles(c)).ToArray();
        var dims = cycles.Select(a => a.Sn.N).ToArray();
        var gens = cycles.Select((a, i) => Padding(dims.Take(i).Sum(), a, dims.Skip(i + 1).Sum())).ToArray();
        var sn = gens[0].Sn;
        return Group.Generate(seq.ToAbString(), sn, gens.ToArray());
    }

    static string MetaCyclicName(int m, int n, int r)
    {
        if (r == 1)
            return $"C{m} x C{n}";
        if (n == 2)
        {
            if (r == m - 1)
                return $"D{2 * m}";
            if (int.IsPow2(m))
                return r == m / 2 - 1 ? $"QD{2 * m}" : $"MM{2 * m}";
        }

        if (n == 4 && r == m - 1 && m % 2 == 1)
            return $"Dic{m}";

        if (IntExt.Gcd(m, n * (r - 1)) == 1)
            return $"F({m}x:{n}){r}";

        return $"M({m}x:{n}){r}";
    }

    public static (Perm a, Perm b) rUmToPerm(int m, int r)
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

    static IEnumerable<(Perm a, Perm b)> ExtendB(Perm a, Perm b, int n)
    {
        var N0 = a.Sn.N;
        foreach (var perms in CyclesSplit(n).AllCombinations().Where(e => e.Length > 0))
        {
            var table = IntExt.PermAndCyclesFromType(perms.SelectMany(e => e.PermType).Order().ToArray()).perm;
            var b0 = new Sn(table.Length).CreateElementTable(table);
            var N = b0.Table.Length + N0;
            var sn = new Sn(N);
            var a1 = sn.CreateElementTable(a.Table.Concat(b0.Table.Length.SeqLazy(N0)).ToArray());
            var b2 = sn.CreateElementTable(b.Table.Concat(b0.Table.Select(i => i + N0).ToArray()).ToArray());
            yield return (a1, b2);
        }
    }

    static ConcreteGroup<Perm> SearchPermGroupMetaCyclic(int m, int n, int r)
    {
        var (a, b) = rUmToPerm(m, r);
        var sn = a.Sn;
        var ar = a ^ r;

        var ord = Group.ElementIsOrder(sn, b, n);
        var c1 = (b ^ n).Equals(sn.Neutral());
        var c2 = (b * a * (b ^ -1)).Equals(ar);
        if (!ord && c1 && c2)
        {
            foreach (var (a0, b0) in ExtendB(a, b, n))
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
            var name = MetaCyclicName(m, n, r);
            var G = Group.Generate(name, sn, a, b);
            if (G.Count() == n * m)
                return G;
        }

        return Group.Generate("G", sn, sn.Neutral());
    }

    static ConcreteGroup<Perm> Dihedral(int n) => SearchPermGroupMetaCyclic(n, 2, n - 1);
    static ConcreteGroup<Perm> SemiDihedral(int n) => SearchPermGroupMetaCyclic(1 << (n - 1), 2, (1 << (n - 2)) - 1);
    static ConcreteGroup<Perm> ModularMax(int n) => SearchPermGroupMetaCyclic(1 << (n - 1), 2, (1 << (n - 2)) + 1);

    static (Perm a, Perm b, string name) QuaternionGens(int n)
    {
        if (!int.IsPow2(n))
            throw new();

        var k1 = n / 2;
        var k2 = n / 4;
        var sn = new Sn(n);
        var a = sn.Cycle(k1.Range(1)) * sn.Cycle(k1.Range(k1 + 1));
        var b = Group.OpSeq(sn, k2.SeqLazy().Select(i => sn.Cycle(i + 1, n - i, i + 1 + k2, n - i - k2)));
        return (a, b, $"Q{n}");
    }

    static (Sn sn, Perm[] gens, string name) DicyclicGens(int n)
    {
        var k = IntExt.PrimesDecomposition(n).Count(i => i == 2);
        if (k != 0)
        {
            var m = n / (1 << k);
            var (b, c, Qname) = QuaternionGens(1 << (k + 2));
            if (m == 1)
                return (b.Sn, [b, c], Qname);

            var (a0, b0) = rUmToPerm(m, m - 1);
            var sn = new Sn(b0.Sn.N + b.Sn.N);
            var a1 = sn.CreateElementTable(a0.Table.Concat(b.Sn.N.SeqLazy(a0.Sn.N)).ToArray());
            var b1 = sn.CreateElementTable(b0.Sn.Neutral().Table.Concat(b.Table.Select(e => e + a0.Sn.N)).ToArray());
            var c1 = sn.CreateElementTable(b0.Table.Concat(c.Table.Select(e => e + a0.Sn.N)).ToArray());

            var H1 = Group.GenerateElements(a1.Sn, a1, b1, c1);
            if (H1.Count == 4 * n)
                return (a1.Sn, [a1, b1, c1], $"C{m} x: {Qname}");

            throw new();
        }

        var (b2, c2) = rUmToPerm(n, n - 1);
        foreach (var (b3, c3) in ExtendB(b2, c2, 4))
        {
            var H1 = Group.GenerateElements(b3.Sn, b3, c3);
            if (H1.Count == 4 * n)
                return (b3.Sn, [b3, c3], $"Dic{n}");
        }

        throw new();
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
        for (int n = 3; n < 33; n++)
        {
            var Dic_m = FG.DicyclicGL2p(n);
            var (sn, gens, name) = DicyclicGens(n);
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
}
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.GenericGroup;
using FastGoat.UserGroup.Perms;

namespace Craft.Craft;

public static class UGCraft
{
    public static GroupAction<Perm, XSet<int>> Image => (g, x) => new(x.Select(i => g.Table[i]));

    public static Perm ConcatPerm(params Perm[] perms)
    {
        var Ns = perms.Select(p => p.Sn.N).ToArray();
        var sn = new Sn(Ns.Sum());
        var cum = Ns.Aggregate(new[] { 0 }, (acc, e) => acc.Append(acc.Last() + e).ToArray());
        var table = perms.Zip(cum).SelectMany(e => e.First.Table.Select(i => i + e.Second)).ToArray();
        return sn.CreateElementTable(table);
    }

    public static Perm PaddingRight(Perm perm, int pad) => pad == 0 ? perm : ConcatPerm(perm, (new Sn(pad)).Neutral());
    public static Perm PaddingLeft(Perm perm, int pad) => pad == 0 ? perm : ConcatPerm((new Sn(pad)).Neutral(), perm);

    public static Perm[] CyclesSplit(int m)
    {
        var dec = IntExt.PrimesDec(m).Select(e => e.Key.Pow(e.Value)).Order().ToArray();
        return dec.Select(e => new Sn(e).Cycle(e.Range(1))).ToArray();
    }

    public static Perm Cycles(int m)
    {
        if (m == 1)
            return new Sn(1).Neutral();

        return ConcatPerm(CyclesSplit(m));
    }

    public static IEnumerable<int[]> CycleWalk(int[] c)
    {
        var n = c.Length;
        return n.SeqLazy(1).Select(i => n.SeqLazy(i).Select(j => c[j % n]).ToArray());
    }

    public static IEnumerable<int[][]> BlocksPermutes(int[][] cycles)
    {
        return cycles.GroupBy(c => c.Length).OrderBy(e => e.Key)
            .Select(a => (ord: a.Key, block: a.ToArray()))
            .Select(a => new Sn(a.block.Length).Select(p => (a.ord, block: p.Apply(a.block))))
            .ToArray()
            .MultiLoop()
            .Select(l => l.SelectMany(m => m.block).ToArray());
    }

    public static IEnumerable<Perm> InnerAut(Perm a, Perm b)
    {
        var sn = a.Sn;
        if (sn.N != b.Sn.N || !a.PermType.SequenceEqual(b.PermType))
            throw new($"a:{sn} a:{a.PermTypeStr} b:{b.Sn} b:{b.PermTypeStr}");

        var aOrbits = a.Orbits.OrderBy(e => e.Length).ToArray();
        var bOrbits = b.Orbits.OrderBy(e => e.Length).ToArray();
        foreach (var blocks in BlocksPermutes(bOrbits))
        {
            foreach (var cycle in blocks.Select(c => CycleWalk(c)).MultiLoop().Select(l => l.ToArray()))
            {
                var autTable = new int[sn.N];
                foreach (var (bOrb, aOrb) in cycle.Zip(aOrbits))
                    for (int i = 0; i < bOrb.Length; i++)
                        autTable[bOrb[i]] = aOrb[i];

                yield return sn.CreateElementTable(autTable);
            }
        }
    }

    public static IEnumerable<Perm> InnerAutCn(Perm a)
    {
        return IntExt.Coprimes(a.Order).Select(i => a ^ i).SelectMany(ar => InnerAut(a, ar));
    }

    public static void TestInnerAutCn(Perm a)
    {
        var sn = a.Sn;
        var allConjs1 = InnerAutCn(a).ToArray();
        var act = Group.ByConjugate(sn);
        var H = Group.Cycle(sn, a);
        var allConjs2 = sn.Where(g => H.Keys.Contains(act(g, a))).ToHashSet();

        Console.WriteLine($"a:{a} conj:{allConjs1.Length}");
        foreach (var x in allConjs1.OrderBy(x => H[act(x, a)]).ThenBy(x => x))
            Console.WriteLine($"    x:{x,-40} x*a*x^-1 = a^{H[act(x, a)],-2} = {act(x, a)}");

        Console.WriteLine($"a:{a} conj:{allConjs2.Count}");
        foreach (var x in allConjs2.OrderBy(x => H[act(x, a)]).ThenBy(x => x))
            Console.WriteLine($"    x:{x,-40} x*a*x^-1 = a^{H[act(x, a)],-2} = {act(x, a)}");

        var check = allConjs2.Count == allConjs1.Length && allConjs2.SetEquals(allConjs1);
        Console.WriteLine($"a:{a} conj:{allConjs1.Length} Set Equal:{check}");
        Console.WriteLine();
        if (!check)
        {
            allConjs2.Except(allConjs1).Println("Missing");
            throw new();
        }
    }

    public static void TestInnerAutCnLazy(Perm a)
    {
        var sn = a.Sn;
        var allConjs = InnerAutCn(a);
        var act = Group.ByConjugate(sn);
        var H = Group.Cycle(sn, a);

        var nb = 0;
        Console.WriteLine($"a:{a}");
        foreach (var x in allConjs)
        {
            ++nb;
            Console.WriteLine($"    x:{x,-40} x*a*x^-1 = a^{H[act(x, a)],-2} = {act(x, a)}");
        }

        Console.WriteLine($"sn:{sn} a:{a} total:{nb}");
        Console.WriteLine();
    }

    public static void TestInnerAut(Perm a, Perm b)
    {
        var sn = a.Sn;
        var allConjs1 = InnerAut(a, b).ToArray();
        var act = Group.ByConjugate(sn);
        var allConjs2 = sn.Where(g => act(g, a).Equals(b)).ToHashSet();
        var digits = -(sn.N + 1).SeqLazy().Sum(e => $"{e},  ".Length);
        var fmt = $"    x:{{0,{digits}}} = {{1,{digits}}} x*a*x^-1 = b = {{2,{digits}}} {{3}}";

        Console.WriteLine($"a:{a} b:{b} conj:{allConjs1.Length}");
        foreach (var x in allConjs1.OrderBy(x => x.Order).ThenBy(x => x.Orbits.Length).ThenBy(x => x))
            Console.WriteLine(fmt, x, $"[{x.Table.Glue(", ")}]", act(x, a), act(x, a).Equals(b));

        Console.WriteLine($"a:{a} b:{b} conj:{allConjs2.Count}");
        foreach (var x in allConjs2.OrderBy(x => x.Order).ThenBy(x => x.Orbits.Length).ThenBy(x => x))
            Console.WriteLine(fmt, x, $"[{x.Table.Glue(", ")}]", act(x, a), act(x, a).Equals(b));

        var check = allConjs2.Count == allConjs1.Length && allConjs2.SetEquals(allConjs1);
        Console.WriteLine($"a:{a} b:{b} conj:{allConjs1.Length} Set Equal:{check}");
        Console.WriteLine();
        if (!check)
            throw new();
    }

    public static void FirstExamples()
    {
        TestInnerAutCn(Cycles(4));
        TestInnerAutCn(Cycles(5));
        TestInnerAutCn(Cycles(6));
        TestInnerAutCn(Cycles(10));
        TestInnerAutCn(ConcatPerm(Cycles(3), Cycles(4)));

        {
            var a = ConcatPerm(Cycles(3), Cycles(2));
            var b = ConcatPerm(Cycles(2), Cycles(3));
            TestInnerAut(a, b);
        }
        {
            var a = ConcatPerm(Cycles(3), Cycles(3));
            var b = ConcatPerm(Cycles(3), Cycles(3));
            TestInnerAut(a, b);
        }
        {
            var a = ConcatPerm(Cycles(3), Cycles(5));
            var b = ConcatPerm(Cycles(5), Cycles(3));
            TestInnerAut(a, b);
        }
    }

    public static void SchreierSims(ConcreteGroup<Perm> G)
    {
        DisplayGroup.HeadElements(G);
        var sn = G.Neutral().Sn;
        var n = sn.N;
        var listGi = n.SeqLazy()
            .Select(i => Group.Generate($"G{i}", sn, G.Where(e => i.SeqLazy().All(j => e.Table[j] == j)).ToArray()))
            .ToArray();
        Console.WriteLine(listGi.Select(g => g.ShortName).Glue(", "));
        var listSi = (n - 1).SeqLazy(1)
            .Select(i => Group.Cosets(listGi[i - 1], listGi[i], CosetType.Left).Values.ToHashSet());
        var gens = listSi.Select(s => s.Select(c => c.X).ToXSet()).ToArray();
        Console.WriteLine(gens.Glue(", "));
        Console.WriteLine();

        Console.WriteLine(SchreierSimsFast(G).Glue(", "));
        Console.WriteLine();
    }

    public static XSet<Perm>[] SchreierSimsFast(ConcreteGroup<Perm> G)
    {
        var sn = G.Neutral().Sn;
        var n = sn.N;
        var listGi = n.SeqLazy().Select(i => G.Where(e => i.SeqLazy().All(j => e.Table[j] == j)).ToArray()).ToArray();
        return listGi.Select((Gi, idx) =>
        {
            XSet<int> i = new(idx);
            var Oi = Gi.Select(g => Image(g, i)).ToHashSet();
            return Oi.Select(j => Gi.OrderBy(e => G.ElementsOrders[e]).ThenBy(e => e)
                    .First(g => Image(g, j).Equals(i)))
                .Where(e => e.Order != 1).ToXSet();
        }).ToArray();
    }

    public static void TestSchreierSims()
    {
        var n = 4;
        var sn = new Sn(n);
        SchreierSims(Group.Generate("Id", sn, sn.Neutral()));
        SchreierSims(Group.Generate("C4", sn, sn[(1, 2, 3, 4)]));
        SchreierSims(Group.Generate("A4", sn, sn[(1, 2, 3)], sn[(2, 3, 4)]));
        SchreierSims(Group.Generate("S4", sn, sn.GetGenerators().ToArray()));
    }
}
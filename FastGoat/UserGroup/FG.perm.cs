using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.GenericGroup;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Perms;

namespace FastGoat.UserGroup;

public static partial class FG
{
    public static Perm ConcatPerm(params Perm[] perms)
    {
        var perms0 = perms.Where(e => e.Sn.N != 1).ToArray();
        if (perms0.Length == 0)
            return new Sn(1).Neutral();

        var Ns = perms0.Select(p => p.Sn.N).ToArray();
        var sn = new Sn(Ns.Sum());
        var cum = Ns.Aggregate(new[] { 0 }, (acc, e) => acc.Append(acc.Last() + e).ToArray());
        var table = perms0.Zip(cum).SelectMany(e => e.First.Table.Select(i => i + e.Second)).ToArray();
        return sn.CreateElementTable(table);
    }

    public static Perm PaddingRight(Perm perm, int pad) => pad == 0 ? perm : ConcatPerm(perm, (new Sn(pad)).Neutral());
    public static Perm PaddingLeft(Perm perm, int pad) => pad == 0 ? perm : ConcatPerm((new Sn(pad)).Neutral(), perm);
    static Perm Padding(int padLeft, Perm perm, int padRight) => PaddingLeft(PaddingRight(perm, padRight), padLeft);

    static Perm[] CyclesSplit(int m)
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

    static (string name, Sn sn, Perm[] gens) AbelianPermGens(params int[] seq)
    {
        var cycles = seq.Select(c => Cycles(c)).ToArray();
        var dims = cycles.Select(a => a.Sn.N).ToArray();
        var name = cycles.Select(c => c.Order).ToAbString();
        var gens = cycles.Select((a, i) => Padding(dims.Take(i).Sum(), a, dims.Skip(i + 1).Sum())).ToArray();
        return (name, gens[0].Sn, gens);
    }

    static (Perm a, Perm b) rUmToPerm(int m, int r)
    {
        var cycles = CyclesSplit(m).Select(c =>
        {
            var N = c.Sn.N;
            if (IntExt.Gcd(N, r) != 1)
                return c;

            var ri = new ZnInt(N, r).Inv();
            return c.Sn.CreateElementTable(N.Range().Select(i => (i * ri).K).ToArray());
        }).ToArray();

        return (Cycles(m), ConcatPerm(cycles.ToArray()));
    }

    static (Perm a, Perm b) OpsCnByUm(int n, Perm a, Perm b)
    {
        var n0 = IntExt.PrimesDec(n)
            .Select(p => p.Key.Pow(p.Value))
            .Where(e => b.Order % e != 0)
            .Aggregate(1, (acc, e) => e * acc);
        var b0 = Cycles(n0);
        return (PaddingRight(a, b0.Sn.N), ConcatPerm(b, b0));
    }

    public static (string name, Sn sn, Perm[] gens) MetaCyclicGens(int m, int n, int r)
    {
        var (a, b) = rUmToPerm(m, r);
        (a, b) = OpsCnByUm(n, a, b);
        var sn = a.Sn;
        var name = MetaCyclicName(m, n, r);
        return (name, sn, [a, b]);
    }

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
        var (b3, c3) = OpsCnByUm(4, b2, c2);
        return ($"Dic{n}", b3.Sn, [b3, c3]);
    }

    public static ConcreteGroup<Perm> PermGroup(string name, int n, params ValueType[] generators)
    {
        var sn = new Sn(n);
        var gi = generators.Select(g => sn.ComposesCycles(Tuple2Array.ComplexTuples(g))).ToArray();
        return Group.Generate(name, sn, gi);
    }

    public static ConcreteGroup<Perm> PermGroup(int n, params ValueType[] generators) => PermGroup("G", n, generators);

    public static ConcreteGroup<Perm> Symmetric(int n) => new Symm(n);

    public static ConcreteGroup<Perm> Alternate(int n)
    {
        if (n < 3)
            throw new GroupException(GroupExceptionType.GroupDef);

        var sn = new Sn(n);
        var gi = (n - 2).Range(3).Select(i => sn[(1, 2, i)]).ToArray();
        return Group.Generate($"Alt{n}", sn, gi);
    }

    public static ConcreteGroup<Perm> AbelianPerm(params int[] seq)
    {
        var (name, sn, gens) = AbelianPermGens(seq);
        return Group.Generate(name, sn, gens.ToArray());
    }

    public static ConcreteGroup<Perm> MetaCyclicPg(int m, int n, int r)
    {
        var (name, sn, gens) = MetaCyclicGens(m, n, r);
        return Group.Generate(name, sn, gens);
    }

    public static List<ConcreteGroup<Perm>> FrobeniusPg(int order)
    {
        return IntExt.Dividors(order).Where(d => d > 1 && d % 2 == 1)
            .SelectMany(m => FrobeniusGetR(m, order / m).Select(r => (m, n: order / m, r)))
            .Select(e => MetaCyclicPg(e.m, e.n, e.r))
            .ToList();
    }

    public static ConcreteGroup<Perm> Dihedral(int n) => MetaCyclicPg(n, 2, n - 1);
    public static ConcreteGroup<Perm> SemiDihedralPg(int n) => MetaCyclicPg(1 << (n - 1), 2, (1 << (n - 2)) - 1);
    public static ConcreteGroup<Perm> ModularMaxPg(int n) => MetaCyclicPg(1 << (n - 1), 2, (1 << (n - 2)) + 1);

    public static ConcreteGroup<Perm> QuaternionPg(int m)
    {
        var (name, sn, gens) = QuaternionGens(m);
        return Group.Generate(name, sn, gens);
    }

    public static ConcreteGroup<Perm> DicyclicPg(int m)
    {
        var (name, sn, gens) = DicyclicGens(m);
        return Group.Generate(name, sn, gens);
    }

    public static List<ConcreteGroup<Perm>> MetaCyclicPg(int ord)
    {
        return IntExt.Dividors(ord).Where(d => d > 1)
            .SelectMany(m => MetaCyclicSdpGetR(m, ord / m).Select(r => (m, n: ord / m, r)))
            .DistinctBy(e => e.r == 1 ? (e.m * e.n, 1, 1) : e)
            .Select(e => MetaCyclicPg(e.m, e.n, e.r))
            .ToList();
    }
    
    private static GroupAction<T, Perm> Action2Perm<T>(ConcreteGroup<T> g, Sn sn) where T : struct, IElt<T>
    {
        var mapEltIdx = g.Select((e, i) => (e, i: i + 1))
            .OrderBy(e => g.ElementsOrders[e.e])
            .ToDictionary(e => e.e, e => e.i);
        var mapIdxElt = mapEltIdx.ToDictionary(e => e.Value, e => e.Key);
        return (e, p) => sn.CreateElement(p.Table.Select(i => mapEltIdx[g.Op(e, mapIdxElt[i + 1])]).ToArray());
    }

    public static (ConcreteGroup<Perm>, Dictionary<T, Perm> mapReg) ToPermGroup<T>(this ConcreteGroup<T> g)
        where T : struct, IElt<T>
    {
        var sn = new Sn(g.Count());
        var act = Action2Perm(g, sn);
        var mapGens = g.GetGenerators().ToDictionary(e => e, e => act(e, sn.Neutral()));
        var mapReg = g.ToDictionary(e => e, e => act(e, sn.Neutral()));
        return (Group.Generate($"{g.Name}", sn, mapGens.Values.ToArray()), mapReg);
    }
    
    public static (ConcreteGroup<Perm>, Dictionary<Perm, int> idx) RegPermAutGroup(ConcreteGroup<Automorphism<Perm>> aut)
    {
        var Dom = aut.Neutral().Domain;
        var sn = Dom.Neutral().Sn;
        if (sn.N != Dom.Count())
            throw new($"{sn} G:{Dom.Count()}");

        var idx = Dom.Index().ToDictionary(e => e.Item, e => e.Index);
        var gens = aut.GetGenerators().Select(e => e.AutMap.ToDictionary(f => idx[f.Key], f => idx[f.Value]))
            .Select(e => sn.CreateElementTable(e.OrderBy(f => f.Key).Select(f => f.Value).ToArray()))
            .ToArray();

        return (Group.Generate(aut.Name, sn, gens), idx);
    }

}
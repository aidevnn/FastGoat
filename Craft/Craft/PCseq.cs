using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.Subgroups;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Integers;

namespace Craft.Craft;

public readonly partial struct PCseq<T> : IElt<PCseq<T>> where T : struct, IElt<T>
{
    public ConcreteGroup<T> PCGroup { get; }
    public T[] X { get; }
    public int[] R { get; }
    public (T xi, int ei)[] XR => X.Zip(R).ToArray();

    public PCseq(ConcreteGroup<T> G, T[] x, int[] r)
    {
        (PCGroup, X, R) = (G, x, r);
        Hash = (PCGroup.Order, R.Length).GetHashCode();
    }

    public bool Equals(PCseq<T> other) =>
        R.Length == other.R.Length
        && R.SequenceEqual(other.R)
        && PCGroup.SetEquals(other.PCGroup);

    public int CompareTo(PCseq<T> other)
    {
        var compLength = R.Length.CompareTo(other.R.Length);
        if (compLength != 0)
            return compLength;

        return R.SequenceCompareTo(other.R);
    }

    public int Hash { get; }

    public override int GetHashCode() => Hash;

    public override string ToString() => $"{X.Zip(R).ToXArray()}";

    public Dictionary<T, (T e, int pow)[]> ExponentVectors()
    {
        var G = PCGroup;
        var x = X.ToArray();
        var expVec = R.Select((k, i) => k.Range().Select(ei => (i, ei)).ToArray())
            .MultiLoop()
            .Select(l => l.Select(e => (xi: x[e.i], e.ei)).ToArray())
            .ToDictionary(l => G.OpSeq(l.Select(e => G.Times(e.xi, e.ei))), l => l);

        if (!G.SetEquals(expVec.Keys))
        {
            expVec.Println(e => $"{e.Key} = {e.Value.ToXArray()}");
            throw new("ExpVec problem");
        }

        return expVec;
    }

    public PCseq<T> InducedPCseq(ConcreteGroup<T> H) => InducedPCseq(H, H.GetGenerators().ToArray());

    public PCseq<T> InducedPCseq(ConcreteGroup<T> H, T[] gensH)
    {
        var expVec = ExponentVectors();
        var _R = R;
        var G = PCGroup;
        var Z = R.Select(_ => G.Neutral()).ToArray();
        var n = Z.Length;
        var gensH0 = new HashSet<T>(gensH);
        GlobalStopWatch.InfiniteLoopBreakerReset();
        while (gensH0.Any())
        {
            GlobalStopWatch.InfiniteLoopBreaker(50);
            var g = gensH0.First();
            gensH0.Remove(g);
            var h = Sift(G, expVec, _R, Z, g);
            var d = Depth(expVec, h).depth;
            if (d < n)
            {
                gensH0.Add(G.Times(h, _R[d]));
                gensH0.UnionWith(Z.Select(z => G.Commute(h, z)));
                Z[d] = h;
            }

            gensH0.Remove(G.Neutral());
        }

        var Y = CanonicalInducedPCSeq(G, expVec, _R, Z.Where(z => !z.Equals(G.Neutral())).ToArray());
        return new(H, Y, Y.Select(yi => _R[Depth(expVec, yi).depth]).ToArray());
    }

    public (ConcreteGroup<T>[] Normals, ConcreteGroup<Coset<T>>[] Factors) NormalGroup()
    {
        var G = PCGroup;
        var _X = X.ToArray();
        var normals = (X.Length + 1).SeqLazy()
            .Select(i => Group.Generate("N", G, _X.Skip(i).Append(G.Neutral()).ToArray()))
            .ToArray();
        var factors = X.Length.SeqLazy().Select(i => normals[i].Over2(normals[i + 1])).ToArray();
        return (normals, factors);
    }
    
    public string Relators()
    {
        var G = PCGroup;
        var _X = X.ToArray();
        var Xis = X.Length.SeqLazy().ToDictionary(i => i, i => $"{(char)('A' + i)}");
        var expVec = ExponentVectors();
        var relators = new Dictionary<string, string>();
        foreach (var i in X.Length.SeqLazy())
        {
            var xi = X[i];
            var ord = G.ElementsOrders[xi];

            var Rii = ord.SeqLazy(1).Select(j => (j, xi_j: G.Times(xi, j)))
                .Select(e => (e.j, e.xi_j, seq: expVec[e.xi_j]))
                .First(e => e.seq.Take(i + 1).All(f => f.pow == 0));
            relators[$"R{i}{i}"] = $"{Xis[i]}{Rii.j}{StrRels(Xis, Rii.seq)}";
        }

        foreach (var i in X.Length.SeqLazy())
        {
            var xi = X[i];
            foreach (var Rij in X.Length.SeqLazy().Where(j => j < i)
                         .Select(j => (j, xj1_xi_xj: G.OpSeq(G.Invert(_X[j]), xi, _X[j])))
                         .Select(e => (e.j, e.xj1_xi_xj, seq: expVec[e.xj1_xi_xj]))
                         .Where(e => e.seq.Take(e.j + 1).All(f => f.pow == 0)))
            {
                relators[$"R{i}{Rij.j}"] = $"{Xis[Rij.j]}-1{Xis[i]}{Xis[Rij.j]}{StrRels(Xis, Rij.seq)}";
            }
        }

        var strRelators = relators.Values.Distinct().Glue(", ");
        Console.WriteLine("PCGroup Relators: {0}", strRelators);
        Console.WriteLine(Xis.ToDictionary(e => _X[e.Key], e => e.Value).GlueMap());
        Console.WriteLine();
        return strRelators;
    }
}
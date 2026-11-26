using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.Subgroups;
using FastGoat.UserGroup.Integers;

namespace Craft.Craft;

public readonly partial struct PCseq<T> : IElt<PCseq<T>> where T : struct, IElt<T>
{
    private static string StrRels(Dictionary<int, string> Xis, (T xi, int ei)[] seq)
    {
        var strRii = seq.Index()
            .Select(e => e.Item.ei == 0 ? "" : e.Item.ei == 1 ? Xis[e.Index] : $"{Xis[e.Index]}{e.Item.ei}")
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
            .Distinct();
    }

    public static IEnumerable<(ConcreteGroup<T> N, T ei, int ri)> DecompStep(ConcreteGroup<T> G)
    {
        var seq = IntExt.PrimesDecomposition(G.Order).Distinct().Order()
            .SelectMany(k => KerHom(G, k).Select(ker => Group.Generate("N", G, ker.ToArray())))
            .Select(N => (N, K: G.Over2(N)))
            .Select(e => (e.N, e: e.K.GetGenerators().First().X, e.K.Order));

        return seq;
    }

    static IEnumerable<(T[] X, int[] R)> DecompInternal(ConcreteGroup<T> G)
    {
        if (G.Order == 1)
            yield return ([], []);
        else
        {
            foreach (var (N, xi, ri) in DecompStep(G))
            {
                foreach (var subSeq in Decomp(N).Distinct())
                    yield return (subSeq.X.Prepend(xi).ToArray(), subSeq.R.Prepend(ri).ToArray());
            }
        }
    }
    
    public static IEnumerable<PCseq<T>> Decomp(ConcreteGroup<T> G) => DecompInternal(G).Select(e => new PCseq<T>(G, e.X, e.R));

    public static (int depth, int lt) Depth(Dictionary<T, (T xi, int ei)[]> expVec, T e)
    {
        var nf = expVec[e];
        var depth = nf.TakeWhile(a => a.ei == 0).Count();
        if (depth == nf.Length)
            return (depth, 0);
        
        var lt = nf[depth].ei;
        return (depth, lt);
    }

    public static T Sift(ConcreteGroup<T> G, Dictionary<T, (T e, int pow)[]> expVec, int[] R, T[] Z, T g)
    {
        var n = Z.Length;
        var h = g;
        var d = Depth(expVec, h).depth;
        while (d < n && !Z[d].Equals(G.Neutral()))
        {
            var k = Z[d];
            var rd = R[d];
            var ledk = Depth(expVec, k).lt;
            var e = IntExt.AmodP(Depth(expVec, h).lt * IntExt.InvModPbez(ledk, rd), rd);
            var h0 = G.Op(G.Times(k, -e), h);
            d = Depth(expVec, h0).depth;
            h = h0;
        }

        return h;
    }

    private static T[] CanonicalInducedPCSeq(ConcreteGroup<T> G, Dictionary<T, (T e, int pow)[]> expVec, int[] R, T[] Y)
    {
        var Y0 = Y.ToArray();
        var ds = Y0.Select(y => Depth(expVec, y).depth).ToArray();
        for (int i = Y0.Length - 1; i >= 0; i--)
        {
            var b0 = expVec[Y0[i]];
            var f = IntExt.InvModPbez(b0[ds[i]].pow, R[ds[i]]);
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
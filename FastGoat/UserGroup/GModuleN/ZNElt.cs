using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Integers;
using QPoly =
    FastGoat.Structures.VecSpace.Polynomial<FastGoat.Structures.VecSpace.EPoly<FastGoat.UserGroup.Integers.ZnInt>,
        FastGoat.Structures.VecSpace.Xi>;

namespace FastGoat.UserGroup.GModuleN;

public readonly struct ZNElt<Tn, Tg> : IElt<ZNElt<Tn, Tg>>
    where Tn : struct, IElt<Tn>
    where Tg : struct, IElt<Tg>
{
    public AbelianDirectSum<Tn> N { get; }
    public Dictionary<Tn, QPoly> Coefs { get; }
    public Indeterminates<Xi> Indeterminates { get; }
    public MapElt<Tg, Automorphism<Tn>> L { get; }

    public ZNElt(EPoly<ZnInt> scalar, Indeterminates<Xi> ind, AbelianDirectSum<Tn> N0, MapElt<Tg, Automorphism<Tn>> L0)
    {
        (N, L) = (N0, L0);
        Indeterminates = ind;
        KZero = scalar.Zero;
        var zero = Zero;
        Coefs = N.Decomp.ToDictionary(e => e.g, _ => zero);
        Hash = (N, ind).GetHashCode();
    }

    public ZNElt(EPoly<ZnInt> scalar, Indeterminates<Xi> ind, AbelianDirectSum<Tn> N0, MapElt<Tg, Automorphism<Tn>> L0, Xi[] xs)
    {
        (N, L) = (N0, L0);
        Indeterminates = ind;
        KZero = scalar.Zero;
        var zero = Zero;
        if (xs.Length != N.Decomp.Count)
            throw new();

        Coefs = N.Decomp.Zip(xs.Select(x => zero.X(x))).ToDictionary(e => e.First.g, e => e.Second);
        Hash = (N, ind).GetHashCode();
    }


    public ZNElt(EPoly<ZnInt> scalar, Indeterminates<Xi> ind, AbelianDirectSum<Tn> N0, MapElt<Tg, Automorphism<Tn>> L0, Dictionary<Tn, QPoly> map)
    {
        (N, L) = (N0, L0);
        Coefs = new(map);
        Indeterminates = ind;
        KZero = scalar.Zero;
        Hash = (N, ind).GetHashCode();
    }
    public EPoly<ZnInt> KZero { get; }

    public QPoly Zero => new(Indeterminates, KZero);

    public int Hash { get; }
    public QPoly this[Tn n] => Coefs.TryGetValue(n, out var e) ? e : Zero;

    public bool Equals(ZNElt<Tn, Tg> other)
    {
        var elt = this;
        return N.Decomp.All(n => other[n.g].Equals(elt[n.g]));
    }

    public int CompareTo(ZNElt<Tn, Tg> other)
    {
        var elt = this;
        return N.Decomp.Select(n => elt[n.g].CompareTo(other[n.g])).FirstOrDefault(k => k != 0, 0);
    }

    public ZNElt<Tn, Tg> Add(ZNElt<Tn, Tg> n)
    {
        var elt = this;
        var decomp = N.Decomp;
        var map = decomp.Select(e => (e, (elt[e.g] + n[e.g])))
            .Where(e => !e.Item2.IsZero()).ToDictionary(e => e.e.g, e => e.Item2);
        return new(KZero, Indeterminates, N, L, map);
    }

    public ZNElt<Tn, Tg> Sub(ZNElt<Tn, Tg> n)
    {
        var elt = this;
        var decomp = N.Decomp;
        var map = decomp.Select(e => (e, (elt[e.g] - n[e.g])))
            .Where(e => !e.Item2.IsZero()).ToDictionary(e => e.e.g, e => e.Item2);
        return new(KZero, Indeterminates, N, L, map);
    }

    public ZNElt<Tn, Tg> Opp()
    {
        var elt = this;
        var decomp = N.Decomp;
        var map = decomp.Select(e => (e, -elt[e.g]))
            .Where(e => !e.Item2.IsZero()).ToDictionary(e => e.e.g, e => e.Item2);
        return new(KZero, Indeterminates, N, L, map);
    }

    public ZNElt<Tn, Tg> Act(int k)
    {
        var elt = this;
        var decomp = N.Decomp;
        var map = decomp.Select(e => (e, k * elt[e.g]))
            .Where(e => !e.Item2.IsZero()).ToDictionary(e => e.e.g, e => e.Item2);
        return new(KZero, Indeterminates, N, L, map);
    }

    public ZNElt<Tn, Tg> Act(Tg g)
    {
        var N0 = N;
        var L0 = L;
        var coefs = Coefs;
        var decomp = N.Decomp.Select(e => e.g).ToArray();
        var map0 = N.Decomp.Select(e =>
                (e, N0.GEltToCan(L0[g][e.g]).Ei.Select((ei, k) => ei.K * coefs[decomp[k]]).Aggregate((a0, a1) => a0 + a1)))
            .ToDictionary(e => e.e.g, e => e.Item2);

        return new(KZero, Indeterminates, N, L, map0);
    }

    public bool IsZero() => Coefs.Values.All(e => e.IsZero() || e.LeadingDetails.lc.IsZero());
    public bool IsKnown() => Coefs.Values.All(e => e.IsZero() || e.LeadingDetails.lm.Degree == 0);

    public override int GetHashCode() => Hash;

    public override string ToString()
    {
        if (IsZero())
            return "0";

        var coefs = Coefs.Where(e => !e.Value.IsZero()).Select(e => $"{e.Key}^({e.Value})");
        return coefs.Glue(" + ");
    }

    public static ZNElt<Tn, Tg> operator +(ZNElt<Tn, Tg> a, ZNElt<Tn, Tg> b) => a.Add(b);

    public static ZNElt<Tn, Tg> operator -(ZNElt<Tn, Tg> a, ZNElt<Tn, Tg> b) => a.Sub(b);

    public static ZNElt<Tn, Tg> operator -(ZNElt<Tn, Tg> a) => a.Opp();

    public static ZNElt<Tn, Tg> operator *(int k, ZNElt<Tn, Tg> b) => b.Act(k);

    public static ZNElt<Tn, Tg> operator *(Tg g, ZNElt<Tn, Tg> b) => b.Act(g);
}
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.VecSpace;

namespace FastGoat.UserGroup.GModuleM;

public readonly struct GZNElt<Tn, Tg> : IGmoduleMElt<Tn, Tg, GZNElt<Tn, Tg>>, IElt<GZNElt<Tn, Tg>>
    where Tg : struct, IElt<Tg>
    where Tn : struct, IElt<Tn>
{
    public ConcreteGroup<Tn> N { get; }
    public ConcreteGroup<Tg> G { get; }
    public SortedList<Tg, ZNElt<Tn>> Coefs { get; }

    public GZNElt(ConcreteGroup<Tn> n, ConcreteGroup<Tg> g, params Xi[] unknowns)
    {
        if (n.GroupType == GroupType.NonAbelianGroup)
            throw new GroupException(GroupExceptionType.OnlyAbelianGroups);

        (N, G) = (n, g);
        var z = new ZNElt<Tn>(N, unknowns);
        Coefs = new(G.ToDictionary(e => e, _ => z.Zero));
        Hash = (N, G).GetHashCode();
    }

    public GZNElt(ConcreteGroup<Tn> n, ConcreteGroup<Tg> g, IDictionary<Tg, ZNElt<Tn>> map)
    {
        (N, G) = (n, g);
        Coefs = new(map);
        Hash = (N, G).GetHashCode();
    }

    public bool Equals(GZNElt<Tn, Tg> other) => Hash == other.Hash && Coefs.All(e => e.Value.Equals(other[e.Key]));

    public int CompareTo(GZNElt<Tn, Tg> other)
    {
        var compUnks = NbUnknownsStrict().CompareTo(other.NbUnknownsStrict());
        if (compUnks != 0)
            return compUnks;

        return Coefs.Values.SequenceCompareTo(other.Coefs.Values);
    }

    public int Hash { get; }

    public ZNElt<Tn> this[Tg g] => Coefs[g];
    public bool IsZero() => Coefs.Values.All(z => z.IsZero());
    public bool IsKnown() => Coefs.Values.All(z => z.IsKnown());
    public int NbUnknowns() => Coefs.Values.Sum(z => z.NbUnknowns());

    public Xi[] GetUnknowns()
    {
        return Coefs.Values.SelectMany(e => e.Coefs.Where(c => !c.Value.IsZero()).Select(c => c.Key)).Distinct().Order().ToArray();
    }

    public int NbUnknownsStrict() => GetUnknowns().Length;

    public GZNElt<Tn, Tg> Zero => new(N, G, Coefs.ToDictionary(e => e.Key, e => e.Value.Zero));
    public ZNElt<Tn> ZnEltZero => this[G.Neutral()].Zero;

    public GZNElt<Tn, Tg> Add(GZNElt<Tn, Tg> a)
    {
        return new(N, G, Coefs.ToDictionary(e => e.Key, e => e.Value + a[e.Key]));
    }

    public GZNElt<Tn, Tg> Sub(GZNElt<Tn, Tg> a)
    {
        return new(N, G, Coefs.ToDictionary(e => e.Key, e => e.Value - a[e.Key]));
    }

    public GZNElt<Tn, Tg> Opp()
    {
        return new(N, G, Coefs.ToDictionary(e => e.Key, e => e.Value.Opp()));
    }

    public GZNElt<Tn, Tg> Act(int k)
    {
        return new(N, G, Coefs.ToDictionary(e => e.Key, e => e.Value.Act(k)));
    }

    public GZNElt<Tn, Tg> Act(Tg g)
    {
        var map = new Dictionary<Tg, ZNElt<Tn>>(G.Count());
        foreach (var (g0, z) in Coefs)
            map[G.Op(g, g0)] = z;

        return new(N, G, map);
    }

    public GZNElt<Tn, Tg> Act(MapElt<Tg, Automorphism<Tn>> L)
    {
        var map = Zero.Coefs;
        foreach (var (g, z) in Coefs)
        {
            if (z.IsKnown())
            {
                if (z.V.Equals(N.Neutral()))
                    continue;

                var v0 = L[g][z.V];
                map[G.Neutral()] = map[G.Neutral()] + new ZNElt<Tn>(N, v0, z.Coefs);
            }
            else
                map[g] = map[g] + z;
        }

        return new(N, G, map);
    }

    public GZNElt<Tn, Tg> Substitute(Xi xi, Tn n)
    {
        var map = Zero.Coefs;
        foreach (var (g, z) in Coefs)
        {
            map[g] = z.Substitute(xi, n);
        }

        return new(N, G, map);
    }

    public GZNElt<Tn, Tg> Substitute((Xi xi, Tn n)[] batch)
    {
        var gz = this;
        return batch.Aggregate(gz, (acc, e) => acc.Substitute(e.xi, e.n));
    }

    public override int GetHashCode() => Hash;

    public override string ToString()
    {
        if (IsZero())
            return "0";

        var n = G.Neutral();
        return Coefs.Where(e => !e.Value.IsZero()).Select(e => e.Key.Equals(n) ? $"{e.Value}" : $"{e.Key}.({e.Value})").Glue(" + ");
    }

    public static GZNElt<Tn, Tg> operator +(GZNElt<Tn, Tg> a, GZNElt<Tn, Tg> b) => a.Add(b);
    public static GZNElt<Tn, Tg> operator -(GZNElt<Tn, Tg> a, GZNElt<Tn, Tg> b) => a.Sub(b);
    public static GZNElt<Tn, Tg> operator -(GZNElt<Tn, Tg> a) => a.Opp();
    public static GZNElt<Tn, Tg> operator *(Tg g, GZNElt<Tn, Tg> a) => a.Act(g);
    public static GZNElt<Tn, Tg> operator *(int k, GZNElt<Tn, Tg> a) => a.Act(k);
}
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Integers;
using FqPoly = FastGoat.Structures.VecSpace.Polynomial<FastGoat.Structures.VecSpace.EPoly<FastGoat.UserGroup.Integers.ZnInt>,
        FastGoat.Structures.VecSpace.Xi>;

namespace FastGoat.UserGroup.GModuleN;

public readonly struct ZNElt<Tn, Tg> : IElt<ZNElt<Tn, Tg>>
    where Tn : struct, IElt<Tn>
    where Tg : struct, IElt<Tg>
{
    public AbelianDirectSum<Tn> Nab { get; }
    public ConcreteGroup<Tn> N => Nab.Ab;
    public Dictionary<Tn, FqPoly> Coefs { get; }
    public Indeterminates<Xi> Indeterminates { get; }
    public MapElt<Tg, Automorphism<Tn>> L { get; }

    public ZNElt(EPoly<ZnInt> scalar, Indeterminates<Xi> ind, AbelianDirectSum<Tn> N0, MapElt<Tg, Automorphism<Tn>> L0)
    {
        (Nab, L) = (N0, L0);
        Indeterminates = ind;
        KZero = scalar.Zero;
        var zero = Zero;
        Coefs = Nab.Decomp.ToDictionary(e => e.g, _ => zero);
        Hash = (N, ind).GetHashCode();
    }

    public ZNElt(EPoly<ZnInt> scalar, Indeterminates<Xi> ind, AbelianDirectSum<Tn> N0, MapElt<Tg, Automorphism<Tn>> L0, Xi[] xs)
    {
        (Nab, L) = (N0, L0);
        Indeterminates = ind;
        KZero = scalar.Zero;
        var zero = Zero;
        if (xs.Length != Nab.Decomp.Count)
            throw new();

        Coefs = Nab.Decomp.Zip(xs.Select(x => zero.X(x))).ToDictionary(e => e.First.g, e => e.Second);
        Hash = (N, ind).GetHashCode();
    }


    public ZNElt(EPoly<ZnInt> scalar, Indeterminates<Xi> ind, AbelianDirectSum<Tn> N0, MapElt<Tg, Automorphism<Tn>> L0, Dictionary<Tn, FqPoly> map)
    {
        (Nab, L) = (N0, L0);
        Coefs = new(map);
        Indeterminates = ind;
        KZero = scalar.Zero;
        Hash = (N, ind).GetHashCode();
    }

    public EPoly<ZnInt> KZero { get; }

    public FqPoly Zero => new(Indeterminates, KZero);

    public int Hash { get; }
    public FqPoly this[Tn n] => Coefs.TryGetValue(n, out var e) ? e : Zero;

    public bool Equals(ZNElt<Tn, Tg> other)
    {
        var elt = this;
        return Nab.Decomp.All(n => other[n.g].Equals(elt[n.g]));
    }

    public int CompareTo(ZNElt<Tn, Tg> other)
    {
        var elt = this;
        return Nab.Decomp.Select(n => elt[n.g].CompareTo(other[n.g])).FirstOrDefault(k => k != 0, 0);
    }

    public ZNElt<Tn, Tg> Add(ZNElt<Tn, Tg> n)
    {
        var elt = this;
        var map = Nab.Decomp.Select(e => (e, elt[e.g] + n[e.g])).ToDictionary(e => e.e.g, e => e.Item2);
        return new(KZero, Indeterminates, Nab, L, map);
    }

    public ZNElt<Tn, Tg> Sub(ZNElt<Tn, Tg> n)
    {
        var elt = this;
        var map = Nab.Decomp.Select(e => (e, elt[e.g] - n[e.g])).ToDictionary(e => e.e.g, e => e.Item2);
        return new(KZero, Indeterminates, Nab, L, map);
    }

    public ZNElt<Tn, Tg> Opp()
    {
        var elt = this;
        var map = Nab.Decomp.Select(e => (e, -elt[e.g])).ToDictionary(e => e.e.g, e => e.Item2);
        return new(KZero, Indeterminates, Nab, L, map);
    }

    public ZNElt<Tn, Tg> Act(int k)
    {
        var elt = this;
        var map = Nab.Decomp.Select(e => (e, k * elt[e.g])).ToDictionary(e => e.e.g, e => e.Item2);
        return new(KZero, Indeterminates, Nab, L, map);
    }

    public ZNElt<Tn, Tg> Act(Tg g)
    {
        var Nab0 = Nab;
        var L0 = L;
        var coefs = Coefs;
        // g.(x1*a+x2*b)
        // x1*g.a+x2*g.b
        // x1*(y1*a+y2*b)+x2*(z1*a+z2*b)
        // (x1*y1+x2*z1)*a+(x1*y2+x2*z2)*b

        // arr
        // a->y1*a+y2*b
        // b->z1*a+z2*b
        var arr = Nab.Decomp.Select(e => (e.g, Nab0.GEltToCanMap(L0[g][e.g]))).ToDictionary(a => a.g, a => a.Item2);
        var map0 = coefs.ToDictionary(
            kv => kv.Key,
            kv => arr.Select(a => coefs[a.Key] * a.Value[kv.Key].K).Aggregate((e0, e1) => e0 + e1)
        );

        return new(KZero, Indeterminates, Nab, L, map0);
    }

    public bool IsZero() => Coefs.Values.All(e => e.IsZero() || e.LeadingDetails.lc.IsZero());
    public bool IsKnown() => Coefs.Values.All(e => e.IsZero() || e.LeadingDetails.lm.Degree == 0);

    public override int GetHashCode() => Hash;

    public override string ToString()
    {
        if (IsZero())
            return "0";

        var coefs = Coefs.OrderBy(e => e.Key).Where(e => !e.Value.IsZero()).Select(e => $"({e.Value})*{e.Key}");
        return coefs.Glue(" + ");
    }

    public static ZNElt<Tn, Tg> operator +(ZNElt<Tn, Tg> a, ZNElt<Tn, Tg> b) => a.Add(b);

    public static ZNElt<Tn, Tg> operator -(ZNElt<Tn, Tg> a, ZNElt<Tn, Tg> b) => a.Sub(b);

    public static ZNElt<Tn, Tg> operator -(ZNElt<Tn, Tg> a) => a.Opp();

    public static ZNElt<Tn, Tg> operator *(int k, ZNElt<Tn, Tg> b) => b.Act(k);

    public static ZNElt<Tn, Tg> operator *(Tg g, ZNElt<Tn, Tg> b) => b.Act(g);
}
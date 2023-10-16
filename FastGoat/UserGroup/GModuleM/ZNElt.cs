using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Integers;

namespace FastGoat.UserGroup.GModuleM;

public readonly struct ZNElt<T> : IGmoduleMElt<T, ZnInt, ZNElt<T>>, IElt<ZNElt<T>> where T : struct, IElt<T>
{
    public SortedList<Xi, ZnInt> Coefs { get; }
    public int Ord { get; }
    public T V { get; }
    public ConcreteGroup<T> N { get; }

    public ZNElt(ConcreteGroup<T> n, params Xi[] unknowns)
    {
        if (n.GroupType == GroupType.NonAbelianGroup)
            throw new GroupException(GroupExceptionType.OnlyAbelianGroups);

        N = n;
        V = N.Neutral();
        var ord = Ord = N.Count();
        Coefs = new(unknowns.ToDictionary(e => e, _ => ZnInt.ZnZero(ord)));
        Hash = N.Hash;
    }

    public ZNElt(ConcreteGroup<T> n, T v, IDictionary<Xi, ZnInt> map)
    {
        N = n;
        V = v;
        Ord = N.Count();
        Coefs = new(map);
        Hash = N.Hash;
    }

    public bool IsZero() => V.Equals(N.Neutral()) && IsKnown();
    public bool IsKnown() => NbUnknowns() == 0;
    public int NbUnknowns() => Coefs.Values.Count(k => !k.IsZero());

    public ZNElt<T> GetUnknown(Xi xi)
    {
        var z0 = Zero.Coefs;
        z0[xi] = ZnInt.ZnZero(Ord).One;
        return new(N, N.Neutral(), z0);
    }

    public bool Equals(ZNElt<T> other)
    {
        return Hash == other.Hash
               && V.Equals(other.V)
               && Coefs.All(e => e.Value.Equals(other[e.Key]));
    }

    public int CompareTo(ZNElt<T> other)
    {
        var compV = V.CompareTo(other.V);
        if (compV != 0)
            return compV;

        var compUnks = NbUnknowns().CompareTo(other.NbUnknowns());
        if (compUnks != 0)
            return compUnks;

        return -Coefs.Values.SequenceCompareTo(other.Coefs.Values);
    }

    public int Hash { get; }

    public ZnInt this[Xi xi] => Coefs[xi];

    public ZNElt<T> Zero
    {
        get
        {
            var ord = Ord;
            return new(N, N.Neutral(), Coefs.ToDictionary(e => e.Key, _ => ZnInt.ZnZero(ord)));
        }
    }

    public ZNElt<T> Add(ZNElt<T> a)
    {
        var map = Coefs.ToDictionary(e => e.Key, e => e.Value + a[e.Key]);
        return new(N, N.Op(V, a.V), map);
    }

    public ZNElt<T> Sub(ZNElt<T> a)
    {
        var map = Coefs.ToDictionary(e => e.Key, e => e.Value - a[e.Key]);
        return new(N, N.Op(V, N.Invert(a.V)), map);
    }

    public ZNElt<T> Act(int k)
    {
        var map = Coefs.ToDictionary(e => e.Key, e => k * e.Value);
        return new(N, N.Times(V, k), map);
    }

    public ZNElt<T> Act(ZnInt a) => Act(a.K);

    public ZNElt<T> Opp() => Act(-1);

    public ZNElt<T> Substitute(Xi xi, T n)
    {
        if (!Coefs.ContainsKey(xi))
            return this;

        var map = Coefs.ToDictionary(e => e.Key, e => e.Value);
        var k = map[xi];
        map[xi] = ZnInt.ZnZero(Ord);
        return new(N, N.Op(V, N.Times(n, k.K)), map);
    }

    public override int GetHashCode() => Hash;

    public override string ToString()
    {
        if (IsZero())
            return "0";

        var str = Coefs.Where(e => !e.Value.IsZero())
            .Select(e => e.Value.K == 1 ? $"{e.Key}" : $"{e.Value}*{e.Key}").Glue("+");

        return V.Equals(N.Neutral()) ? str : $"{V}+{str}";
    }

    public static ZNElt<T> operator +(ZNElt<T> a, ZNElt<T> b) => a.Add(b);
    public static ZNElt<T> operator -(ZNElt<T> a, ZNElt<T> b) => a.Sub(b);
    public static ZNElt<T> operator -(ZNElt<T> a) => a.Opp();
    public static ZNElt<T> operator *(int k, ZNElt<T> a) => a.Act(k);
    public static ZNElt<T> operator *(ZnInt k, ZNElt<T> a) => a.Act(k);
}
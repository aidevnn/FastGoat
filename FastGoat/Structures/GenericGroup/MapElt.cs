using FastGoat.Commons;

namespace FastGoat.Structures.GenericGroup;

public readonly struct MapElt<T1, T2> : IElt<MapElt<T1, T2>>, IMap<T1, T2> where T2 : struct, IElt<T2> where T1 : struct, IElt<T1>
{
    public ConcreteGroup<T2> G2 { get; }
    public Dictionary<T1, T2> map { get; }
    public Dictionary<T1, T2> Map => new(map);
    public MapElt(ConcreteGroup<T1> g1, ConcreteGroup<T2> g2)
    {
        G2 = g2;
        Domain = g1;
        map = Domain.ToDictionary(e => e, _ => g2.Neutral());
        Hash = (Domain.Count(), G2.Count()).GetHashCode();
    }

    public MapElt(ConcreteGroup<T1> g1, ConcreteGroup<T2> g2, Dictionary<T1, T2> map0)
    {
        (Domain, G2) = (g1, g2);
        map = new(map0);
        Hash = (Domain.Count(), G2.Count()).GetHashCode();
    }

    public MapElt<T1, T2> Clone() => new(Domain, G2, Map);

    public bool Equals(IMap<T1, T2>? other)
    {
        return other is not null && IMap<T1, T2>.EqualiltyMap(this, other);
    }

    public int CompareTo(IMap<T1, T2>? other)
    {
        if (other is null)
            return 1;

        return IMap<T1, T2>.CompareMap(this, other);
    }

    public int Hash { get; }
    public ConcreteGroup<T1> Domain { get; }

    public IEnumerable<T1> Kernel()
    {
        var n = G2.Neutral();
        foreach (var e in map.Keys)
        {
            if (this[e].Equals(n))
                yield return e;
        }
    }

    public IEnumerable<T2> Image() => map.Values;

    public T2 this[T1 index] => map[index];

    public bool Equals(MapElt<T1, T2> other) => IMap<T1, T2>.EqualiltyMap(this, other);

    public int CompareTo(MapElt<T1, T2> other) => IMap<T1, T2>.CompareMap(this, other);

    public override int GetHashCode() => Hash;
    public override string ToString()
    {
        return map.AscendingByKey().GlueMap("; ", "{0}->[{1}]");
    }
}

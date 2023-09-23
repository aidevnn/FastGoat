using FastGoat.Commons;

namespace FastGoat.Structures.GenericGroup;

public readonly struct MapGroups<T1, T2> : IElt<MapGroups<T1, T2>>, IMap<T1, T2> where T2 : struct, IElt<T2> where T1 : struct, IElt<T1>
{
    public ConcreteGroup<T2> G2 { get; }
    private Dictionary<T1, T2> map { get; }
    public Dictionary<T1, T2> Map => new(map);
    public MapGroups(ConcreteGroup<T1> g1, ConcreteGroup<T2> g2)
    {
        G2 = g2;
        Domain = g1;
        map = Domain.ToDictionary(e => e, _ => g2.Neutral());
        var hashMap = map.Aggregate(0, (acc, e) => (acc, e.Key, e.Value).GetHashCode());
        Hash = (Domain.Hash, G2.Hash, hashMap).GetHashCode();
    }

    public MapGroups(ConcreteGroup<T1> g1, ConcreteGroup<T2> g2, T2[] arr)
    {
        (Domain, G2) = (g1, g2);
        map = Domain.Select((e, i) => (e, i)).ToDictionary(e => e.e, e => arr[e.i]);
        var hashMap = map.Aggregate(0, (acc, e) => (acc, e.Key, e.Value).GetHashCode());
        Hash = (Domain.Hash, G2.Hash, hashMap).GetHashCode();
    }

    public MapGroups(ConcreteGroup<T1> g1, ConcreteGroup<T2> g2, Dictionary<T1, T2> map0)
    {
        (Domain, G2) = (g1, g2);
        map = new(map0);
        var hashMap = map.Aggregate(0, (acc, e) => (acc, e.Key, e.Value).GetHashCode());
        Hash = (Domain.Hash, G2.Hash, hashMap).GetHashCode();
    }

    public MapGroups<T1, T2> Create(T2[] arr) => new(Domain, G2, arr);

    public bool Equals(IMap<T1, T2>? other)
    {
        Console.WriteLine("!!!!!!!!");
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
        foreach (var e in Domain)
        {
            if (this[e].Equals(n))
                yield return e;
        }
    }

    public IEnumerable<T2> Image() => map.Values;

    public T2 this[T1 index] => map[index];

    public bool Equals(MapGroups<T1, T2> other) => IMap<T1, T2>.EqualiltyMap(this, other);

    public int CompareTo(MapGroups<T1, T2> other) => IMap<T1, T2>.CompareMap(this, other);

    public override int GetHashCode() => Hash;
    public override string ToString()
    {
        return map.GlueMap("; ", "{0}->[{1}]");
    }
}

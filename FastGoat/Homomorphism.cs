namespace FastGoat;

public readonly struct Homomorphism<T1, T2> : IMap<T1, T2> where T1 : struct, IElt<T1> where T2 : struct, IElt<T2>
{
    public IReadOnlyDictionary<T1, T2> HomMap { get; }

    public Homomorphism(ConcreteGroup<T1> domGroup, IDictionary<T1, T2> homMap)
    {
        Domain = domGroup;
        HomMap = new Dictionary<T1, T2>(homMap);
        Hash = HomMap.OrderBy(kp => kp.Key).Aggregate(domGroup.Hash, (acc, kp) => (acc, kp.Key, kp.Value).GetHashCode());
    }

    public bool Equals(IMap<T1, T2>? other) => other?.Hash == Hash;

    public int CompareTo(IMap<T1, T2>? other) => other is null ? 1 : this.CompareMapTo(other);

    public int Count => HomMap.Count;
    public int Hash { get; }
    public ConcreteGroup<T1> Domain { get; }
    public IEnumerable<T1> Kernel()
    {
        var n = this[Domain.Neutral()];
        return HomMap.Where(kp => kp.Value.Equals(n)).Select(kp => kp.Key);
    }

    public IEnumerable<T2> Image() => HomMap.Values.Distinct();

    public T2 this[T1 index] => HomMap[index];

    public override int GetHashCode() => Hash;
    public override string ToString()
    {
        return HomMap.AscendingByKey().GlueMap();
    }
}
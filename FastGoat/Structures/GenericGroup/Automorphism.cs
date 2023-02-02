using FastGoat.Commons;

namespace FastGoat.Structures.GenericGroup;

public readonly struct Automorphism<T> : IMap<T, T>, IElt<Automorphism<T>> where T : struct, IElt<T>
{
    public IReadOnlyDictionary<T, T> AutMap { get; }

    public Automorphism(AutomorphismGroup<T> aut)
    {
        AutMap = aut.G.ToDictionary(e => e, e => e);
        Hash = AutMap.OrderBy(kp => kp.Key).Aggregate(aut.G.Hash, (acc, kp) => (acc, kp.Key, kp.Value).GetHashCode());
        AutGroup = aut;
    }

    public Automorphism(AutomorphismGroup<T> aut, IReadOnlyDictionary<T, T> e)
    {
        AutMap = new Dictionary<T, T>(e);
        Hash = AutMap.OrderBy(kp => kp.Key).Aggregate(aut.G.Hash, (acc, kp) => (acc, kp.Key, kp.Value).GetHashCode());
        AutGroup = aut;
    }

    public bool Equals(Automorphism<T> other) => Hash == other.Hash;

    public int CompareTo(Automorphism<T> other) => IMap<T,T>.CompareMap(this, other);

    public int Hash { get; }
    public AutomorphismGroup<T> AutGroup { get; }

    public ConcreteGroup<T> Domain => AutGroup.G;

    public IEnumerable<T> Kernel()
    {
        var n = Domain.Neutral();
        return AutMap.Where(kp => kp.Value.Equals(n)).Select(kp => kp.Key);
    }

    public IEnumerable<T> Image() => AutMap.Values.Distinct();

    public bool Equals(IMap<T, T>? other) => other?.Hash == Hash;
    public int CompareTo(IMap<T, T>? other) => other is null ? 1 : IMap<T, T>.CompareMap(this, other);
    public override int GetHashCode() => Hash;
    public T this[T index] => AutMap[index];

    public override string ToString()
    {
        return AutMap.AscendingByKey().GlueMap();
    }
}
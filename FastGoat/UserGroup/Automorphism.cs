namespace FastGoat.UserGroup;

public readonly struct Automorphism<T> : IMap<T, T>, IElt<Automorphism<T>> where T : struct, IElt<T>
{
    public Dictionary<T, T> AutMap { get; }

    public Automorphism(AutomorphismGroup<T> aut)
    {
        AutMap = aut.G.ToDictionary(e => e, e => e);
        Hash = AutMap.OrderBy(kp => kp.Key).Aggregate(aut.G.Hash, (acc, kp) => (acc, kp.Key, kp.Value).GetHashCode());
        BaseGroup = aut;
    }

    public Automorphism(AutomorphismGroup<T> aut, Dictionary<T, T> e)
    {
        AutMap = new Dictionary<T, T>(e);
        Hash = AutMap.OrderBy(kp => kp.Key).Aggregate(aut.G.Hash, (acc, kp) => (acc, kp.Key, kp.Value).GetHashCode());
        BaseGroup = aut;
    }

    public bool Equals(Automorphism<T> other) => Hash == other.Hash;

    public int CompareTo(Automorphism<T> other)
    {
        var compKey = AutMap.Keys.Ascending().SequenceCompareTo(other.AutMap.Keys.Ascending());
        if (compKey != 0)
            return compKey;

        return AutMap.OrderBy(a => a.Key).Select(kp => kp.Value)
            .SequenceCompareTo(other.AutMap.OrderBy(a => a.Key).Select(kp => kp.Value));
    }

    public int Hash { get; }
    public IGroup<Automorphism<T>> BaseGroup { get; }

    public HashSet<T> Domain => AutMap.Keys.ToHashSet();
    public HashSet<T> Codomain => AutMap.Values.ToHashSet();
    public bool Equals(IMap<T, T>? other) => other?.Hash == Hash;
    public int CompareTo(IMap<T, T>? other) => -other?.CompareMapTo(this) ?? 1;
    public override int GetHashCode() => Hash;
    public T this[T index] => AutMap[index];

    public override string ToString()
    {
        return AutMap.OrderBy(kp => kp.Key).Select(kp => $"{kp.Key}->{kp.Value}").Glue(", ");
    }
}
namespace FastGoat.UserGroup;

public readonly struct Automorphism<T> : IElt<Automorphism<T>> where T : struct, IElt<T>
{
    public Dictionary<T, T> AutMap { get; }

    public Automorphism(AutomorphismGroup<T> aut)
    {
        AutMap = aut.G.ToDictionary(e => e, e => e);
        BaseGroup = aut;
    }

    public Automorphism(AutomorphismGroup<T> aut, Dictionary<T, T> e)
    {
        AutMap = new Dictionary<T, T>(e);
        BaseGroup = aut;
    }

    public bool Equals(Automorphism<T> other)
    {
        if (!BaseGroup.Equals(other.BaseGroup))
            return false;

        foreach (var kp in AutMap)
            if (!other.AutMap.ContainsKey(kp.Key) || !other.AutMap[kp.Key].Equals(kp.Value))
                return false;

        return true;
    }

    public int CompareTo(Automorphism<T> other)
    {
        var compKey = AutMap.Keys.Ascending().SequenceCompareTo(other.AutMap.Keys.Ascending());
        if (compKey != 0)
            return compKey;

        return AutMap.OrderBy(a => a.Key).Select(kp => kp.Value)
            .SequenceCompareTo(other.AutMap.OrderBy(a => a.Key).Select(kp => kp.Value));
    }

    public int Hash { get; } = 0;
    public IGroup<Automorphism<T>> BaseGroup { get; }
    public override int GetHashCode() => BaseGroup.Hash;

    public T this[T index] => AutMap[index];

    public override string ToString()
    {
        return AutMap.OrderBy(kp => kp.Key).Select(kp => $"{kp.Key}->{kp.Value}").Glue(", ");
    }
}
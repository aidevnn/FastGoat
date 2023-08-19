using System.Collections;

namespace FastGoat.Structures.GenericGroup;

public enum FGroupOp
{
    Additive,
    Multiplicative
}

public class OpGroup<K> : IGroup<K> where K : struct, IFieldElt<K>, IRingElt<K>, IElt<K>
{
    private K neutral { get; }
    public FGroupOp FGroupOp { get; }
    public OpGroup(string name, K scalar, FGroupOp fGroupOp)
    {
        Name = name;
        neutral = fGroupOp == FGroupOp.Additive ? scalar.Zero : scalar.One;
        FGroupOp = fGroupOp;
        Hash = (name, neutral, fGroupOp).GetHashCode();
    }

    public IEnumerator<K> GetEnumerator() => GetElements().GetEnumerator();

    IEnumerator IEnumerable.GetEnumerator()
    {
        return GetEnumerator();
    }

    public bool Equals(IGroup<K>? other) => other?.Hash == Hash;

    public int Hash { get; }
    public string Name { get; }

    public K this[params ValueType[] us]
    {
        get
        {
            if (us.Length != 1)
                throw new GroupException(GroupExceptionType.GroupDef);
            if (us[0] is int k)
                return FGroupOp == FGroupOp.Additive ? neutral + k : neutral * k;
            else
                throw new GroupException(GroupExceptionType.GroupDef);
        }
    }

    public IEnumerable<K> GetElements()
    {
        yield return neutral;
    }

    public IEnumerable<K> GetGenerators()
    {
        yield return neutral;
    }

    public K Neutral() => neutral;

    public K Invert(K e) => FGroupOp == FGroupOp.Additive ? e.Opp() : e.Inv();

    public K Op(K e1, K e2) => FGroupOp == FGroupOp.Additive ? e1.Add(e2) : e1.Mul(e2);
    public override int GetHashCode() => Hash;

    public override string ToString() => Name;
}
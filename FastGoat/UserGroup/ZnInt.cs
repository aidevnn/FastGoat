namespace FastGoat.UserGroup;

public readonly struct ZnInt : IElt<ZnInt>
{
    private Zn Zn { get; }
    public int K { get; }
    private string fmt { get; }

    public ZnInt(Zn zn, int k)
    {
        Zn = zn;
        K = k % zn.Mod;
        if (K < 0)
            K += zn.Mod;

        Hash = (K, zn.Hash).GetHashCode();
        var digits = $"{zn.Mod - 1}".Length;
        fmt = $"{{0,{digits}}}";
    }

    public bool Equals(ZnInt other)
    {
        return other.Hash == Hash;
    }

    public int CompareTo(ZnInt other)
    {
        if (!BaseGroup.Equals(other.BaseGroup))
            throw new GroupException(GroupExceptionType.BaseGroup);

        return K.CompareTo(other.K);
    }

    public int Hash { get; }
    public IGroup<ZnInt> BaseGroup => Zn;

    public override int GetHashCode()
    {
        return Hash;
    }

    public override string ToString()
    {
        return string.Format(fmt, K);
    }
}
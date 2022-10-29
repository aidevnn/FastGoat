using FastGoat.Theory;

namespace FastGoat.UserGroup.Integers;

public readonly struct ZnInt : IElt<ZnInt>
{
    public Zn Zn { get; }
    public int K { get; }

    public ZnInt(Zn zn, int k)
    {
        Zn = zn;
        K = k % zn.Mod;
        if (K < 0)
            K += zn.Mod;

        Hash = (K, zn.Hash).GetHashCode();
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
        return string.Format(Zn.Fmt, K);
    }

    public override bool Equals(object? obj)
    {
        return obj is ZnInt z && Equals(z);
    }

    public static bool operator ==(ZnInt a, ZnInt b) => a.Equals(b);
    public static bool operator !=(ZnInt a, ZnInt b) => !a.Equals(b);

    public static ZnInt operator +(ZnInt a, ZnInt b) => a.BaseGroup.Op(a, b);
    public static ZnInt operator *(ZnInt a, int k) => a.BaseGroup.Times(a, k);
}
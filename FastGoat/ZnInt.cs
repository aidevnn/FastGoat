using FastGoat;

public struct ZnInt : IElt<ZnInt>
{
    public int k { get; }
    public ZnInt()
    {
        throw new Exception();
    }
    public ZnInt(Zn zn)
    {
        k = 0;
        Group = Zn = zn;
        Hash = HashCode.Combine(zn.Hash, k);
    }
    public ZnInt(Zn zn, int k0)
    {
        k = k0 % zn.mod;
        if (k < 0)
            k += zn.mod;

        Group = Zn = zn;
        Hash = HashCode.Combine(zn, k);
    }
    public IGroup<ZnInt> Group { get; }
    public Zn Zn { get; }

    public int Hash { get; }

    public int CompareTo(ZnInt other)
    {
        if (!Group.Equals(other.Group))
            throw new BaseGroupException();

        return k.CompareTo(other.k);
    }

    public bool Equals(ZnInt other) => Group.Equals(other.Group) && k == other.k;
    public override int GetHashCode() => Hash;
    public override string ToString() => string.Format(Zn.fmt, k);
    public override bool Equals(object? obj)
    {
        if (obj is null) return false;
        return this.Equals((ZnInt)obj);
    }
    public static bool operator ==(ZnInt a, ZnInt b) => a.Equals(b);
    public static bool operator !=(ZnInt a, ZnInt b) => !a.Equals(b);
    public static implicit operator ZnInt((Zn zn, int k) p) => new(p.zn, p.k);
    public static ZnInt operator *(ZnInt a, ZnInt b) => a.Group.Op(a, b);
    public static ZnInt operator ^(ZnInt a, int p) => (a.Zn, a.k * p);
}

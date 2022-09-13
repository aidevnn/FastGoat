using FastGoat;

public struct Zn : IGroup<ZnInt>
{
    public string fmt { get; }
    public int mod { get; }
    public int Hash { get; }
    public Zn()
    {
        throw new Exception();
    }
    public Zn(int mod0)
    {
        if (mod0 < 2)
            throw new Exception();

        mod = mod0;
        Hash = mod;
        fmt = $"{{0,{mod.ToString().Length}}}";
    }

    public ZnInt Invert(ZnInt a)
    {
        if (!a.Zn.Equals(this))
            throw new Exception();

        return new(this, -a.k);
    }

    public ZnInt Neutral() => new(this);

    public ZnInt Op(ZnInt a, ZnInt b)
    {
        if (!a.Zn.Equals(this) || !b.Zn.Equals(this))
            throw new Exception();

        return new(this, a.k + b.k);
    }

    public bool Equals(IGroup<ZnInt>? other) => other?.Hash == Hash;

    public ZnInt this[int k] => new(this, k);
    public override int GetHashCode() => Hash;
    public override string ToString() => $"Z/{mod}Z";
}
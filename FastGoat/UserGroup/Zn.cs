using System.Collections;

namespace FastGoat.UserGroup;

public readonly struct Zn : IGroup<ZnInt>
{
    public int Mod { get; }
    public string Name { get; }
    public string Fmt { get; }

    public Zn(int mod)
    {
        if (mod < 2)
            throw new GroupException(GroupExceptionType.GroupDef);

        Mod = mod;
        Name = $"Z{mod}";
        var digits = $"{Mod - 1}".Length;
        Fmt = $"{{0,{digits}}}";
    }

    public bool Equals(IGroup<ZnInt>? other)
    {
        return other?.Hash == Hash;
    }

    public int Hash => Mod;

    public ZnInt Neutral()
    {
        return new(this, 0);
    }

    public ZnInt Invert(ZnInt e)
    {
        if (Mod != e.Zn.Mod)
            throw new GroupException(GroupExceptionType.BaseGroup);

        return new ZnInt(this, -e.K);
    }

    public ZnInt Op(ZnInt e1, ZnInt e2)
    {
        if (Mod != e1.Zn.Mod || Mod != e2.Zn.Mod)
            throw new GroupException(GroupExceptionType.BaseGroup);

        return new ZnInt(this, e1.K + e2.K);
    }

    public ZnInt this[params ValueType[] us]
    {
        get
        {
            try
            {
                var u = Convert.ToInt32(us[0]);
                return new ZnInt(this, u);
            }
            catch
            {
                throw new GroupException(GroupExceptionType.GroupDef);
            }
        }
    }

    public IEnumerable<ZnInt> GetElements()
    {
        yield return Neutral();
    }

    public IEnumerator<ZnInt> GetEnumerator()
    {
        return GetElements().GetEnumerator();
    }

    IEnumerator IEnumerable.GetEnumerator()
    {
        return GetElements().GetEnumerator();
    }

    public override int GetHashCode()
    {
        return Hash;
    }

    public override string ToString()
    {
        return Name;
    }
}
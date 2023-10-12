using System.Collections;
using FastGoat.Structures;

namespace FastGoat.UserGroup.Integers;

public readonly struct Zn : IGroup<ZnInt>
{
    public int Mod { get; }
    public string Name { get; }
    public string Fmt { get; }

    public Zn(int mod)
    {
        if (mod < 1)
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
        return new(Mod, 0);
    }

    public ZnInt Invert(ZnInt e)
    {
        if (Mod != e.Mod)
            throw new GroupException(GroupExceptionType.BaseGroup);

        return e.Opp();
    }

    public ZnInt Op(ZnInt e1, ZnInt e2)
    {
        if (Mod != e1.Mod || Mod != e2.Mod)
            throw new GroupException(GroupExceptionType.BaseGroup);

        return e1.Add(e2);
    }

    public ZnInt this[params ValueType[] us]
    {
        get
        {
            try
            {
                var u = Convert.ToInt32(us[0]);
                return new ZnInt(Mod, u);
            }
            catch
            {
                throw new GroupException(GroupExceptionType.GroupDef);
            }
        }
    }

    public IEnumerable<ZnInt> GetGenerators()
    {
        yield return new ZnInt(Mod, 1);
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
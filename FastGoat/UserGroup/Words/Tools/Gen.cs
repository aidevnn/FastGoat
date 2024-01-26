namespace FastGoat.UserGroup.Words.Tools;

public readonly struct Gen : IEquatable<Gen>, IComparable<Gen>
{
    public const char Id = 'i';
    public char V { get; }

    public Gen(char v0 = Id)
    {
        V = v0;
    }
    public bool Equals(Gen other) => other.V == V;

    private char ToLower() => V == Id ? '0' : char.ToLower(V);
    public int CompareTo(Gen other)
    {
        var comp = ToLower().CompareTo(other.ToLower());
        if (comp != 0)
            return comp;

        return -V.CompareTo(other.V);
    }

    public override int GetHashCode() => V;

    public Gen Invert() => V == Id ? this : new(char.IsLower(V) ? char.ToUpper(V) : char.ToLower(V));
    public override string ToString() => $"{V}";
}
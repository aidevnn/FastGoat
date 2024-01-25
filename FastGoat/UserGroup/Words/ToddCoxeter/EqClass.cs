using System.Diagnostics.CodeAnalysis;

namespace FastGoat.UserGroup.Words.ToddCoxeter;

public struct EqClass : IEquatable<EqClass>, IComparable<EqClass>
{
    public static EqClass Unknown => new();
    public static EqClass One => Unknown.Next;
    int Value { get; } = 0;

    private EqClass(int c)
    {
        Value = c;
    }

    
    public int Hash => Value;
    public EqClass Next => new(Value + 1);
    public bool Equals(EqClass other) => Value == other.Value;
    
    public override bool Equals([NotNullWhen(true)] object? obj)
    {
        return base.Equals(obj);
    }

    public int CompareTo(EqClass other) => Value.CompareTo(other.Value);
    public override int GetHashCode() => Value;
    public override string ToString() => Value == 0 ? " " : $"{Value}";

    public static EqClass Min(EqClass a, EqClass b) => a.CompareTo(b) <= 0 ? a : b;
    public static EqClass Max(EqClass a, EqClass b) => a.CompareTo(b) >= 0 ? a : b;
    public static (EqClass min, EqClass max) MinMax(EqClass a, EqClass b) => (Min(a, b), Max(a, b));
    public static implicit operator EqClass(int i) => new(i);
    public static bool operator ==(EqClass a, EqClass b) => a.Equals(b);
    public static bool operator !=(EqClass a, EqClass b) => !a.Equals(b);
}
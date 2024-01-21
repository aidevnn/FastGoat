using System.Diagnostics.CodeAnalysis;

namespace FastGoat.UserGroup.Words.ToddCoxeter;

public struct Generator : IEquatable<Generator>, IComparable<Generator>
{
    public static Generator Unknown = new();
    public static Generator Id = new('i');
    public char Value { get; }
    char i { get; }
    public int sgn { get; }

    public Generator(char c)
    {
        if (!char.IsLetter(c))
            throw new Exception();

        Value = c;
        i = char.IsLower(c) ? char.ToUpper(c) : char.ToLower(c);
        sgn = char.IsLower(c) ? 1 : -1;
    }

    public char GetLowerCase() => sgn == 1 ? Value : i;
    public Generator Invert() => i == Id.i ? this : new(i);
    public bool Equals(Generator other) => Value == other.Value;
    //
    // public override bool Equals([NotNullWhen(true)] object? obj)
    // {
    //     return base.Equals(obj);
    // }

    public int CompareTo(Generator other)
    {
        var compS = sgn.CompareTo(other.sgn);
        if (compS != 0)
            return -compS;

        return Value.CompareTo(other.Value);
    }

    public override int GetHashCode() => Value;
    public override string ToString() => Value == '\0' ? "?" : $"{Value}";
    public static implicit operator Generator(char c) => new(c);
    public static bool operator ==(Generator a, Generator b) => a.Equals(b);
    public static bool operator !=(Generator a, Generator b) => !a.Equals(b);
}
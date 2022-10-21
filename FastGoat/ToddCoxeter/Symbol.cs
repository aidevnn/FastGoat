using System.Diagnostics.CodeAnalysis;

namespace FastGoat.ToddCoxeter;

public struct Symbol : IEquatable<Symbol>, IComparable<Symbol>
{
    public static Symbol Unknown => new();
    public static Symbol One => Unknown.Next;
    int Value { get; } = 0;
    public Symbol(int c)
    {
        Value = c;
    }
    public Symbol Next => new(Value + 1);
    public bool Equals(Symbol other) => Value == other.Value;
    public override bool Equals([NotNullWhen(true)] object? obj)
    {
        return base.Equals(obj);
    }
    public int CompareTo(Symbol other) => Value.CompareTo(other.Value);
    public override int GetHashCode() => Value;
    public override string ToString() => Value == 0 ? " " : $"{Value}";

    public static Symbol Min(Symbol a, Symbol b) => a.CompareTo(b) <= 0 ? a : b;
    public static Symbol Max(Symbol a, Symbol b) => a.CompareTo(b) >= 0 ? a : b;
    public static (Symbol min, Symbol max) MinMax(Symbol a, Symbol b) => (Min(a, b), Max(a, b));
    public static implicit operator Symbol(int i) => new(i);
    public static bool operator ==(Symbol a, Symbol b) => a.Equals(b);
    public static bool operator !=(Symbol a, Symbol b) => !a.Equals(b);
}

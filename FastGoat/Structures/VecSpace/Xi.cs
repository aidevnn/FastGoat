using System.Collections;
using System.IO.IsolatedStorage;

namespace FastGoat.Structures.VecSpace;

public readonly struct Xi : IElt<Xi>
{
    public char xi { get; }

    public Xi()
    {
        xi = 'X';
        Hash = xi;
    }

    public Xi(int i)
    {
        var c = (char)('a' + i);
        if (!char.IsLetter(c))
            throw new ArgumentException();

        xi = c;
        Hash = xi;
    }

    public Xi(char c)
    {
        if (!char.IsLetter(c))
            throw new ArgumentException();

        xi = c;
        Hash = xi;
    }

    public bool Equals(Xi other) => Hash == other.Hash;

    public int CompareTo(Xi other) => xi.CompareTo(other.xi);

    public override int GetHashCode() => Hash;
    public override string ToString() => $"{xi}";
    public int Hash { get; }
}
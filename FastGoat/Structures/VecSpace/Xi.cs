namespace FastGoat.Structures.VecSpace;

public readonly struct Xi : IElt<Xi>
{
    public string xi { get; }

    public Xi()
    {
        xi = "X";
        Hash = xi.GetHashCode();
    }

    public Xi(int i)
    {
        var c = (char)('a' + i);
        if (!char.IsLetter(c))
            throw new ArgumentException();

        xi = ((char)c).ToString();
        Hash = xi.GetHashCode();
    }

    public Xi(char c)
    {
        if (!char.IsLetter(c))
            throw new ArgumentException();

        xi = c.ToString();
        Hash = xi.GetHashCode();
    }

    public Xi(string expr)
    {
        if (expr.Length == 0)
            throw new ArgumentException();

        xi = expr;
        Hash = xi.GetHashCode();
    }

    public bool Equals(Xi other) => Hash == other.Hash;

    public int CompareTo(Xi other) => String.Compare(xi, other.xi, StringComparison.Ordinal);

    public override int GetHashCode() => Hash;
    public override string ToString() => xi;
    public int Hash { get; }
}
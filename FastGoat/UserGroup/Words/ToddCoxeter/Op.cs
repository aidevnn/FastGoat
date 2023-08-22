namespace FastGoat.UserGroup.Words.ToddCoxeter;

public struct Op : IEquatable<Op>
{
    public static Op Unknown = new();
    public EqClass i { get; }
    public Generator g { get; }
    public EqClass j { get; }
    int hash { get; }

    public Op(EqClass i0, Generator g0, EqClass j0)
    {
        i = i0;
        g = g0;
        j = j0;
        hash = HashCode.Combine(i, g, j);
    }

    public override int GetHashCode() => hash;
    public bool Equals(Op other) => hash == other.hash;
    public Op Invert() => new(j, g.Invert(), i);

    public static Op Create(EqClass i0, Generator g0, EqClass j0) =>
        g0.sgn == 1 ? new(i0, g0, j0) : new(j0, g0.Invert(), i0);

    public bool Contain(EqClass s) => i == s || j == s;
    public override string ToString() => $"({i})·{g}=({j})";
    public static implicit operator Op((int i, char g, int j) p) => new(p.i, p.g, p.j);
}
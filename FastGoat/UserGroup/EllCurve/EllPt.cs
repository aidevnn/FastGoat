using FastGoat.Structures;

namespace FastGoat.UserGroup.EllCurve;

public readonly struct EllPt<T> : IElt<EllPt<T>> where T : struct, IElt<T>, IRingElt<T>, IFieldElt<T>
{
    public T X { get; }
    public T Y { get; }

    public int Hash { get; }
    public bool IsO { get; }

    public EllPt()
    {
        Hash = 0;
        IsO = true;
    }

    public EllPt(T x, T y)
    {
        IsO = false;
        (X, Y) = (x, y);
        Hash = (IsO, X0: X, Y0: Y).GetHashCode();
    }

    public bool Equals(EllPt<T> other)
    {
        return (IsO && other.IsO) || (X.Equals(other.X) && Y.Equals(other.Y));
    }

    public int CompareTo(EllPt<T> other)
    {
        if (Equals(other))
            return 0;

        if (IsO)
            return -1;

        if (other.IsO)
            return 1;

        return (X0: X, Y0: Y).CompareTo((other.X, other.Y));
    }

    public override int GetHashCode() => Hash;

    public override string ToString()
    {
        return IsO ? "O" : $"({X},{Y})";
    }

    public static implicit operator EllPt<T>((T x, T y) P) => new(P.x, P.y);
}
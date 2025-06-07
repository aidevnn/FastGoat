using FastGoat.Commons;
using FastGoat.Structures;

namespace FastGoat.UserGroup.EllCurve;

public readonly struct TriVar : IElt<TriVar>
{
    public TriVar()
    {
        (X3, X2, X1) = (0, 0, 0);
        Hash = (X3, X2, X1).GetHashCode();
    }

    public TriVar(int x3, int x2, int x1)
    {
        (X3, X2, X1) = (x3, x2, x1);
        Hash = (X3, X2, X1).GetHashCode();
    }

    public int X3 { get; }
    public int X2 { get; }
    public int X1 { get; }
    public int Degree => X3 + X2 + X1;
    public int Hash { get; }
    public bool IsOne() => X3 == 0 && X2 == 0 && X1 == 0;
    public TriVar Mul(TriVar e) => new(X3 + e.X3, X2 + e.X2, X1 + e.X1);
    public bool Equals(TriVar other) => X3 == other.X3 && X2 == other.X2 && X1 == other.X1;
    public int CompareTo(TriVar other) => (X3, X2, X1).CompareTo((other.X3, other.X2, other.X1));
    public override int GetHashCode() => Hash;
    public override string ToString() => $"{(X3, X2, X1)}";

    public void Deconstruct(out int x3, out int x2, out int x1)
    {
        (x3, x2, x1) = (X3, X2, X1);
    }

    public int this[int index]
    {
        get
        {
            return index switch
            {
                3 => X3,
                2 => X2,
                _ => X1
            };
        }
    }

    public TriVar Set(int idx, int v)
    {
        return idx switch
        {
            3 => new(v, X2, X1),
            2 => new(X3, v, X1),
            _ => new(X3, X2, v),
        };
    }

    public TriVar GetX3() => new(X3, 0, 0);
    public TriVar GetX2() => new(0, X2, 0);
    public TriVar GetX1() => new(0, 0, X1);
    public TriVar GetX3X2() => new(X3, X2, 0);
    public TriVar GetX3X1() => new(X3, 0, X1);
    public TriVar GetX2X1() => new(0, X2, X1);

    public static (TriVar pa, TriVar pb) Reduce(TriVar a, TriVar b)
    {
        int pax3 = 0, pax2 = 0, pax1 = 0;
        int pbx3 = 0, pbx2 = 0, pbx1 = 0;

        var mx3 = int.Max(a.X3, b.X3);
        if (mx3 != 0 && a.X3 != b.X3)
        {
            if (mx3 == b.X3)
                pax3 = mx3 - a.X3;
            else
                pbx3 = mx3 - b.X3;
        }

        var mx2 = int.Max(a.X2, b.X2);
        if (mx2 != 0 && a.X2 != b.X2)
        {
            if (mx2 == b.X2)
                pax2 = mx2 - a.X2;
            else
                pbx2 = mx2 - b.X2;
        }

        var mx1 = int.Max(a.X1, b.X1);
        if (mx1 != 0 && a.X1 != b.X1)
        {
            if (mx1 == b.X1)
                pax1 = mx1 - a.X1;
            else
                pbx1 = mx1 - b.X1;
        }

        return ((pax3, pax2, pax1), (pbx3, pbx2, pbx1));
    }

    public static implicit operator TriVar((int x3, int x2, int x1) e) => new(e.x3, e.x2, e.x1);
}
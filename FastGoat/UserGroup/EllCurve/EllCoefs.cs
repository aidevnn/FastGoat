using FastGoat.Commons;
using FastGoat.Structures;

namespace FastGoat.UserGroup.EllCurve;

public struct EllCoefs<K> where K : struct, IFieldElt<K>, IRingElt<K>, IElt<K>
{
    public EllCoefs(K a1, K a2, K a3, K a4, K a6)
    {
        Model = (a1, a2, a3, a4, a6);

        var b2 = a1 * a1 + 4 * a2;
        var b4 = a1 * a3 + 2 * a4;
        var b6 = a3 * a3 + 4 * a6;
        var b8 = a1 * a1 * a6 - a1 * a3 * a4 + 4 * a2 * a6 + a2 * a3 * a3 - a4 * a4;
        B_Invariants = (b2, b4, b6, b8);

        var c4 = b2 * b2 - 24 * b4;
        var c6 = -b2.Pow(3) + 36 * b2 * b4 - 216 * b6;
        C_Invariants = (c4, c6);

        var disc = -b2 * b2 * b8 - 8 * b4.Pow(3) - 27 * b6 * b6 + 9 * b2 * b4 * b6;
        var j = c4.Pow(3) / disc;
        (J_Invariant, Disc) = (j, disc);
    }

    public (K a1, K a2, K a3, K a4, K a6) Model { get; }
    public (K b2, K b4, K b6, K b8) B_Invariants { get; }
    public (K c4, K c6) C_Invariants { get; }
    public K J_Invariant { get; }

    public (K a1, K a2, K a3, K a4, K a6, K b2, K b4, K b6, K b8, K c4, K c6, K j) ModelAndInvariants
    {
        get
        {
            var (a1, a2, a3, a4, a6) = Model;
            var (b2, b4, b6, b8) = B_Invariants;
            var (c4, c6) = C_Invariants;
            var j = J_Invariant;
            return (a1, a2, a3, a4, a6, b2, b4, b6, b8, c4, c6, j);
        }
    }

    public K Disc { get; }

    public string Eq
    {
        get
        {
            var (x, y) = Ring.Polynomial(Disc, "x", "y").Deconstruct();
            var (a1, a2, a3, a4, a6) = Model;
            var lhs = y * y + a1 * x * y + a3 * y;
            var rhs = x.Pow(3) + a2 * x * x + a4 * x + a6;
            return $"Ellptic curve {lhs} = {rhs}";
        }
    }

    public string ModelStr => $"[{Model}]".Replace("(", "").Replace(")", "");

    public string B_InvariantsStr
    {
        get
        {
            var (b2, b4, b6, b8) = B_Invariants;
            return $"B Invariants b2={b2} b4={b4} b6={b6} b8={b8}";
        }
    }

    public string C_InvariantsStr
    {
        get
        {
            var (c4, c6) = C_Invariants;
            return $"C Invariants c4={c4} c6={c6}";
        }
    }

    public void Show()
    {
        Console.WriteLine(Eq);
        Console.WriteLine($"Disc = {Disc}");
        Console.WriteLine(B_InvariantsStr);
        Console.WriteLine(C_InvariantsStr);
        Console.WriteLine($"J Invariant j={J_Invariant}");
    }

    public EllCoefs<K> Transform(K r, K s, K t, K u)
    {
        var (a1, a2, a3, a4, a6) = Model;
        var a01 = (a1 + 2 * s) / u;
        var a02 = (a2 - s * a1 + 3 * r - s * s) / u.Pow(2);
        var a03 = (a3 + r * a1 + 2 * t) / u.Pow(3);
        var a04 = (a4 - s * a3 + 2 * r * a2 - (t + r * s) * a1 + 3 * r * r - 2 * s * t) / u.Pow(4);
        var a06 = (a6 + r * a4 + r * r * a2 + r.Pow(3) - t * a3 - t * t - r * t * a1) / u.Pow(6);
        return new(a01, a02, a03, a04, a06);
    }

    public EllCoefs<K> Transform(int r, int s, int t, int u)
    {
        var o = Disc.One;
        return Transform(r * o, s * o, t * o, u * o);
    }

    public EllGroup<K> ToEllGroup()
    {
        var (a1, a2, a3, a4, a6) = Model;
        return new(a1, a2, a3, a4, a6);
    }
}
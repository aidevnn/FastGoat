
namespace FastGoat.Structures.VecSpace;

public class TriVarInd
{
    public MonomOrder MonomOrder { get; } = MonomOrder.Lex;
    public string X1 { get; } = "X";
    public string X2 { get; } = "Y";
    public string X3 { get; } = "Z";
    private int Hash { get; }

    public TriVarInd()
    {
        Hash = (MonomOrder, X1, X2, X3).GetHashCode();
    }

    public TriVarInd(MonomOrder monomOrder)
    {
        MonomOrder = monomOrder;
        Hash = (MonomOrder, X1, X2, X3).GetHashCode();
    }

    public TriVarInd(MonomOrder monomOrder, string x, string y, string z)
    {
        (MonomOrder, X1, X2, X3) = (monomOrder, x, y, z);
        Hash = (MonomOrder, X1, X2, X3).GetHashCode();
    }

    public Comparer<TriVar> GetComparer()
    {
        return MonomOrder switch
        {
            MonomOrder.GrLex => comparerGrLex,
            MonomOrder.GrevLex => comparerGrevLex,
            _ => comparerLex
        };
    }

    public override int GetHashCode() => Hash;
    public override string ToString() => $"[Order={MonomOrder} X1={X1} X2={X2} X3={X3}]";

    static TriVarInd()
    {
        comparerLex = Comparer<TriVar>.Create((a, b) => a.CompareTo(b));
        comparerGrLex = Comparer<TriVar>.Create(
            (a, b) =>
            {
                var compDeg = a.Degree.CompareTo(b.Degree);
                if (compDeg != 0)
                    return compDeg;

                return a.CompareTo(b);
            });
        comparerGrevLex = Comparer<TriVar>.Create(
            (a, b) =>
            {
                var compDeg = a.Degree.CompareTo(b.Degree);
                if (compDeg != 0)
                    return compDeg;

                return -(a.CompareTo(b));
            });
    }

    private static Comparer<TriVar> comparerLex { get; }
    private static Comparer<TriVar> comparerGrLex { get; }
    private static Comparer<TriVar> comparerGrevLex { get; }
}
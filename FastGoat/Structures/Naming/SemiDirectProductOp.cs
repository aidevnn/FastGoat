using FastGoat.Commons;
using FastGoat.Structures.GenericGroup;

namespace FastGoat.Structures.Naming;

public class SemiDirectProductOp : ANameElt
{
    public ANameElt Lhs { get; }
    public ANameElt Rhs { get; }

    public SemiDirectProductOp(ANameElt lhs, ANameElt rhs, ConcreteGroup<TableElt> g)
    {
        ContentGroup = g;
        ContentType = NodeType.SemiDirectProduct;
        (Lhs, Rhs) = (lhs, rhs);
        Name = $"{Lhs.NameParenthesis} x: {Rhs.NameParenthesis}";
        var name2 = CyclicSdp(ContentGroup, Lhs.ContentGroup!, Rhs.ContentGroup!);
        Weight = Lhs.Weight * Rhs.Weight;
        Depth = 1 + int.Max(Lhs.Depth, Rhs.Depth);
        if (!string.IsNullOrEmpty(name2))
        {
            Name = name2;
            Weight = 25;
            Depth = 0;
        }
    }
    
    static Homomorphism<T, Automorphism<T>> GroupAction<T>(ConcreteGroup<T> G, ConcreteGroup<T> K, ConcreteGroup<T> H)
        where T : struct, IElt<T>
    {
        if (!G.SuperSetOf(K) || !G.SuperSetOf(H))
            throw new GroupException(GroupExceptionType.NotSubGroup);

        if (K.Intersect(H).Count() > 1)
            throw new();

        var conj = Group.ByConjugate(G);
        if (G.Grid2D(K).Any(e => !K.Contains(conj(e.t1, e.t2))))
            throw new GroupException(GroupExceptionType.NotNormal);

        var map1 = H.ToDictionary(g => g, g => K.ToDictionary(k => k, k => conj(g, k)));
        var autK = Group.AutBase(K);
        var map2 = map1.ToDictionary(e => e.Key, e => new Automorphism<T>(autK, e.Value));
        return new Homomorphism<T, Automorphism<T>>(H, map2);
    }

    static string GroupActionName<T>(Homomorphism<T, Automorphism<T>> act) where T : struct, IElt<T>
    {
        var K = act.Image().First().Domain;
        var H = act.Domain;

        if (act.Image().Distinct().Count() == 1)
            return $"{K.NameParenthesis()} x {H.NameParenthesis()}";

        var gensK = K.GetGenerators().ToArray();
        var gensH = H.GetGenerators().ToArray();

        var (gh, gk) = (gensH[0], gensK[0]);
        var aut = act[gh];
        var dicK = Group.Cycle(K, gk);
        var (m, n, r) = (K.Count(), H.Count(), dicK[aut[gk]]);

        if (n == 2)
        {
            if (m == 3)
                return "S3";
            else if (r == m - 1)
                return $"D{2 * m}";
            else if (int.IsPow2(m))
            {
                if (r == m / 2 - 1)
                    return $"QD{2 * m}";
                else if (r == m / 2 + 1)
                    return $"MM{2 * m}";
            }
        }

        return IntExt.Gcd(m, n * (r - 1)) == 1 ? $"F({m}x:{n}){r}" : $"M({m}x:{n}){r}";
    }

    static string CyclicSdp<T>(ConcreteGroup<T> G, ConcreteGroup<T> K, ConcreteGroup<T> H) where T : struct, IElt<T>
    {
        var gensK = K.GetGenerators().ToArray();
        var gensH = H.GetGenerators().ToArray();

        if (gensK.Length != 1 || gensH.Length != 1)
            return "";

        var act = GroupAction(G, K, H);
        var name = GroupActionName(act);
        return name;
    }
}
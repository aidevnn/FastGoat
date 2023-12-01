using System.Collections;
using System.ComponentModel;
using System.Diagnostics;
using System.Numerics;
using FastGoat.Commons;
using FastGoat.Examples;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Characters;
using FastGoat.UserGroup.GModuleN;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Matrix;
using FastGoat.UserGroup.Perms;
using FastGoat.UserGroup.Polynoms;
using FastGoat.UserGroup.Words;
using FastGoat.UserGroup.Words.ToddCoxeter;
using static FastGoat.Commons.IntExt;
using static FastGoat.Commons.EnumerableExt;
using FastGoat.UserGroup.Padic;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");

Homomorphism<T, Automorphism<T>> GroupAction<T>(ConcreteGroup<T> G, ConcreteGroup<T> K, ConcreteGroup<T> H) where T : struct, IElt<T>
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

string GroupActionName<T>(Homomorphism<T, Automorphism<T>> act) where T : struct, IElt<T>
{
    if (act.IsNull)
        return "";

    var K = act.Image().First().Domain;
    var H = act.Domain;

    if (act.Image().Distinct().Count() == 1)
        return $"{K.NameParenthesis} x {H.NameParenthesis}";

    var gensK = K.GetGenerators().ToArray();
    var gensH = H.GetGenerators().ToArray();

    if (gensK.Length == 1 && gensH.Length == 1)
    {
        var (gh, gk) = (gensH[0], gensK[0]);
        var aut = act[gh];
        var dicK = Group.Cycle(K, gk);
        var (m, n, r) = (K.Count(), H.Count(), dicK[aut[gk]]);

        if (n == 2)
        {
            if (r == m - 1)
                return $"D{2 * m}";
            else if (int.IsPow2(m))
            {
                if (r == m / 2 - 1)
                    return $"QD{2 * m}";
                else if (r == m / 2 + 1)
                    return $"MM{2 * m}";
            }
        }

        return Gcd(m, n * (r - 1)) == 1 ? $"F({m},{n},{r})" : $"Mc({m},{n},{r})";
    }

    return "";
}

void FindAction()
{
    for (int o = 6; o < 65; o++)
    {
        if (o % 4 == 0)
        {
            var dn = FG.DihedralSdp(o / 2);
            var (k, h) = (dn.Ncan, dn.Gcan);
            var act = GroupAction(dn, k, h);
            Console.WriteLine("{0} => {1}", dn.ShortName, GroupActionName(act));
        }

        var fbs = FG.FrobeniusSdp(o).Cast<SemiDirectProduct<ZnInt, ZnInt>>().ToArray();
        foreach (var fb in fbs)
        {
            var (k, h) = (fb.Ncan, fb.Gcan);
            var act = GroupAction(fb, k, h);
            Console.WriteLine("{0} => {1}", fb.ShortName, GroupActionName(act));
        }

        if (int.IsPow2(o))
        {
            {
                var qdn = FG.SemiDihedralSdp(int.Log2(o));
                var (k, h) = (qdn.Ncan, qdn.Gcan);
                var act = GroupAction(qdn, k, h);
                Console.WriteLine("{0} => {1}", qdn.ShortName, GroupActionName(act));
            }

            {
                var mmn = FG.ModularMaxSdp(int.Log2(o));
                var (k, h) = (mmn.Ncan, mmn.Gcan);
                var act = GroupAction(mmn, k, h);
                Console.WriteLine("{0} => {1}", mmn.ShortName, GroupActionName(act));
            }
        }

        Console.WriteLine();
    }
}

{
    FindAction();
}
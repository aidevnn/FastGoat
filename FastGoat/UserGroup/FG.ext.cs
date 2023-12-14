using FastGoat.Examples;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.UserGroup.GModuleN;

namespace FastGoat.UserGroup;

public record ExtInfos<Tn, Tg>(CrMap<Tn, Tg> c, ConcreteGroup<Ep2<Tn, Tg>> ext, AllSubgroups<Ep2<Tn, Tg>> allSubs)
    where Tg : struct, IElt<Tg>
    where Tn : struct, IElt<Tn>;

public static partial class FG
{
    static IEnumerable<ExtInfos<Tn, Tg>> AllExtensionsInternal<Tn, Tg>(ConcreteGroup<Tn> N, ConcreteGroup<Tg> G, int nbOps)
        where Tg : struct, IElt<Tg>
        where Tn : struct, IElt<Tn>
    {
        var CN = Group.Zentrum(N);
        var autN = Group.AutomorphismGroup(N);
        var ops = Group.AllHomomorphisms(G, autN);
        var set = new HashSet<AllSubgroups<Ep2<Tn, Tg>>>(new IsomorphSubGroupsInfosEquality<Ep2<Tn, Tg>>());
        foreach (var (op, i) in ops.Select((l, i) => (l, i + 1)).Take(nbOps))
        {
            var L = op.ToMapElt(autN);
            var lbl = $"Lbl{i}/{ops.Count}";
            var (cohs, cobs, cocs) = ZNSolver.ReduceCohomologies(CN, G, L, lbl: lbl);
            foreach (var c in cohs)
            {
                var c0 = c.ToMapElt;
                var ext = Group.ExtensionGroup(N, L, c0, G);
                var extBase = ((ExtensionGroupBase<Tn, Tg>)ext.BaseGroup);
                if (extBase.IsGroup)
                {
                    var allSubs = new AllSubgroups<Ep2<Tn, Tg>>(ext);
                    if (set.Add(allSubs))
                        yield return new(c, ext, allSubs);
                }
                else
                {
                    Console.WriteLine("????????????????????? Extension isnt a group"); // TODO FIX
                }
            }

            Console.WriteLine($"Nb Exts:{set.Count}");
        }
    }

    public static IEnumerable<ExtInfos<Tn, Tg>> AllExtensions<Tn, Tg>(params (int nbOps, ConcreteGroup<Tn>, ConcreteGroup<Tg>)[] tuples)
        where Tg : struct, IElt<Tg>
        where Tn : struct, IElt<Tn>
    {
        var set = new HashSet<AllSubgroups<Ep2<Tn, Tg>>>(new IsomorphSubGroupsInfosEquality<Ep2<Tn, Tg>>());
        foreach (var (nbOps, n, g) in tuples)
        {
            foreach (var extInfos in AllExtensionsInternal(n, g, nbOps))
            {
                if (set.Add(extInfos.allSubs))
                    yield return extInfos;
            }

            Console.WriteLine();
            Console.WriteLine($"Total Exts:{set.Count}");
        }
    }

    public static IEnumerable<ExtInfos<Tn, Tg>> AllExtensions<Tn, Tg>(params (ConcreteGroup<Tn>, ConcreteGroup<Tg>)[] tuples)
        where Tg : struct, IElt<Tg>
        where Tn : struct, IElt<Tn>
    {
        foreach (var extInfos in AllExtensions(tuples.Select(e => (10000, e.Item1, e.Item2)).ToArray()))
            yield return extInfos;
    }

    public static IEnumerable<ExtInfos<Tn, Tg>> AllExtensions<Tn, Tg>(int nbOps, params (ConcreteGroup<Tn>, ConcreteGroup<Tg>)[] tuples)
        where Tg : struct, IElt<Tg>
        where Tn : struct, IElt<Tn>
    {
        foreach (var extInfos in AllExtensions(tuples.Select(e => (nbOps, e.Item1, e.Item2)).ToArray()))
            yield return extInfos;
    }
}
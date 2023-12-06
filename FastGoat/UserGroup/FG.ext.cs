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
        var dicExts = new Dictionary<SubGroupsInfos, HashSet<ConcreteGroup<Ep2<Tn, Tg>>>>();
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
                    var infos = allSubs.Infos;
                    if (dicExts.ContainsKey(infos))
                    {
                        if (dicExts[infos].Add(ext))
                            yield return new(c, ext, allSubs);
                    }
                    else
                    {
                        dicExts[infos] = new(new IsomorphEquality<Ep2<Tn, Tg>>()) { ext };
                        yield return new(c, ext, allSubs);
                    }
                }
                else
                {
                    Console.WriteLine("????????????????????? Extension isnt a group"); // TODO FIX
                }
            }

            Console.WriteLine($"Nb Exts:{dicExts.Values.Sum(v => v.Count)}");
        }
    }

    public static IEnumerable<ExtInfos<Tn, Tg>> AllExtensions<Tn, Tg>(int nbOps, params (ConcreteGroup<Tn>, ConcreteGroup<Tg>)[] tuples)
        where Tg : struct, IElt<Tg>
        where Tn : struct, IElt<Tn>
    {
        var dicExts = new Dictionary<SubGroupsInfos, HashSet<ConcreteGroup<Ep2<Tn, Tg>>>>();
        foreach (var (n, g) in tuples)
        {
            foreach (var extInfos in AllExtensionsInternal(n, g, nbOps))
            {
                if (dicExts.ContainsKey(extInfos.allSubs.Infos))
                {
                    if (dicExts[extInfos.allSubs.Infos].Add(extInfos.ext))
                        yield return extInfos;
                }
                else
                {
                    dicExts[extInfos.allSubs.Infos] = new(new IsomorphEquality<Ep2<Tn, Tg>>()) { extInfos.ext };
                    yield return extInfos;
                }
            }

            Console.WriteLine();
            Console.WriteLine($"Total Exts:{dicExts.Values.Sum(v => v.Count)}");
        }
    }

    public static IEnumerable<ExtInfos<Tn, Tg>> AllExtensions<Tn, Tg>(params (ConcreteGroup<Tn>, ConcreteGroup<Tg>)[] tuples)
        where Tg : struct, IElt<Tg>
        where Tn : struct, IElt<Tn>
    {
        foreach (var extInfos in AllExtensions(nbOps: 10000, tuples))
            yield return extInfos;
    }
}
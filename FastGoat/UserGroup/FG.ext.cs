
using FastGoat.Commons;
using FastGoat.Examples;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.Naming;
using FastGoat.UserGroup.GModuleN;

namespace FastGoat.UserGroup;

public record ExtInfos<Tn, Tg>(CrMap<Tn, Tg> c, ConcreteGroup<Ep2<Tn, Tg>> ext, AllSubgroups<Ep2<Tn, Tg>> allSubs)
    : IComparable<ExtInfos<Tn, Tg>>
    where Tg : struct, IElt<Tg>
    where Tn : struct, IElt<Tn>
{
    private (GroupType, int, SubGroupsInfos) ToTuples() => (ext.GroupType, -ext.ElementsOrders.Values.Max(), allSubs.Infos);

    public int CompareTo(ExtInfos<Tn, Tg>? other)
    {
        if (other is null)
            return 1;

        return ToTuples().CompareTo(other.ToTuples());
    }
}

public static partial class FG
{
    static IEnumerable<ExtInfos<Tn, Tg>> AllExtensionsInternal<Tn, Tg>(ConcreteGroup<Tn> N, ConcreteGroup<Tg> G,
        int nbOps, int nbSkip = 0)
        where Tg : struct, IElt<Tg>
        where Tn : struct, IElt<Tn>
    {
        var CN = Group.Zentrum(N);
        var autN = Group.AutomorphismGroup(N);
        var ops = Group.AllHomomorphisms(G, autN);
        var set = new HashSet<AllSubgroups<Ep2<Tn, Tg>>>(new IsomorphSubGroupsInfosEquality<Ep2<Tn, Tg>>());

        foreach (var (op, i) in ops.Select((l, i) => (l, i + 1)).Skip(nbSkip).Take(nbOps))
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

    public static IEnumerable<ExtInfos<Tn, Tg>> AllExtensions<Tn, Tg>(
        params (int nbOps, ConcreteGroup<Tn>, ConcreteGroup<Tg>)[] tuples)
        where Tg : struct, IElt<Tg>
        where Tn : struct, IElt<Tn>
    {
        foreach (var ext in AllExtensions(tuples.Select(e => (e.nbOps, 0, e.Item2, e.Item3)).ToArray()))
            yield return ext;
    }

    public static IEnumerable<ExtInfos<Tn, Tg>> AllExtensions<Tn, Tg>(
        params (int nbOps, int nbSkip, ConcreteGroup<Tn>, ConcreteGroup<Tg>)[] tuples)
        where Tg : struct, IElt<Tg>
        where Tn : struct, IElt<Tn>
    {
        var set = new HashSet<AllSubgroups<Ep2<Tn, Tg>>>(new IsomorphSubGroupsInfosEquality<Ep2<Tn, Tg>>());
        foreach (var (nbOps, nbSkip, n, g) in tuples)
        {
            foreach (var extInfos in AllExtensionsInternal(n, g, nbOps, nbSkip))
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

    public static IEnumerable<ExtInfos<Tn, Tg>> AllExtensions<Tn, Tg>(int nbOps,
        params (ConcreteGroup<Tn>, ConcreteGroup<Tg>)[] tuples)
        where Tg : struct, IElt<Tg>
        where Tn : struct, IElt<Tn>
    {
        foreach (var extInfos in AllExtensions(tuples.Select(e => (nbOps, 0, e.Item1, e.Item2)).ToArray()))
            yield return extInfos;
    }

    public static IEnumerable<AllSubgroups<TableElt>> AppendIsomorphic(this IEnumerable<AllSubgroups<TableElt>> subgs1,
        params IEnumerable<AllSubgroups<TableElt>>[] subs2)
    {
        foreach (var sub in subs2.Prepend(subgs1).SelectMany(e=>e).FilterIsomorphic())
            yield return sub;
    }
    
    public static IEnumerable<AllSubgroups<TableElt>> FilterIsomorphic(this IEnumerable<AllSubgroups<TableElt>> subsgr)
    {
        var set = new HashSet<AllSubgroups<TableElt>>(new IsomorphSubGroupsInfosEquality<TableElt>());
        foreach (var sub in subsgr)
        {
            if (set.Add(sub))
                yield return sub;
        }
    }
    
    public static IEnumerable<(ExtInfos<Tn, Tg> exts, ANameElt[] names)> NamingExts<Tn, Tg>(this IEnumerable<ExtInfos<Tn, Tg>> allExts)
        where Tg : struct, IElt<Tg>
        where Tn : struct, IElt<Tn>
    {
        foreach (var extInfos in allExts)
        {
            var it = NamesTree.BuildName(extInfos.allSubs.ToTable());
            extInfos.ext.Name = it[0].Name;
            yield return (extInfos, it);
        }
    }

    public static IEnumerable<ConcreteGroup<TableElt>> Naming(this IEnumerable<ConcreteGroup<TableElt>> groups)
    {
        foreach (var g in groups)
        {
            var names = NamesTree.BuildName(g);
            g.Name = names[0].Name;
            yield return g;
        }
    }

    public static IEnumerable<(AllSubgroups<TableElt> subsg, ANameElt[] names)> Naming(this IEnumerable<AllSubgroups<TableElt>> subsg)
    {
        foreach (var sub in subsg)
        {
            var names = NamesTree.BuildName(sub);
            sub.Parent.Name = names[0].Name;
            yield return (sub, names);
        }
    }

    public static void DisplayExts<Tn, Tg>(this IEnumerable<(ExtInfos<Tn, Tg> exts, ANameElt[] names)> extsNames, bool lazy = false)
        where Tg : struct, IElt<Tg>
        where Tn : struct, IElt<Tn>
    {
        var nb = 0;
        if (lazy)
        {
            foreach (var (extInfos, names) in extsNames)
            {
                ++nb;
                DisplayExtName(extInfos.ext, extInfos.allSubs.Infos, names);
            }
        }
        else
        {
            var lt = extsNames.ToArray();
            var maxLt = lt.Max(e => e.exts.ext.Name.Length);
            nb = lt.Length;
            foreach (var (extInfos, names) in lt.OrderBy(e => e.exts))
                DisplayExtName(extInfos.ext, extInfos.allSubs.Infos, names, maxLt);
        }

        Console.WriteLine($"Total Exts:{nb}");
    }

    public static void DisplayNames(this IEnumerable<(AllSubgroups<TableElt>subsg, ANameElt[] names)> seq)
    {
        var lt = seq.OrderBy(e => e.subsg.Parent.GroupType)
            .ThenByDescending(e => e.subsg.Parent.ElementsOrders.Values.Max())
            .ThenBy(e => e.subsg.Infos)
            .ToArray();
        var maxLt = lt.Max(e => e.subsg.Parent.Name.Length);
        var nb = lt.Length;
        foreach (var (subsg, names) in lt)
            DisplayExtName(subsg.Parent, subsg.Infos, names, maxLt);
        
        Console.WriteLine($"Total Groups:{nb}");
    }

    public static void DisplayExtsNames<Tn, Tg>(this IEnumerable<ExtInfos<Tn, Tg>> allExts)
        where Tg : struct, IElt<Tg>
        where Tn : struct, IElt<Tn>
    {
        foreach (var extInfos in allExts.Order())
        {
            var names = NamesTree.BuildName(extInfos.allSubs.ToTable());
            DisplayExtName(extInfos.ext, extInfos.allSubs.Infos, names);
        }
    }

    private static void DisplayExtName<T>(ConcreteGroup<T> ext, SubGroupsInfos infos, ANameElt[] names, int maxLt = -1)
        where T : struct, IElt<T>
    {
        var name = ext.Name = names[0].Name;
        maxLt = maxLt == -1 ? name.Length : maxLt;
        var lt = Enumerable.Repeat('#', maxLt + 4).Glue();
        var line = $"#################{lt}#################";
        var fmt = $"#################  {{0,{-maxLt}}}  #################";
        Console.WriteLine(line);
        Console.WriteLine(fmt, ext.Name);
        Console.WriteLine(line);
        DisplayGroup.HeadOrdersNames(ext, infos, names);
    }

    public static ConcreteGroup<TableElt>[] ToTableExts<Tn, Tg>(this IEnumerable<ExtInfos<Tn, Tg>> allExts,
        bool naming = false, string prefix = "Ext")
        where Tg : struct, IElt<Tg>
        where Tn : struct, IElt<Tn>
    {
        var lt = new List<ConcreteGroup<TableElt>>();
        foreach (var (ext, k) in allExts.Select((e, k) => (e, k + 1)))
        {
            lt.Add(ext.ext.ToTable());
            if (naming)
                ext.ext.Name = $"{prefix}{ext.ext.Count()}[{k}]";
        }

        return lt.ToArray();
    }
}
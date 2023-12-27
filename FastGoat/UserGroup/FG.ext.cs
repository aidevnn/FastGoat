
using System.Collections.Concurrent;
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
        var autG = Group.AutomorphismGroup(G);
        var allOps = Group.AllHomomorphisms(G, autN);
        var ops = allOps.ToHashSet(new OpByAutEquality<Tn, Tg>(G, autG, autN));
        
        Console.WriteLine();
        Console.WriteLine($"AutG:{autG.Count()} AutN:{autN.Count()}");
        Console.WriteLine($"AllOps:{allOps.Count} Filtered:{ops.Count}");

        foreach (var (op, i) in ops.Select((l, i) => (l, i + 1)).Skip(nbSkip).Take(nbOps))
        {
            var L = op.ToMapElt(autN);
            var lbl = $"Lbl{i}/{ops.Count}";
            var (cohs, cobs, cocs) = ZNSolver.ReduceCohomologies(CN, G, L, lbl: lbl);
            foreach (var c in cohs)
            {
                var c0 = c.ToMapElt;
                var ext = Group.ExtensionGroup(N, L, c0, G);
                var allSubs = new AllSubgroups<Ep2<Tn, Tg>>(ext);
                yield return new(c, ext, allSubs);
            }
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
        foreach (var (nbOps, nbSkip, n, g) in tuples)
        {
            foreach (var extInfos in AllExtensionsInternal(n, g, nbOps, nbSkip))
                yield return extInfos;
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
        foreach (var sub in subs2.Prepend(subgs1).SelectMany(e => e).FilterIsomorphic())
            yield return sub;
    }
    
    public static IEnumerable<AllSubgroups<TableElt>> FilterIsomorphic(this IEnumerable<AllSubgroups<TableElt>> subsgr)
    {
        var nb = 1;
        var set = new HashSet<AllSubgroups<TableElt>>(1000, new IsomorphSubGroupsInfosEquality<TableElt>());
        Console.WriteLine("## Start New Filter");
        foreach (var sub in subsgr)
        {
            if (set.Add(sub))
            {
                var name = $"    Iso{sub.Parent.Count()} no:{nb++}";
                Console.WriteLine(name);
                yield return sub;
            }
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

    public static (AllSubgroups<TableElt> subsg, ANameElt[] names)[] 
        DisplayNames(this IEnumerable<(AllSubgroups<TableElt>subsg, ANameElt[] names)> seq)
    {
        var lt = seq.OrderBy(e => e.subsg.Parent.GroupType)
            .ThenByDescending(e => e.subsg.Parent.ElementsOrders.Values.Max())
            .ThenBy(e => e.subsg.Infos)
            .ToArray();
        var maxLt = lt.Max(e => e.subsg.Parent.Name.Length);
        var nb = lt.Length;
        foreach (var (subsg, names) in lt)
            DisplayName(subsg.Parent, subsg.Infos, names, maxLt);
        
        Console.WriteLine($"Total Groups:{nb}");
        return lt;
    }

    public static AllSubgroups<TableElt>[] DisplayBox(this IEnumerable<AllSubgroups<TableElt>> seq)
    {
        var list = seq.OrderBy(e => e.Parent.GroupType)
            .ThenByDescending(e => e.Parent.ElementsOrders.Values.Max())
            .ThenBy(e => e.Infos)
            .ToArray();
        var maxLt = list.Max(e => e.Parent.Name.Length);
        var nb = 0;
        var nbSharp = 16;
        foreach (var subsg in list)
        {
            var name = subsg.Parent.Name = $"Grp{subsg.Parent.Count()}[{++nb}]";
            var diff = (maxLt - name.Length) / 2;
            var space = Enumerable.Repeat(' ', diff).Glue();
            var lt = Enumerable.Repeat('#', maxLt + 4).Glue();
            var sharp = Enumerable.Repeat('#', nbSharp).Glue();
            var line = $"{sharp}{lt}{sharp}";
            var fmt = $"{sharp}{space}  {{0,{-maxLt + diff * 2}}}  {space}{sharp}";
            Console.WriteLine(line);
            Console.WriteLine(fmt, name);
            Console.WriteLine(line);
            DisplayGroup.HeadOrders(subsg.Parent);
            Console.WriteLine(subsg.Infos);
            Console.WriteLine();
        }
        
        Console.WriteLine($"Total Groups:{nb}");
        return list;
    }

    private static void DisplayName<T>(ConcreteGroup<T> g, SubGroupsInfos infos, ANameElt[] names, int maxLt = -1)
        where T : struct, IElt<T>
    {
        var nbSharp = 16;
        var name = g.Name = names[0].Name;
        maxLt = int.Max(name.Length, maxLt);
        var diff = (maxLt - name.Length) / 2;
        var space = Enumerable.Repeat(' ', diff).Glue();
        var lt = Enumerable.Repeat('#', maxLt + 4).Glue();
        var sharp = Enumerable.Repeat('#', nbSharp).Glue();
        var line = $"{sharp}{lt}{sharp}";
        var fmt = $"{sharp}{space}  {{0,{-maxLt + diff * 2}}}  {space}{sharp}";
        Console.WriteLine(line);
        Console.WriteLine(fmt, g.Name);
        Console.WriteLine(line);
        DisplayGroup.HeadOrdersNames(g, infos, names);
    }
}

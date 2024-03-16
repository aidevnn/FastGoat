using System.Collections.Concurrent;
using FastGoat.Commons;
using FastGoat.Examples;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.Naming;
using FastGoat.Structures.Subgroups;
using FastGoat.UserGroup.DatabaseSmallGroups;
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
    public static EqualityComparer<Homomorphism<Tg, Automorphism<Tn>>>
        EqOpByAut<Tn, Tg>(ConcreteGroup<Tg> G, ConcreteGroup<Automorphism<Tg>> autG, ConcreteGroup<Automorphism<Tn>> autN)
        where Tn : struct, IElt<Tn>
        where Tg : struct, IElt<Tg>
    {
        var N = ((AutomorphismGroup<Tn>)autN.BaseGroup).G;
        var GN = G.Grid2D(N).Select(e => (g: e.t1, n: e.t2)).ToArray();
        return EqualityComparer<Homomorphism<Tg, Automorphism<Tn>>>.Create(
            (x, y) =>
            {
                return autG.Any(f => GN.All(e => x[e.g][e.n].Equals(y[f[e.g]][e.n])))
                       || autN.Any(f =>
                       {
                           var fi = f.Invert();
                           return GN.All(e => f[x[e.g][fi[e.n]]].Equals(y[e.g][e.n]));
                       });
            },
            obj => obj.Count
        );
    }

    static EqualityComparer<AllSubgroups<T>> EqSubGroups<T>() where T : struct, IElt<T>
    {
        return EqualityComparer<AllSubgroups<T>>.Create(
            (x, y) => x.Parent.IsIsomorphicTo(y.Parent),
            obj => (obj.Parent.Count(), obj.Parent.GroupType, obj.Infos).GetHashCode()
        );
    }
    
    static IEnumerable<ExtInfos<Tn, Tg>> AllExtensionsInternal<Tn, Tg>(ConcreteGroup<Tn> N, ConcreteGroup<Tg> G,
        int nbOps, int nbSkip = 0)
        where Tg : struct, IElt<Tg>
        where Tn : struct, IElt<Tn>
    {
        var CN = Group.Zentrum(N);
        var autN = Group.AutomorphismGroup(N);
        var autG = Group.AutomorphismGroup(G);
        var allOps = Group.AllHomomorphisms(G, autN);
        var ops = allOps.ToHashSet(EqOpByAut(G, autG, autN));
        
        if (Logger.Level != LogLevel.Off)
        {
            Console.WriteLine();
            Console.WriteLine($"AutG:{autG.Count()} AutN:{autN.Count()}");
            Console.WriteLine($"AllOps:{allOps.Count} Filtered:{ops.Count}");
        }

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

    public static IEnumerable<SemiDirectProduct<T1, T2>> AllSDPFilter<T1, T2>(ConcreteGroup<T1> N, ConcreteGroup<T2> G, 
        bool trivial = false)
        where T1 : struct, IElt<T1>
        where T2 : struct, IElt<T2>
    {
        if (Logger.Level != LogLevel.Off)
            Console.WriteLine($"############### AllSDP {N.NameParenthesis()} x: {G.NameParenthesis()}");
        
        var autG = Group.AutomorphismGroup(G);
        var autN = Group.AutomorphismGroup(N);
        var allOps = Group.AllHomomorphisms(G, autN);
        var ops = allOps.Where(kp => trivial || kp.Image().Count() > 1).ToHashSet(EqOpByAut(G, autG, autN));
        
        var nb = ops.Count();
        if (Logger.Level != LogLevel.Off)
        {
            Console.WriteLine($"AutG:{autG.Count()} AutN:{autN.Count()}");
            Console.WriteLine($"AllOps:{allOps.Count} remaining:{nb}");
        }
        var k = 1;
        foreach (var theta in ops)
        {
            if (Logger.Level != LogLevel.Off)
                Console.WriteLine($"  ## {k++,3}/{nb} ##");
            
            yield return Group.SemiDirectProd(N, theta, G);
        }
    }

    public static IEnumerable<SemiDirectProduct<T1, T2>> AllSDPFilterLazy<T1, T2>(ConcreteGroup<T1> N, ConcreteGroup<T2> G,
        bool trivial = false)
        where T1 : struct, IElt<T1>
        where T2 : struct, IElt<T2>
    {
        if (Logger.Level != LogLevel.Off)
            Console.WriteLine($"############### AllSDP {N.NameParenthesis()} x: {G.NameParenthesis()}");
        
        var autG = Group.AutomorphismGroup(G);
        var autN = Group.AutomorphismGroup(N);
        var allOps = Group.AllHomomorphisms(G, autN);
        if (Logger.Level != LogLevel.Off)
            Console.WriteLine($"AutG:{autG.Count()} AutN:{autN.Count()}");
        
        var ops = allOps.Where(kp => trivial || kp.Image().Count() > 1).Distinct(EqOpByAut(G, autG, autN));
        var k = 1;
        foreach (var theta in ops)
        {
            if (Logger.Level != LogLevel.Off)
                Console.WriteLine($"  ##   {k++,3}   ##");
            
            yield return Group.SemiDirectProd(N, theta, G);
        }
    }

    public static IEnumerable<AllSubgroups<WElt>> AppendIsomorphic(this IEnumerable<AllSubgroups<WElt>> subgs1,
        params IEnumerable<AllSubgroups<WElt>>[] subs2)
    {
        foreach (var sub in subs2.Prepend(subgs1).SelectMany(e => e).FilterIsomorphic())
            yield return sub;
    }

    public static IEnumerable<AllSubgroups<T>> FilterIsomorphic<T>(this IEnumerable<AllSubgroups<T>> subsgr)
        where T : struct, IElt<T>
    {
        var dic = new Dictionary<int, int>();
        var set = new HashSet<AllSubgroups<T>>(3000, EqSubGroups<T>());
        var nbSubs = new Dictionary<int, Dictionary<SubGroupsInfos, int>>(130);
        if (Logger.Level != LogLevel.Off)
            Console.WriteLine("## Start New Filter");
        foreach (var sub in subsgr)
        {
            var og = sub.Parent.Count();
            if (!dic.ContainsKey(og))
                dic[og] = 0;

            if (!nbSubs.ContainsKey(og))
                nbSubs[og] = new Dictionary<SubGroupsInfos, int>(og * og);

            if (!nbSubs[og].ContainsKey(sub.Infos))
                nbSubs[og][sub.Infos] = 0;

            if (nbSubs[og][sub.Infos] == NbSubGroups(og, sub.Infos))
                continue;

            if (set.Add(sub))
            {
                nbSubs[og][sub.Infos]++;
                if (Logger.Level != LogLevel.Off)
                {
                    var ids = allIds[og].Where(e => e.Infos == sub.Infos).Select(e => e.No).Glue(",", "{0:000}");
                    var name = $"    Iso{sub.Parent.Count()} no:{++dic[og]}/{GroupExt.A000001[og]} [{ids}]:{nbSubs[og][sub.Infos]}";
                    Console.WriteLine(name);
                }

                yield return sub;
            }
        }
    }

    public static IEnumerable<ExtInfos<Tn, Tg>> FilterIsomorphic<Tn, Tg>(this IEnumerable<ExtInfos<Tn, Tg>> subsgr)
        where Tn : struct, IElt<Tn>
        where Tg : struct, IElt<Tg>
    {
        var nb = 1;
        var set = new HashSet<AllSubgroups<Ep2<Tn, Tg>>>(3000, EqSubGroups<Ep2<Tn, Tg>>());
        Console.WriteLine("## Start New Filter");
        foreach (var sub in subsgr)
        {
            if (set.Add(sub.allSubs))
            {
                var name = $"    Iso{sub.allSubs.Parent.Count()} no:{nb++}";
                Console.WriteLine(name);
                yield return sub;
            }
        }
    }

    public static IEnumerable<(AllSubgroups<T> subsg, ANameElt[] names)> Naming<T>(this IEnumerable<ConcreteGroup<T>> groups,
        bool rename = true)
        where T : struct, IElt<T>
    {
        foreach (var g in groups)
        {
            var subs = g.AllSubgroups();
            var names = NamesTree.BuildName(subs.ToGroupWrapper());
            if (rename)
                g.Name = names[0].Name;
            yield return (subs, names);
        }
    }

    public static IEnumerable<(AllSubgroups<T> subsg, ANameElt[] names)> Naming<T>(this IEnumerable<AllSubgroups<T>> subsg,
        bool rename = true)
        where T : struct, IElt<T>
    {
        foreach (var sub in subsg)
        {
            if (sub is AllSubgroups<WElt> subTb)
            {
                var names = NamesTree.BuildName(subTb);
                if (rename)
                    sub.Parent.Name = names[0].Name;
                yield return (sub, names);
            }
            else
            {
                var names = NamesTree.BuildName(sub.ToGroupWrapper());
                if (rename)
                    sub.Parent.Name = names[0].Name;
                yield return (sub, names);
            }
        }
    }

    public static (AllSubgroups<T> subsg, ANameElt[] names)[]
        DisplayNames<T>(this IEnumerable<(AllSubgroups<T>subsg, ANameElt[] names)> seq, bool rename = false) where T : struct, IElt<T>
    {
        var lt = seq.OrderBy(e=>e.subsg.Parent.Count())
            .ThenBy(e => e.subsg.Parent.GroupType)
            .ThenBy(e => e.names[0])
            .ThenByDescending(e => e.subsg.Parent.ElementsOrders.Values.Max())
            .ThenBy(e => e.subsg.Infos)
            .ToArray();
        var dicOrd = lt.Select(e => e.subsg.Parent.Count()).Distinct().ToDictionary(k => k, _ => 0);
        var maxLt = rename ? lt.Max(e => e.names[0].Name.Length) : lt.Max(e => e.subsg.Parent.Name.Length);
        var nb = lt.Length;
        foreach (var (subsg, names) in lt)
        {
            var o = subsg.Parent.Count();
            Console.WriteLine($"Group{o}[{++dicOrd[o]}]");
            DisplayName(subsg.Parent, subsg.Infos, names, rename, maxLt);
        }

        Console.WriteLine($"Total Groups:{nb}");
        return lt;
    }

    public static AllSubgroups<T>[] DisplayBoxes<T>(this IEnumerable<AllSubgroups<T>> seq, bool rename = false)
        where T : struct, IElt<T>
    {
        var list = seq.OrderBy(e => e.Parent.GroupType)
            .ThenByDescending(e => e.Parent.ElementsOrders.Values.Max())
            .ThenBy(e => e.Infos)
            .ToArray();

        var ord = list.Select(e => e.Parent.Count()).First();
        var maxLt = rename ? ($"Grp{ord}[{list.Length}]").Length + 2 : list.Max(e => e.Parent.Name.Length + 2);
        var nb = 0;
        foreach (var subsg in list)
            DisplayBox(subsg, ++nb, rename, maxLt);

        Console.WriteLine($"Total Groups:{nb}");
        Console.WriteLine();
        return list;
    }

    public static void DisplayBox<T>(AllSubgroups<T> subsg, int nb, bool rename = false, int maxLt = -1) where T : struct, IElt<T>
    {
        var nbSharp = 16;
        if (rename)
            subsg.Parent.Name = $"Grp{subsg.Parent.Count()}[{nb}]";

        var name = subsg.Parent.Name;
        var max = int.Max(maxLt, name.Length);
        var diff = (max - name.Length) / 2;
        var space = Enumerable.Repeat(' ', diff).Glue();
        var lt = Enumerable.Repeat('#', max + 4).Glue();
        var sharp = Enumerable.Repeat('#', nbSharp).Glue();
        var line = $"{sharp}{lt}{sharp}";
        var fmt = $"{sharp}{space}  {{0,{-max + diff * 2}}}  {space}{sharp}";
        Console.WriteLine(line);
        Console.WriteLine(fmt, name);
        Console.WriteLine(line);
        DisplayGroup.HeadOrders(subsg.Parent);
        Console.CursorTop--;
        Console.WriteLine(subsg.Infos);
        var gapInfos = FindIdGroup(subsg.Parent, subsg.Infos);
        var s = gapInfos.Length > 1 ? " (TODO)" : "";
        foreach (var e in gapInfos)
            Console.WriteLine($"{$"Gap SmallGroup({e.Order},{e.No})",-24} Name:{e.Name}{s}");

        Console.WriteLine();
    }

    public static void DisplayName<T>(ConcreteGroup<T> g, SubGroupsInfos infos, ANameElt[] names, bool rename = false, int maxLt = -1)
        where T : struct, IElt<T>
    {
        var nbSharp = 16;
        if (rename)
            g.Name = names[0].Name;

        var name = g.Name;
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
        Console.CursorTop--;
        var gapInfos = FindIdGroup(g, infos);
        var s = gapInfos.Length > 1 ? " (TODO)" : "";
        foreach (var e in gapInfos)
            Console.WriteLine($"{$"Gap SmallGroup({e.Order},{e.No})",-24} Name:{e.Name}{s}");

        Console.WriteLine();
    }

    private static Dictionary<int, IdGroup[]> allIds { get; }
    private static Dictionary<int, Dictionary<SubGroupsInfos, int>> nbSubGroupsDetails { get; }

    public static int NbSubGroups(int ord, SubGroupsInfos infos)
    {
        return nbSubGroupsDetails[ord][infos];
    }

    public static IdGroup[] AllIds(int o) => allIds.ContainsKey(o) ? allIds[o] : [];

    public static IdGroup[] FindIdGroup<T>(ConcreteGroup<T> g, SubGroupsInfos infos) where T : struct, IElt<T>
    {
        var ord = g.Count();
        if (ord == 32)
        {
            if (infos.ToTuples() == (42, 30, 20))
            {
                var g3244 = Group.AllSemiDirectProd(SemiDihedralSdp(4), Abelian(2))
                    .First(e => e.AllSubgroups().Infos.ToTuples() == (42, 30, 20));
                if (g3244.IsIsomorphicTo(g))
                    return [AllIds(ord).First(e => e.No == 44)];
                else
                    return [AllIds(ord).First(e => e.No == 33)];
            }

            if (infos.ToTuples() == (26, 18, 14))
            {
                var g3213 = MetaCyclicSdp(8, 4, 3);
                if (g3213.IsIsomorphicTo(g))
                    return [AllIds(ord).First(e => e.No == 13)];
                else
                    return [AllIds(ord).First(e => e.No == 14)];
            }
        }

        if (ord == 42)
        {
            if (infos.ToTuples() == (20, 8, 6))
            {
                var g424 = Product.Generate(Abelian(3), DihedralSdp(7));
                if (g424.IsIsomorphicTo(g))
                    return [AllIds(ord).First(e => e.No == 4)];
                else
                    return [AllIds(ord).First(e => e.No == 2)];
            }
        }

        if (ord == 72)
        {
            if (infos.ToTuples() == (24, 24, 24))
            {
                var g729Type = Group.AbelianGroupType(g);
                if (g729Type.SequenceEqual([36, 2]))
                    return [AllIds(ord).First(e => e.No == 9)];
                else
                    return [AllIds(ord).First(e => e.No == 14)];
            }

            if (infos.ToTuples() == (48, 48, 48))
            {
                var g729Type = Group.AbelianGroupType(g);
                if (g729Type.SequenceEqual([18, 2, 2]))
                    return [AllIds(ord).First(e => e.No == 18)];
                else
                    return [AllIds(ord).First(e => e.No == 36)];
            }
        }

        // TODO ord 64, 72, 96, 110, 114

        return AllIds(ord).Where(e => e.Infos == infos).ToArray();
    }
}
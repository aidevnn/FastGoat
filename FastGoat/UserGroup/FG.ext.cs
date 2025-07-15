using System.Collections.Concurrent;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.Naming;
using FastGoat.Structures.Subgroups;
using FastGoat.UserGroup.DatabaseSmallGroups;
using FastGoat.UserGroup.GModuleN;
using FastGoat.UserGroup.Words.Tools;

namespace FastGoat.UserGroup;

public record ExtInfos<Tn, Tg>(CrMap<Tn, Tg> c, ConcreteGroup<Ep2<Tn, Tg>> ext, AllSubgroups<Ep2<Tn, Tg>> allSubs)
    : IComparable<ExtInfos<Tn, Tg>>
    where Tg : struct, IElt<Tg>
    where Tn : struct, IElt<Tn>
{
    private (GroupType, int, SubGroupsInfos) ToTuples() =>
        (ext.GroupType, -ext.ElementsOrders.Values.Max(), allSubs.Infos);

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
        EqOpByAut<Tn, Tg>(ConcreteGroup<Tg> G, ConcreteGroup<Automorphism<Tg>> autG,
            ConcreteGroup<Automorphism<Tn>> autN)
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

    public static IEnumerable<ExtInfos<Tn, Tg>> AllExtensions<Tn, Tg>(
        params (ConcreteGroup<Tn>, ConcreteGroup<Tg>)[] tuples)
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

    public static IEnumerable<SemiDirectProduct<T1, T2>> AllSDPFilterLazy<T1, T2>(ConcreteGroup<T1> N,
        ConcreteGroup<T2> G,
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
                    if (og <= GroupExt.A000001.Length)
                    {
                        var ids = allIds[og].Where(e => e.Infos == sub.Infos).Select(e => e.No).Glue(",", "{0:000}");
                        var name =
                            $"    Iso{sub.Parent.Count()} no:{++dic[og]}/{GroupExt.A000001[og]} [{ids}]:{nbSubs[og][sub.Infos]}";
                        Console.WriteLine(name);
                    }
                    else
                    {
                        var name = $"    Iso{sub.Parent.Count()} no:{++dic[og]}";
                        Console.WriteLine(name);
                    }
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

    public static IEnumerable<(AllSubgroups<T> subsg, ANameElt[] names)> Naming<T>(
        this IEnumerable<ConcreteGroup<T>> groups,
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

    public static IEnumerable<(AllSubgroups<T> subsg, ANameElt[] names)> Naming<T>(
        this IEnumerable<AllSubgroups<T>> subsg,
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
        DisplayNames<T>(this IEnumerable<(AllSubgroups<T> subsg, ANameElt[] names)> seq, bool rename = false,
            bool showBasegroup = true, bool showGenerators = true)
        where T : struct, IElt<T>
    {
        var lt = seq.OrderBy(e => e.subsg.Parent.Count())
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
            var s = GroupExt.A000001.Length >= o ? $"/{GroupExt.A000001[o]}" : "";
            Console.WriteLine($"Group{o}[{++dicOrd[o]}{s}]");
            DisplayName(subsg.Parent, subsg, names, rename, showBasegroup, showGenerators, maxLt);
        }

        Console.WriteLine($"Total Groups:{nb}");
        return lt;
    }

    public static AllSubgroups<T>[] DisplayBoxes<T>(this IEnumerable<AllSubgroups<T>> seq, bool rename = false,
        bool showGenerators = false)
        where T : struct, IElt<T>
    {
        var list = seq.OrderBy(e => e.Parent.Count())
            .ThenBy(e => e.Parent.GroupType)
            .ThenByDescending(e => e.Parent.ElementsOrders.Values.Max())
            .ThenBy(e => e.Infos)
            .ToArray();

        var ord = list.Select(e => e.Parent.Count()).First();
        var maxLt = rename ? ($"Grp{ord}[{list.Length}]").Length + 2 : list.Max(e => e.Parent.Name.Length + 2);
        var dicOrd = list.Select(e => e.Parent.Count()).Distinct().ToDictionary(k => k, _ => 0);
        var nb = list.Length;
        var k = 0;
        foreach (var subsg in list)
        {
            var o = subsg.Parent.Count();
            var s = GroupExt.A000001.Length >= o ? $"/{GroupExt.A000001[o]}" : "";
            Console.WriteLine($"Group{o}[{++dicOrd[o]}{s}]");
            DisplayBox(subsg, ++k, rename, false, showGenerators, maxLt);
        }

        Console.WriteLine($"Total Groups:{nb}");
        return list;
    }

    public static void DisplayDetails<T>(AllSubgroups<T> subsg, bool showBasegroup, bool showGenerators = true)
        where T : struct, IElt<T>
    {
        var g = subsg.Parent;
        var derived = subsg.GetDerivedSerie();
        var upper = subsg.GetUpperSerie();
        var lower = subsg.GetLowerSerie();
        var digits = derived.Content.Concat(upper.Content).Concat(lower.Content).Max(s => s.Representative.Name.Length);
        derived = new(derived.Content, derived.SerieType, digits, " --> ");
        upper = new(upper.Content, upper.SerieType, digits, " --> ");
        lower = new(lower.Content, lower.SerieType, digits, " --> ");
        var isSolvable = derived.Content.Last().Order == 1;
        var isNilpotent = upper.Content[0].Order == g.Count();
        var simplicity = subsg.IsSimple() ? "Simple" : "NotSimple";
        var nilpotency = isNilpotent ? "Nilpotent" : "NotNilpotent";
        var solubility = isSolvable ? "Solvable" : "NotSolvable";

        Console.WriteLine(g.ShortName);
        Console.WriteLine($"{simplicity}, {g.GroupType}, {nilpotency}, {solubility}");
        if (showBasegroup)
            Console.WriteLine($"BaseGroup {g.BaseGroup}");
        Console.WriteLine();

        DisplayGroup.Orders(g);
        Console.WriteLine(subsg.Infos);

        Console.WriteLine("Lower Serie");
        Console.WriteLine(lower);
        Console.WriteLine("Upper Serie");
        Console.WriteLine(upper);
        Console.WriteLine("Derived Serie");
        Console.WriteLine(derived);
        Console.WriteLine($"Zentrum  Z(G) = {subsg.ZentrumSubGroup()}");
        Console.WriteLine($"Frattini Φ(G) = {subsg.FrattiniSubGroup()}");
        Console.WriteLine($"Fitting  F(G) = {subsg.FittingSubGroup()}");
        Console.WriteLine();

        var lvl = Logger.Level;
        Logger.Level = LogLevel.Off;
        var rels = Graph.DefiningRelatorsOfGroup(g);
        var gens = rels.Where(c => char.IsLetter(c)).Distinct().Select(c => char.ToLower(c)).Order().ToArray();
        var def = $"< {gens.Glue(",")} | {rels.Replace(" ", "").Replace(",", ", ").Replace("=", " = ")} >";
        Console.WriteLine("Word Group");
        Console.WriteLine(def);
        Console.WriteLine();

        if (showGenerators)
            DisplayGroup.Generators(g);

        var gapInfos = FindIdGroup(g, subsg.Infos);
        var s = gapInfos.Length > 1 ? " (TODO)" : "";
        foreach (var e in gapInfos)
            Console.WriteLine($"{$"Gap SmallGroup({e.Order},{e.No})",-24} Name:{e.Name}{s}");

        Logger.Level = lvl;
        Console.WriteLine();
    }

    public static void DisplayBox<T>(AllSubgroups<T> subsg, int nb, bool rename = false,
        bool showBasegroup = true, bool showGenerators = true, int maxLt = -1) where T : struct, IElt<T>
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

        DisplayDetails(subsg, showBasegroup, showGenerators);
    }

    public static void DisplayName<T>(ConcreteGroup<T> g, AllSubgroups<T> subsg, ANameElt[] names, bool rename = false,
        bool showBasegroup = true, bool showGenerators = true, int maxLt = -1)
        where T : struct, IElt<T>
    {
        if (rename)
            subsg.Parent.Name = g.Name = names[0].Name;
        else
            subsg.Parent.Name = g.Name;

        DisplayBox(subsg, 0, rename, showBasegroup, showGenerators, maxLt);
        names.Println("Group names");
        Console.WriteLine();
    }

    public static (AllSubgroups<T> subsg, ANameElt[] names)[]
        CheckMissings<T>(this IEnumerable<(AllSubgroups<T> subsg, ANameElt[] names)> seq)
        where T : struct, IElt<T>
    {
        var seq0 = seq.ToArray();
        var rem = seq0.Select(s => (s.subsg.Parent.Count(), s.subsg.Parent, s.subsg))
            .Where(e => e.Item1 < 129)
            .GroupBy(e => e.Item1)
            .Select(e => (e.Key, e.SelectMany(f => FindIdGroup(f.Parent, f.subsg.Infos)).ToArray()))
            .ToDictionary(e => e.Key, e => AllIds(e.Key).Except(e.Item2).ToArray());

        var rem2 = rem.Where(e => e.Value.Length != 0).ToArray();
        foreach (var (ord, ids) in rem2)
            ids.Println(e => e.FullName, $"Missing Order:{ord}");

        Console.WriteLine($"Total Missings:{rem2.Sum(e => e.Value.Length)}");

        return seq0;
    }

    private static Dictionary<int, IdGroup[]> allIds { get; }
    private static Dictionary<int, Dictionary<SubGroupsInfos, int>> nbSubGroupsDetails { get; }

    public static int NbSubGroups(int ord, SubGroupsInfos infos)
    {
        return nbSubGroupsDetails.ContainsKey(ord) ? nbSubGroupsDetails[ord][infos] : -1;
    }

    public static IdGroup[] AllIds(int o) => allIds.ContainsKey(o) ? allIds[o] : new IdGroup[0];

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
                    return new[] { AllIds(ord).First(e => e.No == 44) };
                else
                    return new[] { AllIds(ord).First(e => e.No == 33) };
            }

            if (infos.ToTuples() == (26, 18, 14))
            {
                var g3213 = MetaCyclicSdp(8, 4, 3);
                if (g3213.IsIsomorphicTo(g))
                    return new[] { AllIds(ord).First(e => e.No == 13) };
                else
                    return new[] { AllIds(ord).First(e => e.No == 14) };
            }
        }

        if (ord == 42)
        {
            if (infos.ToTuples() == (20, 8, 6))
            {
                var g424 = Product.Generate(Abelian(3), DihedralSdp(7));
                if (g424.IsIsomorphicTo(g))
                    return new[] { AllIds(ord).First(e => e.No == 4) };
                else
                    return new[] { AllIds(ord).First(e => e.No == 2) };
            }
        }

        if (ord == 64)
        {
            if (infos.ToTuples() == (37, 27, 21))
            {
                var g6416 = MetaCyclicSdp(8, 8, 7);
                if (g6416.IsIsomorphicTo(g))
                    return new[] { AllIds(ord).First(e => e.No == 16) };
                else
                    return new[] { AllIds(ord).First(e => e.No == 15) };
            }

            if (infos.ToTuples() == (45, 23, 17))
            {
                var g6447 = MetaCyclicSdp(16, 4, 15);
                if (g6447.IsIsomorphicTo(g))
                    return new[] { AllIds(ord).First(e => e.No == 47) };
                else
                    return new[] { AllIds(ord).First(e => e.No == 48) };
            }

            if (infos.ToTuples() == (57, 32, 19))
            {
                var g6413 = WordGroup("b4, a4b-2, a2ba2b-1, ababababab-1ab-1a-1b-1a-1b-1");
                if (g6413.IsIsomorphicTo(g))
                    return new[] { AllIds(ord).First(e => e.No == 13) };
                else
                    return new[] { AllIds(ord).First(e => e.No == 14) };
            }

            if (infos.ToTuples() == (89, 52, 27))
            {
                var g64164 = WordGroup("b4, c2, cbcb-1, ca2ca-2, a2ba2b-1, ab2ca-1c, a3ba-1b-1");
                if (g64164.IsIsomorphicTo(g))
                    return new[] { AllIds(ord).First(e => e.No == 164) };
                else
                    return new[] { AllIds(ord).First(e => e.No == 165) };
            }

            if (infos.ToTuples() == (97, 65, 49))
            {
                var g64106 = WordGroup("b4, c2, caca-1, cbcb-1, a2ba2b-1, a3ba-1b-1");
                if (g64106.IsIsomorphicTo(g))
                    return new[] { AllIds(ord).First(e => e.No == 106) };
                else
                    return new[] { AllIds(ord).First(e => e.No == 107) };
            }

            if (infos.ToTuples() == (113, 57, 27))
            {
                var g64162 = WordGroup("b4, c2, abcabc, cbcb-1, ca2ca-2, a2ba2b-1, a3ba-1b-1");
                if (g64162.IsIsomorphicTo(g))
                    return new[] { AllIds(ord).First(e => e.No == 162) };
                else
                    return new[] { AllIds(ord).First(e => e.No == 161) };
            }

            if (infos.ToTuples() == (101, 54, 29))
            {
                var g64157 = WordGroup("b4, b2c2, bcbc, caca-1, a3ba-1b-1");
                if (g64157.IsIsomorphicTo(g))
                    return new[] { AllIds(ord).First(e => e.No == 157) };
                else
                    return new[] { AllIds(ord).First(e => e.No == 155) };
            }

            if (infos.ToTuples() == (77, 48, 29))
            {
                var g64156 = WordGroup("b4, c4, caca-1, ababc-2, a2ba2b-1, cbc-1b-1, a3ba-1b-1");
                if (g64156.IsIsomorphicTo(g))
                    return new[] { AllIds(ord).First(e => e.No == 156) };
                else
                    return new[] { AllIds(ord).First(e => e.No == 158) };
            }

            if (infos.ToTuples() == (81, 49, 33))
            {
                var g64179 = WordGroup("b4, b2c2, bcbc-1, bab-1a-1, a3ca-1c-1");
                if (g64179.IsIsomorphicTo(g))
                    return new[] { AllIds(ord).First(e => e.No == 179) };
                else
                    return new[] { AllIds(ord).First(e => e.No == 181) };
            }

            if (infos.ToTuples() == (145, 133, 121))
            {
                var g64248 = WordGroup("b2, c2, d2, cdcd, a3cac, a3dad, baba-1, abcba-1c, abdba-1d");
                if (g64248.IsIsomorphicTo(g))
                    return new[] { AllIds(ord).First(e => e.No == 248) };
                else
                    return new[] { AllIds(ord).First(e => e.No == 200) };
            }

            if (infos.ToTuples() == (165, 93, 39))
            {
                var g6480 = WordGroup("a4, b4, c4, baba-1, caca-1, b2cb2c-1, cbc2b-1c, a2bcb-1c-1");
                if (g6480.IsIsomorphicTo(g))
                    return new[] { AllIds(ord).First(e => e.No == 80) };
                else
                    return new[] { AllIds(ord).First(e => e.No == 77) };
            }

            if (infos.ToTuples() == (181, 118, 75))
            {
                var g64236 =
                    WordGroup("a4, b4, c2, d2, adad, a2bdbd, a2ca2c, a2cdcd, b2db2d, a2bcb-1c, ab2ca-1c, bab-1a-1");
                if (g64236.IsIsomorphicTo(g))
                    return new[] { AllIds(ord).First(e => e.No == 236) };
                else
                    return new[] { AllIds(ord).First(e => e.No == 240) };
            }

            if (infos.ToTuples() == (69, 45, 27))
            {
                var g6425 = WordGroup("a8, a4b4, bab2a-1b, ba2b-1a-2, ababab-1ab-1");
                if (g6425.IsIsomorphicTo(g))
                    return new[] { AllIds(ord).First(e => e.No == 25) };
                else
                    return new[] { AllIds(ord).First(e => e.No == 160) };
            }
        }

        if (ord == 72)
        {
            if (infos.ToTuples() == (24, 24, 24))
            {
                var g729Type = Group.AbelianGroupType(g);
                if (g729Type.SequenceEqual([36, 2]))
                    return new[] { AllIds(ord).First(e => e.No == 9) };
                else
                    return new[] { AllIds(ord).First(e => e.No == 14) };
            }

            if (infos.ToTuples() == (48, 48, 48))
            {
                var g729Type = Group.AbelianGroupType(g);
                if (g729Type.SequenceEqual([18, 2, 2]))
                    return new[] { AllIds(ord).First(e => e.No == 18) };
                else
                    return new[] { AllIds(ord).First(e => e.No == 36) };
            }
        }

        // TODO ord 64, 72, 96, 110, 114

        return AllIds(ord).Where(e => e.Infos == infos).ToArray();
    }
}
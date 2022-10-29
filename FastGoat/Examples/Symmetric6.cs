using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.GenericGroup;
using FastGoat.UserGroup.Perms;

namespace FastGoat.Examples;

public static class Symmetric6
{
    static void AddToSubGroups<T>(ConcreteGroup<T> g0, ConcreteGroup<T> g1, HashSet<HashSet<T>> all)
        where T : struct, IElt<T>
    {
        var allHoms = Group.AllHomomorphisms(g0, g1);
        // Console.WriteLine($"{g0}, {g1}, allHoms:{allHoms.Count()}");
        foreach (var m in allHoms)
        {
            var im = m.Image().ToHashSet();
            var ker = m.Kernel().ToHashSet();
            if (im.Count > 0 && im.All(g1.Contains))
                all.Add(im);

            if (ker.Count > 0 && ker.All(g1.Contains))
                all.Add(ker);
        }
    }

    // # All SubGroups of Symm5 : 156 Time:816 ms
    // [1, 1][2, 25][3, 10][4, 35][5, 6][6, 30][8, 15][10, 6][12, 15][20, 6][24, 5][60, 1][120, 1]
    public static void SubGroupsA5()
    {
        GlobalStopWatch.Restart();
        var gb = new Symm(5);
        var g0 = gb;
        var allSubGroups = new HashSet<HashSet<Perm>>(10000, new SetEquality<Perm>());
        foreach (var e1 in g0)
        {
            foreach (var e2 in g0)
            {
                var g1 = Group.GenerateElements(gb, e1, e2).ToHashSet();
                allSubGroups.Add(g1);
            }
        }

        GlobalStopWatch.Show($"All SubGroups of {g0.Name} : {allSubGroups.Count}");
        var allSorted = allSubGroups.GroupBy(sg => sg.Count).ToDictionary(a => a.Key, a => a.ToList());
        Console.WriteLine(allSorted.AscendingByKey().Select(p => $"[{p.Key}, {p.Value.Count}]").Glue());
    }

    // # All SubGroups of A6 : 501 Time:19739 ms
    // [1, 1][2, 45][3, 40][4, 75][5, 36][6, 120][8, 45][9, 10][10, 36][12, 30][18, 10][24, 30][36, 10][60, 12][360, 1]
    public static void SubGroupsA6()
    {
        GlobalStopWatch.Restart();
        var gb = new Symm(6);
        var g0 = Group.Generate("A6", gb, gb[(4, 5, 6)], gb[(1, 2, 3, 4, 5)]);
        var allSubGroups = new HashSet<HashSet<Perm>>(10000, new SetEquality<Perm>());
        foreach (var e1 in g0)
        {
            foreach (var e2 in g0)
            {
                var g1 = Group.GenerateElements(gb, e1, e2).ToHashSet();
                allSubGroups.Add(g1);
            }
        }

        GlobalStopWatch.Show($"All SubGroups of {g0.Name} : {allSubGroups.Count}");
        var allSorted = allSubGroups.GroupBy(sg => sg.Count).ToDictionary(a => a.Key, a => a.ToList());
        Console.WriteLine(allSorted.AscendingByKey().Select(p => $"[{p.Key}, {p.Value.Count}]").Glue());

        var g36 = Group.Generate("E36", g0, gb[(4, 5, 6)], gb[(1, 4), (2, 5, 3, 6)]);
        var allSubGroupsE36 = new HashSet<HashSet<Perm>>(new SetEquality<Perm>()) { g36.ToHashSet() };
        Console.WriteLine(g36.All(g0.Contains));
        DisplayGroup.Head(g36);
        AddToSubGroups(g36, g36, allSubGroupsE36);
        GlobalStopWatch.Show($"All SubGroups of {g0.Name} : {allSubGroups.Count}");

        var allSortedE36 = allSubGroupsE36.GroupBy(sg => sg.Count).ToDictionary(a => a.Key, a => a.ToList());
        Console.WriteLine(allSortedE36.AscendingByKey().Select(p => $"[{p.Key}, {p.Value.Count}]").Glue());
        var g18 = Group.Generate("E18", g0, allSubGroupsE36.First(sg => sg.Count == 18).ToArray());
        DisplayGroup.Head(g18);
        Console.WriteLine(g18.PseudoGenerators.Glue());

        AddToSubGroups(g18, g0, allSubGroups);
        GlobalStopWatch.Show($"All SubGroups of {g0.Name} : {allSubGroups.Count}");
        allSorted = allSubGroups.GroupBy(sg => sg.Count).ToDictionary(a => a.Key, a => a.ToList());
        Console.WriteLine(allSorted.AscendingByKey().Select(p => $"[{p.Key}, {p.Value.Count}]").Glue());
    }

    // # All SubGroups of Symm6 : 1455 Time:127297 ms
    // [1, 1][2, 75][3, 40][4, 255][5, 36][6, 280][8, 255][9, 10][10, 36][12, 150][16, 45][18, 50][20, 36][24, 90][36, 30][48, 30][60, 12][72, 10][120, 12][360, 1][720, 1]
    public static void SubGroupsS6()
    {
        GlobalStopWatch.Restart();
        var gb = new Symm(6);
        var g0 = gb;
        var allSubGroups = new HashSet<HashSet<Perm>>(10000, new SetEquality<Perm>());
        foreach (var e1 in g0)
        {
            foreach (var e2 in g0)
            {
                var g1 = Group.GenerateElements(gb, e1, e2).ToHashSet();
                allSubGroups.Add(g1);
            }
        }

        GlobalStopWatch.Show($"All SubGroups of {g0.Name} : {allSubGroups.Count}");
        var allSorted = allSubGroups.GroupBy(sg => sg.Count).ToDictionary(a => a.Key, a => a.ToList());
        Console.WriteLine(allSorted.AscendingByKey().Select(p => $"[{p.Key}, {p.Value.Count}]").Glue());

        var set36 = allSubGroups.First(s => s.Count == 36);
        var g36 = Group.Generate("E36", g0, set36.ToArray());
        var allSubGroupsE36 = new HashSet<HashSet<Perm>>(new SetEquality<Perm>()) { g36.ToHashSet() };
        Console.WriteLine(g36.All(g0.Contains));
        DisplayGroup.Head(g36);
        AddToSubGroups(g36, g36, allSubGroupsE36);

        var allSortedE36 = allSubGroupsE36.GroupBy(sg => sg.Count).ToDictionary(a => a.Key, a => a.ToList());
        Console.WriteLine(allSortedE36.AscendingByKey().Select(p => $"[{p.Key}, {p.Value.Count}]").Glue());
        var g18 = Group.Generate("E18", g0, allSubGroupsE36.First(sg => sg.Count == 18).ToArray());
        DisplayGroup.Head(g18);
        Console.WriteLine(g18.PseudoGenerators.Glue());
        AddToSubGroups(g18, g0, allSubGroups);
        GlobalStopWatch.Show($"All SubGroups of {g0.Name} : {allSubGroups.Count}");

        var g16 = Group.Generate("E16", g0, g0[(1, 2, 3, 4)], g0[(1, 4), (2, 3)], g0[(5, 6)]);
        AddToSubGroups(g16, g0, allSubGroups);
        GlobalStopWatch.Show($"All SubGroups of {g0.Name} : {allSubGroups.Count}");
        allSorted = allSubGroups.GroupBy(sg => sg.Count).ToDictionary(a => a.Key, a => a.ToList());
        Console.WriteLine(allSorted.AscendingByKey().Select(p => $"[{p.Key}, {p.Value.Count}]").Glue());
    }
}
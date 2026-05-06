using System.Text;
using Craft;
using Examples;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.EllCurve;
using FastGoat.UserGroup.Floats;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Matrix;
using FastGoat.UserGroup.Perms;
using FastGoat.UserGroup.Polynoms;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");
Perm.Style = DisplayPerm.CyclesComma;

IEnumerable<int> ActOnSet(Perm m, IEnumerable<int> l) => l.Select(i => m.Table[i]).ToHashSet();
GroupAction<Perm, XSet<int>> Image = (g, x) => new(ActOnSet(g, x));

IEnumerable<T2> Orbits<T1, T2>(HashSet<T1> gens, GroupAction<T1, T2> act, T2 x, HashSet<T2> set)
    where T1 : struct, IElt<T1>
    where T2 : struct, IElt<T2>
{
    yield return x;
    Queue<T2> q = new Queue<T2>(gens.Count);
    q.Enqueue(x);
    while (q.Count != 0)
    {
        var x0 = q.Dequeue();
        foreach (var g in gens)
        {
            var x1 = act(g, x0);
            if (set.Add(x1))
            {
                q.Enqueue(x1);
                yield return x1;
            }
        }
    }
}

IEnumerable<XSet<int>> AllSets(int n, int k = 1) => IntExt.YieldCombsKinN(k, n)
    .Select(l => l.Zip(n.Range()).Where(e => e.First).Select(e => e.Second))
    .Select(l => new XSet<int>(l));

IEnumerable<XSet<XSet<int>>> SetPartitions(int n)
{
    var set = AllSets(n).ToArray();
    var combs = set.AllCombinations().Where(l => l.Length > 0)
        .Select(l => new XSet<int>(l.SelectMany(e => e))).ToArray();
    var table = combs.GroupBy(e => e.Count).ToDictionary(e => e.Key, e => e.ToArray());

    foreach (var part in IntExt.Partitions32[n])
    {
        var r = new List<XSet<int>>(n);
        var t = 0;
        foreach (var g in part.Order())
        {
            r.Add(table[g].First(l => l.All(i => i >= t)));
            t += g;
        }

        yield return new(r);
    }
}

// {
//     GlobalStopWatch.Restart();
//     var n = 10;
//     var sn = new Sn(n);
//     var fmt = $"{{0,-3}} {{1,{-4 * n}}} type:{{2,{-2 * n - 2}}} count:{{3}}";
//     var classes = new Dictionary<Perm, HashSet<Perm>>();
//     foreach (var (l, c) in SetPartitions(n).Select(l => (l, c: IntExt.Lcm(l.Select(e => e.Count).ToArray())))
//                  .OrderBy(e => e.c))
//     {
//         var perm = sn.CreateElementTable(IntExt.PermAndCyclesFromType(l.Select(e => e.Count).ToArray()).perm);
//         var orbx = Orbits(sn.GetGenerators().ToHashSet(), Group.ByConjugate(sn), perm, new() { perm });
//         Console.WriteLine(fmt, c, perm, $"[{perm.PermType.Glue(",")}]", $"{ orbx.Count()}");
//     }
//
//     GlobalStopWatch.Show();
// }

Perm Cycles(int m)
{
    var dec = IntExt.PrimesDec(m).Select(e => e.Key.Pow(e.Value)).ToArray();
    var perm = IntExt.PermAndCyclesFromType(dec).perm;
    return (new Sn(perm.Length)).CreateElementTable(perm);
}

Perm[] CyclesSplit(int m)
{
    var dec = IntExt.PrimesDec(m).Select(e => e.Key.Pow(e.Value)).ToArray();
    return dec.Select(e => new Sn(e).Cycle(e.Range(1))).ToArray();
}

Dictionary<string, string> Subs(Dictionary<string, string> map, Dictionary<string, string> subs)
{
    var map0 = new Dictionary<string, string>();
    foreach (var (k, v) in map)
    {
        if (subs.ContainsKey(k))
            map0[subs[k]] = subs.ContainsKey(v) ? subs[v] : v;
        else if (subs.ContainsKey(v))
            map0[k] = subs[v];
        else
            map0[k] = v;
    }

    return map0;
}

int Score(KeyValuePair<string, string> e)
{
    var (x, y) = (e.Key.Contains('x') ? 0 : 1, e.Value.Contains('x') ? 0 : 1);
    var bonus = e.Key.Equals(e.Value) ? 0 : 100;
    return 10 * y + x + bonus;
}

var (count, success, all, maxIter) = (0, 0, 0, 0);

IEnumerable<Dictionary<string, string>> SolveAction(Dictionary<string, string> map, Dictionary<string, string> f,
    Dictionary<string, string> subs)
{
    ++count;
    if (subs.Count == map.Count)
    {
        yield return subs;
        yield break;
    }
    
    if (subs.Count < map.Count)
    {
        // Solving the cycles
        var (k, v) = map.OrderByDescending(e => Score(e)).First(e => e.Key.Contains('x') || e.Value.Contains('x'));
        var (kx, vx) = (k.Contains('x'), v.Contains('x'));
        if (kx)
        {
            if (vx)
            {
                var remK = f.Keys.Where(e => !subs.ContainsValue(e)).ToArray();
                foreach (var s in remK)
                {
                    var subs0 = subs.ToDictionary(e => e.Key, e => e.Value);
                    subs0[k] = s;
                    var map0 = Subs(map, subs0);
                    foreach (var sol in SolveAction(map0, f, subs0))
                        yield return sol;
                }
            }
            else
            {
                var subs0 = subs.ToDictionary(e => e.Key, e => e.Value);
                var _k = f.First(e => e.Value == v).Key;
                if (subs0.Values.Contains(_k))
                {
                    var remK = f.Keys.Where(e => !subs.ContainsValue(e)).ToArray();
                    foreach (var s in remK)
                    {
                        var subs1 = subs.ToDictionary(e => e.Key, e => e.Value);
                        subs1[k] = s;
                        var map0 = Subs(map, subs1);
                        foreach (var sol in SolveAction(map0, f, subs1))
                            yield return sol;
                    }
                }
                else
                {
                    subs0[k] = _k;
                    var map0 = Subs(map, subs0);
                    foreach (var sol in SolveAction(map0, f, subs0))
                        yield return sol;
                }
            }
        }
        else
        {
            if (vx)
            {
                var subs0 = subs.ToDictionary(e => e.Key, e => e.Value);
                var _k = f.First(e => e.Key == k).Value;
                subs0[v] = _k;
                var map0 = Subs(map, subs0);
                foreach (var sol in SolveAction(map0, f, subs0))
                    yield return sol;
            }
        }
    }
}

(Perm a, Perm b) FindAB(int m, int r)
{
    var a0 = Cycles(m);
    var cycles = IntExt.PermutationToCycles(a0.Sn.N, a0.Table);
    var seq = new List<Perm>() { a0.Sn.Neutral() };
    foreach (var c in cycles)
    {
        var dic = c.Index().ToDictionary(e => e.Index, e => e.Item);
        var rev = dic.ToDictionary(e => e.Value, e => e.Key);
        var N = c.Length;
        var sn = new Sn(N);
        var a = sn.Cycle(N.Range(1));
        var ar = a ^ r;

        var d_a = a.Table.Index().ToDictionary(e => $"{e.Index}", e => $"{e.Item}");
        var d_ar = ar.Table.Index().ToDictionary(e => $"{e.Index}", e => $"{e.Item}");
        var d_bi = N.SeqLazy(1).Index().ToDictionary(e => $"{e.Index}", e => $"x{e.Item}");
        var d_b = d_bi.ToDictionary(e => e.Value, e => e.Key);
        var d_ba = d_b.ToDictionary(e => e.Key, e => d_a[e.Value]);
        var d_babi = d_ba.ToDictionary(e => e.Key, e => d_bi[e.Value]);

        var sol = SolveAction(d_babi, d_ar, new()).First();
        var d_bf = Subs(d_b, sol);
        var table = N.SeqLazy().ToDictionary(e => e, e => d_bf.ContainsKey($"{e}") ? int.Parse(d_bf[$"{e}"]) : e)
            .OrderBy(e => e.Key)
            .Select(e => e.Value).ToArray();
        var b = sn.CreateElementTable(table);
        
        seq.Add(a0.Sn.CreateElementTable(
                a0.Sn.N.SeqLazy()
                    .Select(i => rev.ContainsKey(i) ? dic[b.Table[rev[i]]] : i)
                    .ToArray()
            )
        );
    }

    return (a0, Group.OpSeq(a0.Sn, seq));
}

IEnumerable<(Perm a, Perm b)> ExtendB(Perm a, Perm b, int n)
{
    var N0 = a.Sn.N;
    foreach (var perms in CyclesSplit(n).AllCombinations().Where(e => e.Length > 0))
    {
        var table = IntExt.PermAndCyclesFromType(perms.SelectMany(e => e.PermType).Order().ToArray()).perm;
        var b0 = new Sn(table.Length).CreateElementTable(table);
        var N = b0.Table.Length + N0;
        var sn = new Sn(N);
        var a1 = sn.CreateElementTable(a.Table.Concat(b0.Table.Length.SeqLazy(N0)).ToArray());
        var b2 = sn.CreateElementTable(b.Table.Concat(b0.Table.Select(i => i + N0).ToArray()).ToArray());
        yield return (a1, b2);
    }
}

ConcreteGroup<Perm> SearchPermGroupMetaCyclic<T>(ConcreteGroup<T> G0, int m, int n, int r)
    where T : struct, IElt<T>
{
    Console.WriteLine($"############ Search {G0.ShortName} ############");
    var (a, b) = FindAB(m, r);
    var sn = a.Sn;
    var ar = a ^ r;

    var ord = Group.ElementIsOrder(sn, b, n);
    var c1 = (b ^ n).Equals(sn.Neutral());
    var c2 = (b * a * (b ^ -1)).Equals(ar);
    if (!ord && c1 && c2)
    {
        foreach (var (a0, b0) in ExtendB(a, b, n))
        {
            var H1 = Group.GenerateElements(a0.Sn, a0, b0);
            if (H1.Count == m * n)
            {
                (sn, a, ar, b) = (a0.Sn, a0, a0 ^ r, b0);
                break;
            }
        }

        ord = Group.ElementIsOrder(sn, b, n);
        c1 = (b ^ n).Equals(sn.Neutral());
        c2 = (b * a * (b ^ -1)).Equals(ar);
    }

    if (ord && c1 && c2)
    {
        var G = Group.Generate(G0.Name, sn, a, b);
        if (G.Count() == n * m)
        {
            Console.WriteLine($"count:{count}");
            Console.WriteLine($"a:{a}");
            Console.WriteLine($"b:{b}");
            Console.WriteLine($"b * a * b^-1 = a^{r} = {ar}");
            Console.WriteLine("############ CANDIDATE ############");
            DisplayGroup.HeadGenerators(G);
            DisplayGroup.AreIsomorphics(G0, G);
            ++success;
            Console.WriteLine();
            return G;
        }
    }

    return Group.Generate("G", sn, sn.Neutral());
}

void SolveMetaCyclic(int m, int n, int r)
{
    ++all;
    count = 0;
    var G = FG.MetaCyclicSdp(m, n, r);
    G.Name = IntExt.Gcd(m, n * (r - 1)) == 1 ? $"F({m}x:{n}){r}" : $"M({m}x:{n}){r}";
    if (SearchPermGroupMetaCyclic(G, m, n, r).Count() == m * n)
        maxIter = int.Max(maxIter, count);
    else
    {
        Console.Beep();
        Console.WriteLine($"count:{count}");
        Console.WriteLine("############ NOT FOUND ############");
        var G1 = G.ToPermGroup();
        DisplayGroup.HeadOrdersGenerators(G1.Item1);
    }
}

IEnumerable<(int m, int n, int r)> MetaCyclic(int ord)
{
    return IntExt.Dividors(ord).Where(d => d > 1)
        .SelectMany(m => FG.MetaCyclicSdpGetR(m, ord / m).Select(r => (m, n: ord / m, r)))
        .Select(e => (e.m, e.n, e.r))
        .Where(e => e.r > 1)
        .OrderBy(e => e);
}

void FirstExamples()
{
    SolveMetaCyclic(4, 4, 3);
    SolveMetaCyclic(3, 8, 2);
    SolveMetaCyclic(5, 10, 4);

    SolveMetaCyclic(12, 4, 5);
    SolveMetaCyclic(12, 4, 7);
    SolveMetaCyclic(24, 2, 5);
    SolveMetaCyclic(5, 12, 2);
    SolveMetaCyclic(5, 12, 4);
}

void MetaCyclicPermFromUpTo128()
{
    (success, all, maxIter) = (0, 0, 0);
    GlobalStopWatch.Restart();
    for (int i = 5; i <= 128; i++)
    {
        foreach (var (m, n, r) in MetaCyclic(i))
            SolveMetaCyclic(m, n, r);
    }

    GlobalStopWatch.Show($"Success:{success}/{all} maxIter:{maxIter}");
}

{
    MetaCyclicPermFromUpTo128(); // # Success:390/390 maxIter:65 Time:1.292s
    // FirstExamples();
    // |F(3x:8)2| = 24
    // Type        NonAbelianGroup
    // BaseGroup   S11
    // 
    // Generators of F(3x:8)2
    // gen1 of order 3
    // [(1, 2, 3)]
    // gen2 of order 8
    // [(2, 3), (4, 5, 6, 7, 8, 9, 10, 11)]
    // 
    // F(3x:8)2 IsIsomorphicTo F(3x:8)2 : True
    // 
    // 
    // |F(5x:12)2| = 60
    // Type        NonAbelianGroup
    // BaseGroup   S8
    // 
    // Generators of F(5x:12)2
    // gen1 of order 5
    // [(1, 2, 3, 4, 5)]
    // gen2 of order 12
    // [(2, 4, 5, 3), (6, 7, 8)]
    // 
    // F(5x:12)2 IsIsomorphicTo F(5x:12)2 : True
    // 
    // 
    // |F(5x:12)4| = 60
    // Type        NonAbelianGroup
    // BaseGroup   S12
    // 
    // Generators of F(5x:12)4
    // gen1 of order 5
    // [(1, 2, 3, 4, 5)]
    // gen2 of order 12
    // [(2, 5), (3, 4), (6, 7, 8), (9, 10, 11, 12)]
    // 
    // F(5x:12)4 IsIsomorphicTo F(5x:12)4 : True
    // 
    // 
}
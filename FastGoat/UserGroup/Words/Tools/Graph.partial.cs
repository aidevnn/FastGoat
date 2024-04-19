using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.GenericGroup;

namespace FastGoat.UserGroup.Words.Tools;

public partial class Graph
{
    private Graph(Gen[] gens, int capacity = 10000)
    {
        Gens = gens;
        NbGens = gens.Length;
        End = false;
        Relators = new();
        Subgroup = new Gen[0];
        var Null = new Class(0, this);
        Classes = new(capacity) { Null };
    }

    private bool STDone { get; set; }
    private bool RelatorsComplete { get; set; } = true;

    private void SpanningTree()
    {
        var clI = Classes[1];
        var clL = Classes[1];
        var n = Classes.Count - 1;
        var k = 1;
        var gi = new Gen();
        var gens = Gens.OrderByDescending(c => char.IsLower(c.V)).ToArray();
        var W = Classes.Select(c0 => (c0, c1: Class.Null, c2: Class.Null, g: gi)).ToArray();
        W[1].c2 = clI;
        while (true)
        {
            foreach (var gJ in gens)
            {
                var clT = clI![gJ];
                if (W[clT.V].c1 is not null || clT.V == clL.V)
                    continue;

                ++k;
                W[clT.V].c2 = clI;
                W[clT.V].g = gJ;
                W[clL.V].c1 = clT;
                clL = clT;
                if (k == n)
                    break;
            }

            if (k == n)
                break;

            clI = W[clI!.V].c1;
        }

        foreach (var (c0, _, c2, g) in W.Skip(1))
        {
            if (!g.Equals(gi))
            {
                c0.STClass = c2!;
                c0.STGen = g;
                c2!.Color(g);
            }
        }
    }

    private void GenerateWords()
    {
        foreach (var @class in Classes.Skip(1))
        {
            var tmp = @class;
            while (tmp.V != 1)
            {
                var g = tmp.STGen;
                tmp = tmp.STClass!;
                @class.Word.Insert(0, g);
            }

            @class.WordInv.AddRange(@class.Word.Select(g => g.Invert()).Reverse());
            if (Logger.Level != LogLevel.Off)
                if (Classes.Count < 33)
                    Console.WriteLine($"({@class}) = ({@class.Word.Glue(" * ")})");
        }
    }

    private void NextRelator()
    {
        var gensPos = Gens.OrderByDescending(g => char.IsLower(g.V))
            .ThenBy(g => g)
            .Select((g, i) => (g, i))
            .ToDictionary(e => e.g, e => e.i);
        var comp = Comparer<IEnumerable<Gen>>.Create((a, b) => a.Select(e => gensPos[e]).SequenceCompareTo(b.Select(e => gensPos[e])));
        var rem = Classes.SelectMany(e => e.GetEdges().Select(f => (i: e, s: f.Key, j: f.Value)))
            .Where(e => !e.i.Coloured[e.s])
            .Select(e => (e, e.i.Word.Append(e.s).Concat(e.j.WordInv).ToArray()))
            .OrderBy(e => e.Item2.Length)
            .ThenBy(e => e.Item2.Select(c => char.ToLower(c.V)).Distinct().Count())
            .ThenByDescending(e => char.IsLower(e.e.s.V))
            .ThenBy(e => e.Item2, comp)
            .ToArray();

        if (rem.Length == 0)
            return;

        var (e0, rel) = rem[0];

        e0.i.Color(e0.s);
        var relator = new Relator(Relators.Count, rel);
        Relators.Add(relator);
        foreach (var @class in Classes.Skip(1))
        {
            var circuit = new Circuit(@class, relator);
            var (_, err) = circuit.UpdateCircuit();
            if (err is not null)
                throw new();

            @class.Circuits.Add(circuit);
        }

        StepColorAll();
    }

    private void StepColorAll()
    {
        var eq = EqualityComparer<(Class, Gen)>.Create(
            (a, b) => (a.Item1.V, b.Item2.V).Equals((b.Item1.V, b.Item2.V)),
            obj => (obj.Item1.V, obj.Item2.V).GetHashCode());
        while (true)
        {
            var stop = true;
            var circuits = Classes.SelectMany(c => c.Circuits).ToArray();
            var edges = circuits.Select(c => c.NotColoured)
                .Where(c => c.Count == 1)
                .SelectMany(c => c)
                .Distinct(eq)
                .ToArray();

            foreach (var (@class, gen) in edges)
            {
                stop = false;
                @class.Color(gen);
            }

            if (stop)
                break;
        }
    }

    private Dictionary<Gen, Gen> GensConv { get; set; } = new();

    private void DefineRelators()
    {
        if (!STDone)
        {
            Relators.Clear();
            RelatorsComplete = false;

            SpanningTree();
            GenerateWords();

            while (true)
            {
                var sz = Relators.Count;
                NextRelator();
                if (Relators.Count == sz)
                    break;

                if (Logger.Level != LogLevel.Off)
                {
                    var rel = Relators.Last();
                    Console.WriteLine($"Rel[{Relators.Count}]:{StringExt.ReducedWordForm1(rel.Gens.Glue())} Length:{rel.Length}");
                }
            }

            STDone = true;
        }

        var gensConv = GensConv;
        if (Logger.Level != LogLevel.Off)
            Relators.Select(c => c.Gens.Select(g => gensConv[g]))
                .Select(c => c.Glue())
                .Select(StringExt.ReducedWordForm1)
                .OrderBy(s => s.Length)
                .ThenBy(s => s)
                .Println("All Relators");
    }

    public IEnumerable<char> Rewrite(IEnumerable<char> w)
    {
        var one = Classes[1];
        var cl = w.Aggregate(one, (acc, c) => acc[new(c)]);
        return cl.Word.Select(g => g.V);
    }

    public IEnumerable<char> Generators() => Gens.Select(c => c.V).Where(c => char.IsLower(c));
    public IEnumerable<IEnumerable<char>> Words() => Classes.Skip(1).Select(c => c.Word.Select(e => e.V));

    public bool CheckHomomorphism<T>(ConcreteGroup<T> g, Dictionary<char, T> map) where T : struct, IElt<T>
    {
        foreach (var rel in Relators)
        {
            var r = g.OpSeq(rel.Gens.Select(e => char.IsLower(e.V) ? map[e.V] : g.Invert(map[char.ToLower(e.V)])));
            if (!r.Equals(g.Neutral()))
                return false;
        }

        return true;
    }

    public static Graph Run(string sg, string rel)
    {
        var g = new Graph(sg, rel);
        g.Build();
        g.SpanningTree();
        g.GenerateWords();
        return g;
    }

    public static Graph Run(string rel) => Run("i", rel);

    private static Graph ClassesFromGroupSubgroup<T>(ConcreteGroup<T> g, ConcreteGroup<T> h) where T : struct, IElt<T>
    {
        var E = g.GetGenerators().OrderByDescending(e => g.ElementsOrders[e]).ToArray();
        var r = E.Length;
        var cosets = Group.Cosets(g, h, CosetType.Right);
        var repr = cosets.ToDictionary(a => a.Key, a => a.Value.X);

        var gens = r.Range().Select(i => (char)('a' + i)).Select(i => new Gen(i)).Prepend(new Gen()).ToArray();
        var graph = new Graph(gens.Skip(1).SelectMany(c => new[] { c, c.Invert() }).ToArray(), g.Count() * 2);

        graph.GensConv = gens.Skip(1).Select((e, i) => (e, i))
            .SelectMany(e => new[]
            {
                (e.e, e.e),
                g.ElementsOrders[E[e.i]] == 2
                    ? (e.e.Invert(), e.e)
                    : (e.e.Invert(), e.e.Invert())
            })
            .ToDictionary(e => e.Item1, e => e.Item2);

        var F = new List<T>() { g.Neutral() };
        graph.Classes.Add(new(1, graph));
        var m = 1;
        graph.Classes[m].Coloured = graph.Gens.ToDictionary(e => e, _ => false);
        var idx = new Dictionary<T, int>(g.Count());
        idx[g.Neutral()] = 1;
        for (int i = 1; i <= m; ++i)
        {
            for (int j = 1; j <= r; ++j)
            {
                var x = repr[g.Op(F[i - 1], E[j - 1])];
                if (idx.TryGetValue(x, out var l))
                    graph.Classes[i][gens[j]] = graph.Classes[l];
                else
                {
                    F.Add(x);
                    var cl = new Class(++m, graph);
                    idx[x] = m;
                    cl.Coloured = graph.Gens.ToDictionary(e => e, _ => false);
                    graph.Classes.Add(cl);
                    graph.Classes[i][gens[j]] = cl;
                }
            }
        }

        idx.Clear();
        F.Clear();
        return graph;
    }

    public static string DefiningRelatorsOfGroup<T>(ConcreteGroup<T> g) where T : struct, IElt<T>
    {
        var h = Group.Generate("()", g, g.Neutral());

        if (Logger.Level != LogLevel.Off)
            GlobalStopWatch.AddLap();

        var graph = ClassesFromGroupSubgroup(g, h);

        if (Logger.Level != LogLevel.Off)
            Console.WriteLine(g.ShortName);

        if (Logger.Level != LogLevel.Off)
            if (g.Count() < 33)
            {
                graph.DisplayTableOps();
                Console.WriteLine();
            }

        graph.DefineRelators();

        if (Logger.Level != LogLevel.Off)
            GlobalStopWatch.Show();

        return graph.Relators.Select(c => c.Gens.Select(g0 => graph.GensConv[g0]))
            .Select(c => c.Glue())
            .Select(StringExt.ReducedWordForm1)
            .OrderBy(s => s.Length)
            .ThenBy(s => s)
            .Glue(",");
    }
}
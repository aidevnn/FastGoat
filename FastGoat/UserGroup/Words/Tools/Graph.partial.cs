using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.GenericGroup;

namespace FastGoat.UserGroup.Words.Tools;

public partial class Graph
{
    private Graph(Gen[] gens)
    {
        Gens = gens;
        NbGens = gens.Length;
        End = false;
        Relators = new();
        Subgroup = new Gen[0];
        var Null = new Class(0, this);
        Classes = new() { Null };
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

    private void GenerateWords(bool details = false)
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

            if (details && Classes.Count < 33)
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

    private void DefineRelators(bool details = true)
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
            }

            STDone = true;
        }

        var gensConv = GensConv;
        if (details)
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

    public static Graph Run(string sg, string rel)
    {
        var g = new Graph(sg, rel);
        g.Build(false, false);
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
        var graph = new Graph(gens.Skip(1).SelectMany(c => new[] { c, c.Invert() }).ToArray());

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
        for (int i = 1; i <= m; ++i)
        {
            for (int j = 1; j <= r; ++j)
            {
                var x = repr[g.Op(F[i - 1], E[j - 1])];
                var l = F.FindIndex(e => e.Equals(x)) + 1;
                if (l != 0)
                    graph.Classes[i][gens[j]] = graph.Classes[l];
                else
                {
                    F.Add(x);
                    var cl = new Class(++m, graph);
                    cl.Coloured = graph.Gens.ToDictionary(e => e, _ => false);
                    graph.Classes.Add(cl);
                    graph.Classes[i][gens[j]] = cl;
                }
            }
        }

        return graph;
    }

    public static string DefiningRelatorsOfGroup<T>(ConcreteGroup<T> g, bool details = true) where T : struct, IElt<T>
    {
        var h = Group.Generate("()", g, g.Neutral());

        if (details)
            GlobalStopWatch.AddLap();
        
        var graph = ClassesFromGroupSubgroup(g, h);

        if (details)
            Console.WriteLine(g.ShortName);

        if (details && g.Count() < 33)
        {
            graph.DisplayTableOps();
            Console.WriteLine();
        }

        graph.DefineRelators(details);
        
        if (details)
            GlobalStopWatch.Show();
        
        return graph.Relators.Select(c => c.Gens.Select(g0 => graph.GensConv[g0]))
            .Select(c => c.Glue())
            .Select(StringExt.ReducedWordForm1)
            .OrderBy(s => s.Length)
            .ThenBy(s => s)
            .Glue(",");
    }
}
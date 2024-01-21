using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.GenericGroup;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Words;

namespace FastGoat.Examples;

// DISCRETE MATHEMATICS 5 (1973) 105 -129. North-Molland Publishing Company
// CONSTRUCTION OF DEFINING RELATORS
// FOR FINITE GROUPS
// John J. CANNON
public static class DefiningRelators
{
    readonly struct Table<T> where T : struct, IEquatable<T>
    {
        private T[,] table { get; }

        public Table(int m, int n)
        {
            table = new T[m, n];
        }

        private Table(Table<T> table0)
        {
            table = new T[table0.W, table0.H];
            Array.Copy(table0.table, table, table0.W * table0.H);
        }

        public Table<T> Clone => new(this);

        public IEnumerable<(int I, int S, T T)> GetContent()
        {
            var table0 = table;
            return Enumerable.Range(1, W).Grid2D(Enumerable.Range(1, H))
                .Select(e => (e.t1, e.t2, table0[e.t1 - 1, e.t2 - 1]));
        }

        public int W => table.GetLength(0);
        public int H => table.GetLength(1);

        public T this[int i, int j]
        {
            get => table[i - 1, j - 1];
            set => table[i - 1, j - 1] = value;
        }

        public override string ToString()
        {
            return Ring.Matrix2String(table);
        }
    }

    static (T1[] E, List<T1> F, Dictionary<T1, T1> repr, Table<int> T) EdgesTable<T1>(ConcreteGroup<T1> g, ConcreteGroup<T1> h)
        where T1 : struct, IElt<T1>
    {
        var E = g.GetGenerators().OrderBy(e => g.ElementsOrders[e]).ToArray();
        var r = E.Length;
        var hi = g.Count() / h.Count();
        var cosets = Group.Cosets(g, h, CosetType.Right);
        var repr = cosets.ToDictionary(a => a.Key, a => a.Value.X);
        var T = new Table<int>(hi, 2 * r);

        var F = new List<T1>() { g.Neutral() };
        var m = 1;
        for (int i = 1; i <= m; ++i)
        {
            for (int j = 1; j <= r; ++j)
            {
                var x = repr[g.Op(F[i - 1], E[j - 1])];
                var l = F.FindIndex(e => e.Equals(x)) + 1;
                if (l != 0)
                {
                    T[i, j] = l;
                    T[l, r + j] = i;
                }
                else
                {
                    F.Add(x);
                    ++m;
                    T[i, j] = m;
                    T[m, r + j] = i;
                }
            }
        }

        return (E, F, repr, T);
    }

    static Table<int> SpanningTree(Table<int> T)
    {
        var (n, s) = (T.W, T.H);
        var W = new Table<int>(n, 3);

        W[1, 2] = 1;
        int i = 1, k = 1, l = 1;
        while (true)
        {
            for (int j = 1; j <= s; j++)
            {
                var Tij = T[i, j];
                if (W[Tij, 1] != 0 || Tij == l)
                    continue;

                ++k;
                W[Tij, 2] = i;
                W[Tij, 3] = j;
                W[l, 1] = Tij;
                l = Tij;

                if (k == n)
                    return W;
            }

            if (k == n)
                break;

            i = W[i, 1];
        }

        return W;
    }

    static (Dictionary<int, LinkedList<int>> words, HashSet<(int I, int S, int J)> allWedges) Words(int r, Table<int> W)
    {
        var rgJ = r.Range(1);
        var inv = rgJ.SelectMany(k => new[] { (k, r + k), (r + k, k) }).ToDictionary(e => e.Item1, e => e.Item2);
        var words = new Dictionary<int, LinkedList<int>>() { [1] = new() };

        var allWedges = new HashSet<(int I, int S, int J)>();
        for (int p = 1; p <= W.W; p++)
        {
            if (W[p, 3] == 0)
                continue;

            var q = W[p, 2];
            var sj = W[p, 3];
            var sji = inv[sj];
            allWedges.Add((q, sj, p));
            allWedges.Add((p, sji, q));
        }

        var tmpWedges = allWedges.ToList();
        while (tmpWedges.Count != 0)
        {
            var (i, s, j) = tmpWedges.Where(e => words.ContainsKey(e.I)).MinBy(e => e.I);
            var l0 = new LinkedList<int>(words[i]);
            if (l0.Count == 0 || l0.Last!.Value != inv[s])
                l0.AddLast(s);
            else
                l0.RemoveLast();

            words[j] = l0;
            tmpWedges.Remove((i, s, j));
        }

        return (words, allWedges);
    }

    static void StepColorAll(int r, List<LinkedList<int>> rels, Table<int> edges, Table<bool> coloured)
    {
        var inv = r.Range(1).SelectMany(i => new[] { (i, r + i), (r + i, i) }).ToDictionary(e => e.Item1, e => e.Item2);
        var circuits = new List<List<(int I, int S, int J)>>();
        while (true)
        {
            var stop = true;
            foreach (var rel in rels)
            {
                circuits.Clear();
                var start = true;
                foreach (var s in rel)
                {
                    if (start)
                    {
                        start = false;
                        circuits = edges.W.Range(1).Select(i => new List<(int I, int S, int J)>() { (i, s, edges[i, s]) }).ToList();
                        continue;
                    }

                    var tmp = circuits.ToList();
                    circuits.Clear();
                    foreach (var circuit0 in tmp)
                    {
                        var j0 = circuit0.Last().J;
                        var circuit1 = circuit0.Append((I: j0, S: s, J: edges[j0, s])).ToList();
                        if (circuit1.Count(e => !coloured[e.I, e.S]) <= 1)
                            circuits.Add(circuit1);
                    }
                }

                foreach (var circuit in circuits)
                {
                    var edge = circuit.Where(e => !coloured[e.I, e.S]).ToArray();
                    if (edge.Length == 1)
                    {
                        var e = edge[0];
                        coloured[e.I, e.S] = coloured[e.J, inv[e.S]] = true;
                        stop = false;
                    }
                }
            }

            if (stop)
                break;
        }
    }

    static void NextRelator(int r, List<LinkedList<int>> rels, Table<int> edges, Table<bool> coloured,
        Dictionary<int, LinkedList<int>> words)
    {
        var comp = Comparer<IEnumerable<int>>.Create((a, b) => a.SequenceCompareTo(b));
        var inv = r.Range(1).SelectMany(i => new[] { (i, r + i), (r + i, i) }).ToDictionary(e => e.Item1, e => e.Item2);
        var rem = edges.GetContent().Where(e => !coloured[e.I, e.S] && e.S <= r)
            .Select(e => (e, new LinkedList<int>(words[e.I].Append(e.S).Concat(words[e.T].Reverse().Select(i => inv[i])))))
            .OrderBy(e => e.Item2.Count)
            .ThenBy(e => e.Item2, comp)
            .ToArray();

        if (rem.Length == 0)
            return;

        var (e0, rel) = rem[0];
        coloured[e0.I, e0.S] = coloured[e0.T, inv[e0.S]] = true;
        rels.Insert(0, rel);
        StepColorAll(r, rels, edges, coloured);
    }

    static char[][] FindRelators<T1>(ConcreteGroup<T1> g, bool details = false)
        where T1 : struct, IElt<T1>
    {
        var h = Group.Generate("()", g, g.Neutral());
        var (E, F, repr, T) = EdgesTable(g, h);
        var Ei = E.Concat(E.Select(e => g.Invert(e))).ToArray();
        var r = E.Length;
        var rgJ = r.Range(1);
        var tmp = rgJ.Select(k => (char)('a' + k - 1)).ToArray();
        var gens = tmp.Concat(tmp.Select((a, i) => g.ElementsOrders[Ei[i]] == 2 ? a : char.ToUpper(a))).Select((a, i) => (a, i + 1))
            .ToDictionary(e => e.Item2, e => e.a);

        var W = SpanningTree(T);
        var (words, allWedges) = Words(r, W);

        if (details)
        {
            E.Select((c, i) => ((tmp[i], c), i + 1))
                .ToDictionary(e => e.Item2, e => e.Item1)
                .OrderBy(e => e.Key)
                .Select(e => $"{e.Key}. {e.Value.Item1} -> {e.Value.c}")
                .Println("Generators");

            if (g.Count() < 33)
            {
                Console.WriteLine("Words");
                var table = words.OrderBy(e => e.Key)
                    .SelectMany(e => new[]
                        { $"(1)*({e.Value.Select(k => gens[k]).Glue("*")})", "=", $"({e.Key})", "=", $"{F[e.Key - 1]}" })
                    .ToArray();
                var fmt = Ring.MatrixDisplayForm;
                Ring.MatrixDisplayForm = Ring.MatrixDisplay.TableLeft;
                Ring.DisplayMatrix(Ring.Matrix(words.Count, table), "  ");
                Ring.MatrixDisplayForm = fmt;
                Console.WriteLine();
            }
        }

        var coloured = new Table<bool>(T.W, T.H);
        foreach (var e in allWedges)
            coloured[e.I, e.S] = true;

        var nbEdges = T.W * T.H;
        if (details)
            Console.WriteLine("Coloured:{0}/{1}", coloured.GetContent().Count(e => e.T), nbEdges);
        var rels = new List<LinkedList<int>>();
        while (true)
        {
            NextRelator(r, rels, T, coloured, words);
            var nbColoured = coloured.GetContent().Count(e => e.T);
            if (details)
            {
                Console.WriteLine("Coloured:{0}/{1}", nbColoured, nbEdges);
                Console.WriteLine($"Relator:{rels[0].Select(i => gens[i]).Glue("*")}");
                Console.WriteLine($"ToddCoxeterVerification:{ToddCoxeterVerification(rels, T)}");
            }
            if (nbColoured == nbEdges)
                break;
        }

        return rels.Reverse<LinkedList<int>>().Select(rel => rel.Select(i => gens[i]).ToArray()).ToArray();
    }

    public static void ShowRelators<T>(ConcreteGroup<T> g) where T : struct, IElt<T>
    {
        GlobalStopWatch.Restart();
        Console.WriteLine($"G:{g.ShortName}");
        var rels = FindRelators(g, details: true);
        var strRels = rels.Select(rel => StringExt.ReducedWordForm1(rel.Glue())).OrderBy(s => s.Length).ThenAscending().ToArray();
        var gens = rels.SelectMany(rel => rel).Select(c => char.ToLower(c)).Distinct().Order().ToArray();
        Console.WriteLine($"G:{g.ShortName}");
        Console.WriteLine($"G:<({gens.Glue(",")}) | {strRels.Glue(", ")}>");
        GlobalStopWatch.Show();
        Console.WriteLine();
    }

    static bool ToddCoxeterVerification(List<LinkedList<int>> rels, Table<int> edges)
    {
        var cls = edges.W.Range(1);
        var allRels = rels.SelectMany(rel => rel).ToArray();
        var r = edges.H / 2;
        var coloured = new Table<bool>(edges.W, edges.H);
        var inv = r.Range(1).SelectMany(i => new[] { (i, r + i), (r + i, i) }).ToDictionary(e => e.Item1, e => e.Item2);
        var nb = 0;
        foreach (var i0 in cls)
        {
            var i = i0;
            foreach (var s in allRels)
            {
                var j = edges[i, s];
                if (!coloured[i, s])
                {
                    ++nb;
                    coloured[i, s] = true;
                }
                
                i = j;
            }
        }

        return nb == edges.W * edges.H;
    }

    public static void Example1Abelian()
    {
        ShowRelators(FG.Abelian(2, 2));
        ShowRelators(FG.Abelian(2, 4));
        ShowRelators(FG.Abelian(2, 3));
        ShowRelators(FG.Abelian(2, 6));
        ShowRelators(FG.Abelian(3, 6));
        ShowRelators(FG.Abelian(2, 4, 8));
    }

    public static void Example2Dihedral()
    {
        ShowRelators(FG.Dihedral(3));
        ShowRelators(FG.Dihedral(4));
        ShowRelators(FG.Dihedral(5));
        ShowRelators(FG.Dihedral(6));
        ShowRelators(FG.Dihedral(7));
    }

    public static void Example3Dicyclic()
    {
        Ring.MatrixDisplayForm = Ring.MatrixDisplay.SquareBracket;
        ShowRelators(FG.DiCyclic(2));
        ShowRelators(FG.DiCyclic(3));
        ShowRelators(FG.DiCyclic(4));
        ShowRelators(FG.DiCyclic(5));
        ShowRelators(FG.DiCyclic(6));

        ShowRelators(FG.DiCyclicSdp(2));
        ShowRelators(FG.DiCyclicSdp(3));
        ShowRelators(FG.DiCyclicSdp(4));
        ShowRelators(FG.DiCyclicSdp(5));
        ShowRelators(FG.DiCyclicSdp(6)); // TODO fix it
    }

    public static void Example4PermGroup()
    {
        ShowRelators(FG.Alternate(3));
        ShowRelators(FG.Symmetric(3));
        ShowRelators(FG.Alternate(4));
        ShowRelators(FG.Symmetric(4));
        ShowRelators(FG.Alternate(5));
        ShowRelators(FG.Symmetric(5));
        ShowRelators(FG.Alternate(6));
        ShowRelators(FG.Symmetric(6));
        ShowRelators(FG.Alternate(7)); // Time:1.836s
        ShowRelators(FG.Symmetric(7)); // Time:1.769s
    }

    public static void Example5MetaCyclicSdp()
    {
        var nb = 33.Range(1).Except(IntExt.Primes10000).ToArray();
        foreach (var g in nb.SelectMany(k => FG.MetaCyclicSdp(k)))
            ShowRelators(g);
    }
}
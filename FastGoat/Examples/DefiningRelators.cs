using System.Runtime;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.GenericGroup;
using FastGoat.UserGroup;

namespace FastGoat.Examples;

// DISCRETE MATHEMATICS 5 (1973) 105 -129. North-Molland Publishing Company
// CONSTRUCTION OF DEFINING RELATORS
// FOR FINITE GROUPS
// John J. CANNON
public static class DefiningRelators
{
    readonly struct Edge : IElt<Edge>
    {
        public Edge(int s, int i, int j)
        {
            (S, I, J) = (s, i, j);
            Hash = (S, I, J).GetHashCode();
        }

        public int S { get; }
        public int I { get; }
        public int J { get; }
        public int Hash { get; }

        public void Deconstruct(out int s0, out int i0, out int j0)
        {
            (s0, i0, j0) = (S, I, J);
        }

        public bool Equals(Edge other) => S == other.S && I == other.I && J == other.J;
        public int CompareTo(Edge other) => (I, J, S).CompareTo((other.I, other.J, other.S));
        public override int GetHashCode() => Hash;
        public override string ToString() => $"{(S, I, J)}";
    }

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

        public IEnumerable<(int i, int j, T t)> GetContent()
        {
            var table0 = table;
            var d = default(T);
            return Enumerable.Range(1, W).Grid2D(Enumerable.Range(1, H))
                .Select(e => (e.t1, e.t2, table0[e.t1 - 1, e.t2 - 1]))
                .Where(e => !e.Item3.Equals(d));
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

    static char[][] FindRelators<T1>(ConcreteGroup<T1> g, bool details = false)
        where T1 : struct, IElt<T1>
    {
        var h = Group.Generate("()", g, g.Neutral());
        var (E, F, repr, T) = EdgesTable(g, h);
        var Ei = E.Concat(E.Select(e => g.Invert(e))).ToArray();
        var (r, rows, cols) = (E.Length, T.W, T.H);
        var rgJ = r.Range(1);
        var inv = rgJ.SelectMany(k => new[] { (k, r + k), (r + k, k) }).ToDictionary(e => e.Item1, e => e.Item2);
        var W = SpanningTree(T);

        var tmp = rgJ.Select(k => (char)('a' + k - 1)).ToArray();
        var gens = tmp.Concat(tmp.Select((a, i) => g.ElementsOrders[Ei[i]] == 2 ? a : char.ToUpper(a))).Select((a, i) => (a, i + 1))
            .ToDictionary(e => e.Item2, e => e.a);

        var allTedges = new HashSet<Edge>();
        for (int i0 = 1; i0 <= rows; i0++)
        {
            var i = i0;
            allTedges.UnionWith(rgJ.SelectMany(si => new[]
            {
                new Edge(si, i, T[i, si]), 
                new Edge(r + si, i, T[i, r + si])
            }));
        }

        var allWedges = new HashSet<Edge>();
        for (int p = 1; p <= rows; p++)
        {
            if (W[p, 3] == 0)
                continue;

            var q = W[p, 2];
            var sj = W[p, 3];
            var sji = inv[sj];
            allWedges.UnionWith(new[] { new Edge(sj, q, p), new Edge(sji, p, q) });
        }

        var tmpWedges = allWedges.ToList();
        var words = new Dictionary<int, LinkedList<int>>() { [1] = new() };
        while (tmpWedges.Count != 0)
        {
            var (s, i, j) = tmpWedges.Where(e => words.ContainsKey(e.I)).MinBy(e => e.I);
            var l0 = new LinkedList<int>(words[i]);
            if (l0.Count == 0 || l0.Last!.Value != inv[s])
                l0.AddLast(s);
            else
                l0.RemoveLast();

            if (words.TryGetValue(j, out var l1))
            {
                if (!l1.SequenceEqual(l0))
                    throw new($"({i})*{gens[s]} = ({j}); [{l0.Glue(" ")}] != [{l1.Glue(" ")}]");
            }
            else
                words[j] = l0;

            tmpWedges.Remove(new Edge(s, i, j));
        }

        if (details)
        {
            E.Select((c, i) => ((tmp[i], c), i + 1))
                .ToDictionary(e => e.Item2, e => e.Item1)
                .OrderBy(e => e.Key)
                .Select(e => $"{e.Key}. {e.Value.Item1} -> {e.Value.c}")
                .Println("Generators");
            
            foreach (var (s, i, j) in allWedges.Union(allTedges))
            {
                var (a, e1, e2) = (Ei[s - 1], F[i - 1], F[j - 1]);
                if (!repr[g.Op(e1, a)].Equals(e2))
                    throw new();
            }

            if (words.Count <= 32)
            {
                var nb = words.Count;
                var table = words.OrderBy(e => e.Key)
                    .SelectMany(e => new[]
                        { $"(1)*({e.Value.Select(k => gens[k]).Glue("*")})", "=", $"({e.Key})", "=", $"{F[e.Key - 1]}" })
                    .ToArray();

                Console.WriteLine("Words");
                var fmt = Ring.MatrixDisplayForm;
                Ring.MatrixDisplayForm = Ring.MatrixDisplay.TableLeft;
                Ring.DisplayMatrix(Ring.Matrix(nb, table), "  ");
                Ring.MatrixDisplayForm = fmt;
            }

            foreach (var (k, op) in words)
            {
                var e1 = F[k - 1];
                var e2 = op.Aggregate(g.Neutral(), (acc, i) => repr[g.Op(acc, Ei[i - 1])]);
                if (!e1.Equals(e2))
                    throw new($"[(1)*({op.Select(k0 => gens[k0]).Glue("*")}) = {e2}] != [({k}) = {e1}]");
            }

            if (words.Count != F.Count)
                throw new($"ops:{words.Count} != F:{F.Count}");
        }

        var coloured = allWedges.ToHashSet();
        var rels = new List<LinkedList<int>>();
        while (true)
        {
            var (rels0, newcols) = AddRelator(r, rels, allTedges, coloured, words);
            if (rels0.Count == rels.Count) break;
            rels = rels0;
            coloured.UnionWith(newcols);
        }

        return rels.Reverse<LinkedList<int>>().Select(rel => rel.Select(i => gens[i]).ToArray()).ToArray();
    }
    
    static (List<LinkedList<int>> rels, HashSet<Edge> col) AddRelator(int r, List<LinkedList<int>> rels, HashSet<Edge> edges, 
        HashSet<Edge> coloured, Dictionary<int, LinkedList<int>> words)
    {
        var comp = Comparer<IEnumerable<int>>.Create((a, b) => a.SequenceCompareTo(b));
        var inv = r.Range(1).SelectMany(i => new[] { (i, r + i), (r + i, i) }).ToDictionary(e => e.Item1, e => e.Item2);
        var rem = edges.Except(coloured)
            .Where(e => e.S <= r) 
            .Select(e => (e, new LinkedList<int>(words[e.I].Append(e.S).Concat(words[e.J].Reverse().Select(i => inv[i])))))
            .OrderBy(e => e.Item2.Count)
            .ThenBy(e => e.Item2, comp)
            .ToArray();

        if (rem.Length == 0)
            return (rels, coloured);

        var (es0, rel) = rem[0];
        var es1 = new Edge(inv[es0.S], es0.J, es0.I);
        var circuits = new List<List<Edge>>();
        var cols = coloured.ToHashSet();
        cols.UnionWith(new[] { es0, es1 });
        var rels0 = rels.Prepend(rel).ToList();
        while (true)
        {
            var sz = cols.Count;
            foreach (var rel0 in rels0)
            {
                circuits.Clear();
                var start = false;
                foreach (var s0 in rel0)
                {
                    if (!start)
                    {
                        circuits = edges.Where(e => e.S == s0).Select(e => new List<Edge>() { e }).ToList();
                        start = true;
                        continue;
                    }

                    var tmp = circuits.ToList();
                    circuits.Clear();
                    foreach (var circuit in tmp)
                    {
                        var j0 = circuit.Last().J;
                        foreach (var e1 in edges.Where(e => e.S == s0 && e.I == j0))
                        {
                            var nc = circuit.Append(e1).ToList();
                            circuits.Add(nc);
                        }
                    }
                }

                var final = circuits.Where(c => c.Count(e => !cols.Contains(e)) == 1).ToList();
                cols.UnionWith(final.SelectMany(e => e.SelectMany(f => new[] { f, new Edge(inv[f.S], f.J, f.I) })));
            }

            if (cols.Count == sz)
                break;
        }

        return (rels0, cols);
    }

    static void ShowRelators<T>(ConcreteGroup<T> g) where T : struct, IElt<T>
    {
        Console.WriteLine($"G:{g.ShortName}");
        var rels = FindRelators(g, details: true);
        rels.Select(rel => $"{rel.Glue("*")} = 1").Println("Relators");
        Console.WriteLine();
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
    
    public static void Example2Dihedrals()
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
    }
}
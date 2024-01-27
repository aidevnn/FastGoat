using FastGoat.Commons;

namespace FastGoat.UserGroup.Words.Tools;

public partial class Class
{
    public int V { get; set; }
    public Graph G { get; }

    public bool IsComplete { get; set; }
    public List<Circuit> Circuits { get; }

    public Class(int v, Graph g)
    {
        V = v;
        G = g;
        Edges = new(G.NbGens * 2);
        Circuits = G.Relators.Select(rel => new Circuit(this, rel)).ToList();
        if (v == 0)
        {
            IsComplete = true;
            Circuits.Clear();
        }
    }
    private Dictionary<Gen, Class> Edges { get; }
    public bool HasEdge(Gen g) => Edges.ContainsKey(g);
    public IEnumerable<KeyValuePair<Gen, Class>> GetEdges() => Edges;

    public Class this[Gen index]
    {
        get => Edges[index];
        set
        {
            var cl = Edges[index] = value;
            cl.Edges.TryAdd(index.Invert(), this);
        }
    }

    public (Class? cl1, Class? cl2) UpdateClass()
    {
        if (IsComplete)
            return (null, null);

        IsComplete = false;
        foreach (var circuit in Circuits)
        {
            var e = circuit.UpdateCircuit();
            if (e.cl1 is not null)
                return e;
        }

        IsComplete = Circuits.All(c => c.IsComplete);
        return (null, null);
    }
    
    public void Substitute(Class cl, Class err)
    {
        foreach (var (g, cl0) in Edges.ToArray())
        {
            if (cl0.V == err.V)
                Edges[g] = cl;
        }

        foreach (var circuit in Circuits)
        {
            circuit.Substitute(cl, err);
        }
    }
    
    public (Class? cl, Gen g) FindCandidate()
    {
        return Circuits.Where(c => !c.IsComplete)
            .Select(c => c.FindCandidate())
            .FirstOrDefault(e => e.Item1 is not null, (null, new()));
    }
    
    public string Display(int digits)
    {
        var fmt = $"{{0,{digits + 1}}}";
        var s0 = Circuits.SelectMany(c => c.Content.SkipLast(1)).Append(this).Glue("", fmt);
        var s1 = (s0 + " ").ToArray();
        foreach (var k in G.Separators)
            s1[(digits + 1) * k] = s1[(digits + 1) * (k + 1)] = '│';

        var s2 = s1.Glue();
        return s2;
    }

    public override string ToString() => $"{V}";

    public static Class Min(Class a, Class b) => a.V.CompareTo(b.V) <= 0 ? a : b;
    public static Class Max(Class a, Class b) => a.V.CompareTo(b.V) >= 0 ? a : b;
    public static (Class min, Class max) MinMax(Class a, Class b) => (Min(a, b), Max(a, b));

}
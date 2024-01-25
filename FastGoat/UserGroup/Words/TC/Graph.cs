using FastGoat.Commons;
using FastGoat.Structures;

namespace FastGoat.UserGroup.Words.TC;

public partial class Graph
{
    public Gen[] Subgroup { get; set; }
    public List<Relator> Relators { get; set; }
    public Gen[] Gens { get; }
    public int NbGens { get; }
    private int Step { get; set; }
    private bool End { get; set; }
    private List<Class> Classes { get; }
    private Graph(string sg, string rels)
    {
        Subgroup = CreateHeader(sg.Split(',', StringSplitOptions.TrimEntries)).SelectMany(c => c).ToArray();
        Relators = CreateHeader(rels.Split(',', StringSplitOptions.TrimEntries))
            .Select((gen, k) => new Relator(k, gen)).ToList();
        Gens = Relators.SelectMany(r => r.Gens)
            .Concat(Subgroup)
            .Select(c => char.ToLower(c.V))
            .Distinct()
            .Order()
            .Where(c => c != Gen.Id)
            .SelectMany(c => new[] { new Gen(c), new Gen(c).Invert() })
            .ToArray();
        NbGens = Gens.Length;

        var Null = new Class(0, this);
        var One = new Class(1, this);
        foreach (var gen in Subgroup)
            One[gen] = One;

        Classes = new() { Null, One };
        End = false;
    }

    private (Class? cl1, Class? cl2) UpdateGraph()
    {
        foreach (var @class in Classes.Where(e => !e.IsComplete))
        {
            var e = @class.UpdateClass();
            if (e.cl1 is not null)
                return e;
        }

        return (null, null);
    }

    private (Class? cl, Gen g) FindCandidate()
    {
        return Classes.OrderBy(c => c.V)
            .Where(c => !c.IsComplete)
            .Select(c => c.FindCandidate())
            .FirstOrDefault(e => e.Item1 is not null, (null, new()));
    }

    private string Coincidence(Class cl1, Class cl2)
    {
        (cl1, cl2) = Class.MinMax(cl1, cl2);
        Classes.RemoveAll(c => c.V == cl2.V);
        var (v1, v2) = (cl1.V, cl2.V);
        foreach (var c in Classes)
            c.Substitute(cl1, cl2);

        foreach (var (c, k) in Classes.Select((c, k) => (c, k)))
        {
            c.IsComplete = false;
            c.V = k;
        }

        var op = $"Replace:({v2}) ~-> ({v1})";
        // Console.WriteLine(op);
        return op;
    }

    private void Display(string action)
    {
        Console.WriteLine($"### Step:{Step} {action}");
        var digits = $"{Classes.Count - 1}".Length;
        var fmt = $"{{0,{digits + 1}}}";
        Console.WriteLine(Header(fmt));
        foreach (var cl in Classes.Skip(1))
            Console.WriteLine(cl.Format(fmt));

        Console.WriteLine();
        var space = string.Format(fmt, ' ');
        Console.WriteLine(space + Gens.Glue("", fmt));
        foreach (var cl in Classes.Skip(1))
            Console.WriteLine(string.Format(fmt, cl.V) + Gens.Select(g => cl.HasEdge(g) ? string.Format(fmt, cl[g]) : space).Glue());
    }

    private void DisplayTable(string fmt)
    {
        var space = string.Format(fmt, ' ');
        Console.WriteLine(space + Gens.Glue("", fmt));
        foreach (var cl in Classes.Skip(1))
            Console.WriteLine(string.Format(fmt, cl.V) + Gens.Select(g => cl.HasEdge(g) ? string.Format(fmt, cl[g]) : space).Glue());
        
    }

    private void Build(bool details = true, bool time = true)
    {
        if (time)
            GlobalStopWatch.AddLap();
        
        if (details)
            Display("Start");

        while (!End)
        {
            if (time && Step >= 50 && Step % 50 == 0)
                Console.WriteLine($"Step:{Step} NbClasses:{Classes.Count - 1}");
            
            var (cl1, cl2) = UpdateGraph();
            if (cl1 is not null)
            {
                var op = Coincidence(cl1, cl2!);
                ++Step;
                if (details)
                    Display(op);
            }
            else
            {
                var (cl3, g) = FindCandidate();
                if (cl3 is not null)
                {
                    var v = Classes.Max(c => c.V) + 1;
                    var cl4 = new Class(v, this);
                    Classes.Add(cl4);
                    cl3[g] = cl4;
                    cl4[g.Invert()] = cl3;
                    var op = $"Add Op:({cl3}) * {g} = ({cl4})";
                    ++Step;
                    if (details)
                        Display(op);
                }
                else
                    End = true;
            }
        }

        if (details)
        {
            Display("End");
            Console.WriteLine();
        }

        if (time)
        {
            Console.WriteLine($"Step:{Step} NbClasses:{Classes.Count - 1}");
            GlobalStopWatch.Show();
            Console.WriteLine();
        }
    }

    private static Gen[][] CreateHeader(params string[] gens)
    {
        // The character 'i' is reserved for neutral subgroup
        var gi = gens.ToList().FindAll(g => g.Contains(Gen.Id));
        if (gi.Count > 1 || (gi.Count == 1 && gi[0].Length != 1))
            throw new GroupException(GroupExceptionType.GroupDef);

        return gens.Select(StringExt.ExpandRelator).OrderBy(w => w.Length).ThenBy(w => w)
            .Select(w => w.Select(c => new Gen(c)).ToArray()).ToArray();
    }

    private string Header(string fmt)
    {
        var s = string.Format(fmt, ' ');
        var s1 = s.Substring(s.Length / 2);
        return s1 + Relators.Select(rel => rel.Format(fmt)).Glue();
    }

    public static void RunTC(string sg, string rel, bool details = true, bool time = true) => new Graph(sg, rel).Build(details, time);
    public static void RunTC(string rel, bool details = true, bool time = true) => RunTC("i", rel, details, time);

    public static void RunTCtable(string rels)
    {
        var graph = new Graph("i", rels);
        graph.Build(false, false);
        var digits = $"{graph.Classes.Count - 1}".Length;
        var fmt = $"{{0,{digits + 1}}}";
        graph.DisplayTable(fmt);
    }
}
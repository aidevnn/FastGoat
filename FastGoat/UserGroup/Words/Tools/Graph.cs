using FastGoat.Commons;
using FastGoat.Structures;

namespace FastGoat.UserGroup.Words.Tools;

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

    private void Build()
    {
        if (Logger.Level != LogLevel.Off)
            GlobalStopWatch.AddLap();
        
        if (Logger.Level == LogLevel.Level2)
            DisplayFancy("Start");

        while (!End)
        {
            if (Logger.Level != LogLevel.Off)
                if (Step >= 50 && Step % 50 == 0)
                    Console.WriteLine($"Step:{Step} NbClasses:{Classes.Count - 1}");
            
            var (cl1, cl2) = UpdateGraph();
            if (cl1 is not null)
            {
                var op = Coincidence(cl1, cl2!);
                ++Step;
                if (Logger.Level == LogLevel.Level2)
                    DisplayFancy(op);
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
                    if (Logger.Level == LogLevel.Level2)
                        DisplayFancy(op);
                }
                else
                    End = true;
            }
        }

        if (Logger.Level == LogLevel.Level2)
        {
            DisplayFancy("End");
            Console.WriteLine();
        }

        if (Logger.Level != LogLevel.Off)
        {
            Console.WriteLine($"Step:{Step} NbClasses:{Classes.Count - 1}");
            GlobalStopWatch.Show();
            Console.WriteLine();
        }
    }

    #region Display Fancy
    public int[] Separators
    {
        get
        {
            List<int> seps = new() { 0 };
            var arrGens = Relators.ToArray();
            foreach (var g in arrGens)
                seps.Add(seps.Last() + g.Length);

            return seps.ToArray();
        }
    }
    void DisplayFancy(string action)
    {
        Console.WriteLine($"### Step:{Step} {action}");
        DisplayTableRelators();
        DisplayTableOps();
    }

    void DisplayTableRelators()
    {
        var digits = $"{Classes.Count - 1}".Length;
        var fmt = $"{{0,{digits + 1}}}";
        Console.WriteLine("# Relators table");
        Console.WriteLine($" {Relators.Select(r => r.Format(fmt)).Glue()}");
        DisplayLineDown(digits);
        foreach (var kv in Classes.Skip(1))
        {
            if (kv.V > 1 && kv.V % 40 == 1)
            {
                DisplayLineUp(digits);
                Console.WriteLine($" {Relators.Select(r => r.Format(fmt)).Glue()}");
                DisplayLineDown(digits);
            }
            Console.WriteLine(kv.Display(digits));
        }

        DisplayLineUp(digits);
        Console.WriteLine();
    }
    
    void DisplayLineUp(int digits)
    {
        var s1 = Enumerable.Repeat('─', (digits + 1) * (Relators.Sum(r => r.Length) + 1)).Append(' ').Glue().ToArray();
        foreach (var k in Separators)
            s1[(digits + 1) * k] = s1[(digits + 1) * (k + 1)] = '┴';

        s1[0] = '└';
        s1[^1] = '┘';
        Console.WriteLine(s1.Glue());
    }
    
    void DisplayLineDown(int digits)
    {
        var s1 = Enumerable.Repeat('─', (digits + 1) * (Relators.Sum(r => r.Length) + 1)).Append(' ').Glue().ToArray();
        foreach (var k in Separators)
            s1[(digits + 1) * k] = s1[(digits + 1) * (k + 1)] = '┬';

        s1[0] = '┌';
        s1[^1] = '┐';
        Console.WriteLine(s1.Glue());
    }
    
    public void DisplayTableOps()
    {
        Console.WriteLine("# Classes table");
        var gens = Gens.OrderByDescending(c => char.IsLower(c.V)).ThenBy(c => c).ToArray();
        var digits = $"{Classes.Count - 1}".Length;
        var fmt = $"{{0,{digits + 1}}}";
        var head = string.Format(fmt, "") + "│" + gens.Glue("│", fmt) + "│";
        var lineTop = Enumerable.Range(0, head.Length - (digits + 2)).Select(i =>
            i % (digits + 2) == digits + 1 ? (i == head.Length - digits - 3 ? '┐' : '┬') : '─').Glue();
        var lineMid = Enumerable.Range(0, head.Length)
            .Select(i => i % (digits + 2) == digits + 1 ? (i == head.Length - 1 ? '┤' : '┼') : '─').Glue();
        var lineEnd = Enumerable.Range(0, head.Length)
            .Select(i => i % (digits + 2) == digits + 1 ? (i == head.Length - 1 ? '┘' : '┴') : '─').Glue();

        var rows = new List<string>();
        foreach (var i in Classes.Skip(1))
        {
            var r = gens.Select(g => i.HasEdge(g) ? $"{i[g]}" : "").Prepend($"{i}").Glue("│", fmt);
            rows.Add(r.Glue());
        }

        Console.WriteLine(" " + string.Format(fmt, "") + "┌" + lineTop);
        Console.WriteLine(" " + head);
        Console.WriteLine("┌" + lineMid);
        rows.ForEach(r => Console.WriteLine("│" + r + "│"));
        Console.WriteLine("└" + lineEnd);
        Console.WriteLine();
    }
    #endregion

    private static Gen[][] CreateHeader(params string[] gens)
    {
        // The character 'i' is reserved for neutral subgroup
        var gi = gens.ToList().FindAll(g => g.Contains(Gen.Id));
        if (gi.Count > 1 || (gi.Count == 1 && gi[0].Length != 1))
            throw new GroupException(GroupExceptionType.GroupDef);

        return gens.Select(StringExt.ExpandRelator).OrderBy(w => w.Length).ThenBy(w => w)
            .Select(w => w.Select(c => new Gen(c)).ToArray()).ToArray();
    }

    public static void RunToddCoxeterAlgo(string sg, string rel) => new Graph(sg, rel).Build();
    public static void RunToddCoxeterAlgo(string rel) => RunToddCoxeterAlgo("i", rel);
}
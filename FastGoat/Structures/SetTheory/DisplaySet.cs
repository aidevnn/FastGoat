namespace FastGoat.Structures.SetTheory;

public class DisplaySet<U> where U : struct, IElt<U>
{
    protected static List<string> GenLetters(int n)
    {
        if (n > 50)
            return Enumerable.Range(1, n).Select(a => $"E{a,2:000}").ToList();

        return "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ".Take(n).Select(c => $"{c}").ToList();
    }

    static string MyFormat(char c, string g, List<char> l) => string.Format("{0, 2}|{1}", c, string.Join(g, l));
    static string MyFormat(string c, string g, List<char> l) => string.Format("{0, 2}|{1}", c.Substring(0, 2), string.Join(g, l));

    public DisplaySet(SubSet<U> sub)
    {
        G = sub;
        SpecialElts = new Dictionary<U, char>(new EltEquality<U>());
    }

    SubSet<U> G { get; }
    Dictionary<U, char> SpecialElts { get; }
    public void SpecialChar(params (char, U)[] spec)
    {
        foreach (var (c, e) in spec) SpecialElts[e] = c;
    }

    public string Name { get; set; } = "G";
    string infos = "";
    string fmt = "|{0}| = {1} {2}";
    public void SetDetails(string? name = null, string? fmt = null, string? infos = null)
    {
        Name = name ?? Name;
        this.fmt = fmt ?? this.fmt;
        this.infos = infos ?? this.infos;
    }

    public virtual void DisplayHead()
    {
        Console.WriteLine(fmt, Name, G.Count, infos);
    }

    protected virtual string DisplayElement(U e, string name) => string.Format("{0} = {1}", name, e);

    public virtual void DisplayElements()
    {
        var elts = G.AllElements.ToList();
        if (elts.Count == 0)
        {
            Console.WriteLine("Empty Set");
            return;
        }

        DisplayHead();
        if (elts.Count > 300)
        {
            Console.WriteLine("TOO BIG");
            return;
        }

        elts.Sort(G.CompareElt);
        var word = new Queue<string>(GenLetters(elts.Count));
        foreach (var e in elts)
        {
            var w = SpecialElts.ContainsKey(e) ? SpecialElts[e].ToString() : word.Dequeue();
            Console.WriteLine(DisplayElement(e, w));
        }

        Console.WriteLine();
    }

    public virtual void DisplayTable(char symb, Func<U, U, U> Op)
    {
        var elts = G.AllElements.ToList();
        if (elts.Count == 0)
        {
            Console.WriteLine("Empty Set");
            return;
        }

        DisplayHead();

        if (elts.Count > 50)
        {
            Console.WriteLine("TOO BIG");
            return;
        }

        elts.Sort(G.CompareElt);
        var word = new Queue<char>(string.Join("", GenLetters(elts.Count)));
        Dictionary<char, U> ce = new Dictionary<char, U>();
        Dictionary<U, char> ec = new Dictionary<U, char>(new EltEquality<U>());

        foreach (var e in elts)
        {
            var c = SpecialElts.ContainsKey(e) ? SpecialElts[e] : word.Dequeue();
            ce[c] = e;
            ec[e] = c;
        }

        var rows = new List<(char, List<char> row)>();
        foreach (var e0 in elts)
        {
            var l = new List<char>();
            rows.Add((ec[e0], l));
            foreach (var e1 in elts)
            {
                var e2 = Op(e0, e1);
                if (!G.Contains(e2))
                    return;

                l.Add(ec[e2]);
            }
        }

        var head = MyFormat(symb, " ", rows.Select(a => a.Item1).ToList());
        var line = MyFormat("--", "", Enumerable.Repeat('-', rows.Count * 2).ToList());
        Console.WriteLine(head);
        Console.WriteLine(line);

        foreach (var (e0, row) in rows)
            Console.WriteLine(MyFormat(e0, " ", row));

        Console.WriteLine();
    }
}
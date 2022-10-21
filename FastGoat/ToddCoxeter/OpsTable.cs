namespace FastGoat.ToddCoxeter;

public class OpsTable
{
    Header sgHeader { get; }
    Header relHeader { get; }
    SortedDictionary<OpKey, Symbol> opsTable { get; }

    public OpsTable(Header sgheader, Header relheader)
    {
        if (!sgheader.Generators.IsSubsetOf(relheader.Generators))
            throw new Exception();

        sgHeader = new(sgheader);
        relHeader = new(relheader);
        sgTable = new(sgHeader);
        rTable = new(relHeader);
        opsTable = new();
    }

    public OpsTable(OpsTable table)
    {
        Previous = table ?? throw new Exception();
        sgHeader = table.sgHeader;
        relHeader = table.relHeader;
        sgTable = new(table.sgTable);
        rTable = new(table.rTable);
        opsTable = new(table.opsTable);
    }

    public SubGroupTable sgTable { get; }
    public RelatorsTable rTable { get; }
    public OpsTable? Previous { get; }

    public void BuildTable()
    {
        int sz = 0;
        HashSet<Op> newOps = new();
        (Symbol, Symbol) err = new();
        do
        {
            newOps.Clear();
            sz = sgTable.CountUnknown + rTable.CountUnknown;
            err = sgTable.ApplyOp(opsTable, newOps);
            if (Substitute(err.Item1, err.Item2))
                foreach (var op in newOps)
                    ApplyOp(op);

            newOps.Clear();
            err = rTable.ApplyOp(opsTable, newOps);
            if (Substitute(err.Item1, err.Item2))
                foreach (var op in newOps)
                    ApplyOp(op);
        } while (newOps.Count != 0 || sz != sgTable.CountUnknown + rTable.CountUnknown);
    }

    bool Substitute(Symbol s0, Symbol s1)
    {
        if (s0 == Symbol.Unknown)
            return true;

        opsTable.Clear();

        sgTable.Subtitute(s0, s1);
        rTable.SubtituteRemove(s0, s1);

        while (rTable.ContainsKey(s1.Next))
        {
            s0 = s1;
            s1 = s0.Next;
            sgTable.Subtitute(s0, s1);
            rTable.SubtituteWithKey(s0, s1);
        }

        return false;
    }

    public void ApplyOp(Op op)
    {
        var opKey = new OpKey(op.i, op.g);
        var opiKey = new OpKey(op.j, op.g.Invert());
        if (!opsTable.ContainsKey(opKey) && !opsTable.ContainsKey(opiKey))
        {
            opsTable[opKey] = op.j;
            opsTable[opiKey] = op.i;
        }
        else
        {
            if (opsTable.ContainsKey(opKey))
            {
                var (s0, s1) = Symbol.MinMax(op.j, opsTable[opKey]);
                Substitute(s0, s1);
            }

            if (opsTable.ContainsKey(opiKey))
            {
                var (s0, s1) = Symbol.MinMax(op.i, opsTable[opiKey]);
                Substitute(s0, s1);
            }
        }
    }

    public Op FirstOp()
    {
        if (opsTable.Count == 0)
            return new();

        var fop = opsTable.First();
        return new(fop.Key.i, fop.Key.g, fop.Value);
    }

    public Op NewOp()
    {
        var j = opsTable.SelectMany(e => new[] { e.Key.i, e.Value }).Descending().FirstOrDefault().Next;

        var sgOp = sgTable.GetOps().FirstOrDefault(op =>
            op.i != Symbol.Unknown && op.g != Generator.Unknown && op.j == Symbol.Unknown);
        if (sgOp.i != Symbol.Unknown && sgOp.g != Generator.Unknown && sgOp.j == Symbol.Unknown)
            return new Op(sgOp.i, sgOp.g, j);

        var rOp = rTable.GetOps().FirstOrDefault(op =>
            op.i != Symbol.Unknown && op.g != Generator.Unknown && op.j == Symbol.Unknown);
        if (rOp.i != Symbol.Unknown && rOp.g != Generator.Unknown && rOp.j == Symbol.Unknown)
            return new Op(rOp.i, rOp.g, j);

        return new();
    }

    string TableFind(OpKey opk)
    {
        if (opsTable.ContainsKey(opk))
            return opsTable[opk].ToString();

        return " ";
    }

    private Dictionary<Symbol, IEnumerable<char>> elementsTable;
    private char[] generators;

    void ChainSymbol(Symbol s, List<OpKey> chain)
    {
        if (s.Equals(Symbol.One))
            return;

        var prev = opsTable.First(kp => kp.Value.Equals(s));
        chain.Add(prev.Key);
        var s0 = prev.Key.i;
        ChainSymbol(s0, chain);
    }

    public char[] Generators() => generators;
    public IEnumerable<IEnumerable<char>> Words() => elementsTable.Values;
    public void GenerateWords()
    {
        generators = opsTable.Keys.Select(k => k.g.GetLowerCase()).Distinct().Where(k => k != Generator.Id.Value)
            .Ascending().ToArray();
        var symbs = opsTable.Keys.Select(k => k.i).Distinct().Ascending().ToArray();
        elementsTable = new();
        foreach (var symbol in symbs)
        {
            List<OpKey> chain = new();
            ChainSymbol(symbol,chain);
            elementsTable[symbol] = chain.Select(o => o.g.Value);
        }
    }

    public IEnumerable<char> Rewrite(IEnumerable<char> w)
    {
        if (!w.Select(char.ToLower).All(generators.Contains))
            throw new GroupException(GroupExceptionType.GroupDef);
        
        var symbol = w.Reverse().Select(c => new Generator(c)).Aggregate(Symbol.One, (i, g) => opsTable[new OpKey(i, g)]);
        return elementsTable[symbol];
    }

    public void DisplayTable(int digits)
    {
        var fmt = $"{{0,{digits + 1}}}";
        Console.WriteLine("# Classes table");
        var gens = opsTable.Keys.Select(k => k.g).Distinct().Ascending();
        var symbs = opsTable.Keys.Select(k => k.i).Distinct().Ascending();

        var head = string.Format(fmt, "") + "│" + gens.Glue("│", fmt) + "│";
        var lineTop = Enumerable.Range(0, head.Length - (digits + 2)).Select(i =>
            i % (digits + 2) == digits + 1 ? (i == head.Length - digits - 3 ? '┐' : '┬') : '─').Glue();
        var lineMid = Enumerable.Range(0, head.Length)
            .Select(i => i % (digits + 2) == digits + 1 ? (i == head.Length - 1 ? '┤' : '┼') : '─').Glue();
        var lineEnd = Enumerable.Range(0, head.Length)
            .Select(i => i % (digits + 2) == digits + 1 ? (i == head.Length - 1 ? '┘' : '┴') : '─').Glue();

        var rows = new List<string>();
        foreach (var i in symbs)
        {
            var r = gens.Select(g => TableFind(new(i, g))).Prepend(i.ToString()).Glue("│", fmt);
            rows.Add(r.Glue());
        }

        Console.WriteLine(" " + string.Format(fmt, "") + "┌" + lineTop);
        Console.WriteLine(" " + head);
        Console.WriteLine("┌" + lineMid);
        rows.ForEach(r => Console.WriteLine("│" + r + "│"));
        Console.WriteLine("└" + lineEnd);
        Console.WriteLine();
    }

    public void Display()
    {
        var digits = opsTable.Count == 0 ? 2 : opsTable.Max(p0 => Symbol.Max(p0.Key.i, p0.Value).ToString().Length);

        sgTable.Display(digits);
        rTable.Display(digits);
        DisplayTable(digits);
        Console.WriteLine();
    }

    public void DisplayOps()
    {
        var digits = opsTable.Count == 0 ? 2 : opsTable.Max(p0 => Symbol.Max(p0.Key.i, p0.Value).ToString().Length);
        DisplayTable(digits);
        Console.WriteLine();
    }

    public static Header CreateHeader(params string[] gens)
    {
        // The character 'i' is reserved for neutral subgroup
        var gi = gens.ToList().FindAll(g => g.Contains(Generator.Id.Value));
        if (gi.Count > 1 || (gi.Count == 1 && gi[0].Length != 1))
            throw new GroupException(GroupExceptionType.GroupDef);
        
        var head = gens.Select(Group.ExpandRelator).OrderBy(w => w.Length).ThenBy(w => w)
            .Select(w => w.Select(c => new Generator(c))).ToArray();

        return new Header(head);
    }
}
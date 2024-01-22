using FastGoat.Commons;
using FastGoat.Structures;

namespace FastGoat.UserGroup.Words.ToddCoxeter;

public class OpsTable
{
    Header sgHeader { get; }
    Header relHeader { get; }
    SortedDictionary<OpKey, EqClass> opsTable { get; }

    public OpsTable(Header sgheader, Header relheader)
    {
        if (!sgheader.Generators.IsSubsetOf(relheader.Generators))
            throw new Exception();

        sgHeader = new(sgheader);
        relHeader = new(relheader);
        sgTable = new(sgHeader);
        rTable = new(relHeader);
        opsTable = new();

        generators = Array.Empty<char>();
        elementsTable = new();
    }

    public OpsTable(OpsTable table)
    {
        Previous = table ?? throw new Exception();
        sgHeader = table.sgHeader;
        relHeader = table.relHeader;
        sgTable = new(table.sgTable);
        rTable = new(table.rTable);
        opsTable = new(table.opsTable);

        generators = Array.Empty<char>();
        elementsTable = new();
    }

    public SubGroupTable sgTable { get; }
    public RelatorsTable rTable { get; }
    public OpsTable? Previous { get; }

    public void BuildTable(bool details = false)
    {
        int sz = 0;
        HashSet<Op> newOps = new();
        do
        {
            newOps.Clear();
            sz = sgTable.CountUnknown + rTable.CountUnknown;
            var errSg = sgTable.ApplyOp(opsTable, newOps);
            if (NoCoincidence(errSg.Item1, errSg.Item2, details))
                foreach (var op in newOps)
                    ApplyOp(op);

            newOps.Clear();
            var errRel = rTable.ApplyOp(opsTable, newOps);
            if (NoCoincidence(errRel.Item1, errRel.Item2, details))
                foreach (var op in newOps)
                    ApplyOp(op);
        } while (newOps.Count != 0 || sz != sgTable.CountUnknown + rTable.CountUnknown);
    }
    
    bool NoCoincidence(EqClass s0, EqClass s1, bool details = false)
    {
        if (s0 == EqClass.Unknown)
            return true;
        
        var digits = opsTable.Values.Concat(opsTable.Keys.Select(e => e.i)).Distinct().Max(e => $"{e}".Length);
        if (details)
        {
            Console.WriteLine("# Coincidence detected");
            rTable.Display(digits);
            Console.WriteLine($"({s1}) ~-> ({s0})");
        }
        
        sgTable.Subtitute(s0, s1);
        rTable.SubtituteRemove(s0, s1);
        Substitute(s0, s1);
        
        if (details)
            rTable.Display(digits);

        while (rTable.ContainsKey(s1.Next))
        {
            (s0, s1) = (s1, s1.Next);
            if (details)
            {
                Console.WriteLine($"({s1}) ~-> ({s0})");
                rTable.Display(digits);
            }
            
            sgTable.Subtitute(s0, s1);
            rTable.SubtituteAddRemove(s0, s1);
            Substitute(s0, s1);
        }
        
        return false;
    }
    
    void Substitute(EqClass s0, EqClass s1)
    {
        foreach (var op in opsTable.Where(e => e.Value.Equals(s1) || e.Key.i.Equals(s1)).ToArray())
        {
            opsTable.Remove(op.Key);
            if (op.Key.i.Equals(s1) && op.Value.Equals(s1))
                opsTable[new OpKey(s0, op.Key.g)] = s0;
            else if (op.Key.i.Equals(s1) && !op.Value.Equals(s1))
                opsTable[new OpKey(s0, op.Key.g)] = op.Value;
            else if (!op.Key.i.Equals(s1) && op.Value.Equals(s1))
                opsTable[new OpKey(op.Key.i, op.Key.g)] = s0;
        }
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
            if (opsTable.TryGetValue(opKey, out var s))
            {
                var (s0, s1) = EqClass.MinMax(op.j, s);
                NoCoincidence(s0, s1);
            }

            if (opsTable.TryGetValue(opiKey, out var si))
            {
                var (s0, s1) = EqClass.MinMax(op.i, si);
                NoCoincidence(s0, s1);
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
            op.i != EqClass.Unknown && op.g != Generator.Unknown && op.j == EqClass.Unknown);
        if (sgOp.i != EqClass.Unknown && sgOp.g != Generator.Unknown && sgOp.j == EqClass.Unknown)
            return new Op(sgOp.i, sgOp.g, j);

        var rOp = rTable.GetOps().FirstOrDefault(op =>
            op.i != EqClass.Unknown && op.g != Generator.Unknown && op.j == EqClass.Unknown);
        if (rOp.i != EqClass.Unknown && rOp.g != Generator.Unknown && rOp.j == EqClass.Unknown)
            return new Op(rOp.i, rOp.g, j);

        return new();
    }

    string TableFind(OpKey opk)
    {
        if (opsTable.TryGetValue(opk, out var sk))
            return sk.ToString();

        return " ";
    }

    private Dictionary<EqClass, IEnumerable<char>> elementsTable;
    private char[] generators;

    void ChainSymbol(EqClass s, List<OpKey> chain)
    {
        if (s.Equals(EqClass.One))
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
            ChainSymbol(symbol, chain);
            elementsTable[symbol] = chain.Select(o => o.g.Value);
        }
    }

    public IEnumerable<char> Rewrite(IEnumerable<char> w)
    {
        if (!w.Select(char.ToLower).All(generators.Contains))
            throw new GroupException(GroupExceptionType.GroupDef);

        var symbol = w.Reverse().Select(c => new Generator(c))
            .Aggregate(EqClass.One, (i, g) => opsTable[new OpKey(i, g)]);
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
        var digits = opsTable.Count == 0 ? 2 : opsTable.Max(p0 => EqClass.Max(p0.Key.i, p0.Value).ToString().Length);

        sgTable.Display(digits);
        rTable.Display(digits);
        DisplayTable(digits);
        Console.WriteLine();
    }

    public void DisplayOps()
    {
        var digits = opsTable.Count == 0 ? 2 : opsTable.Max(p0 => EqClass.Max(p0.Key.i, p0.Value).ToString().Length);
        DisplayTable(digits);
        Console.WriteLine();
    }

    public static Header CreateHeader(params string[] gens)
    {
        // The character 'i' is reserved for neutral subgroup
        var gi = gens.ToList().FindAll(g => g.Contains(Generator.Id.Value));
        if (gi.Count > 1 || (gi.Count == 1 && gi[0].Length != 1))
            throw new GroupException(GroupExceptionType.GroupDef);

        var head = gens.Select(StringExt.ExpandRelator).OrderBy(w => w.Length).ThenBy(w => w)
            .Select(w => w.Select(c => new Generator(c))).ToArray();

        return new Header(head);
    }
}
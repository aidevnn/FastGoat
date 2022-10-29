using FastGoat.Commons;

namespace FastGoat.ToddCoxeter;

public class RelatorsTable
{
    Dictionary<Symbol, Line> table { get; }
    Header header { get; }

    public RelatorsTable(Header head)
    {
        table = new();
        header = head;
        table[Symbol.One] = NewLine(Symbol.One);
    }

    public RelatorsTable(RelatorsTable rTable)
    {
        header = rTable.header;
        table = rTable.table.ToDictionary(s => s.Key, l => new Line(l.Value));
    }

    public bool ContainsKey(Symbol s) => table.ContainsKey(s);
    Line NewLine(Symbol s) => new(s, header);
    public IEnumerable<Op> GetOps() => table.SelectMany(kv => kv.Value.GetOps());
    public int CountUnknown => table.Sum(r => r.Value.CountUnknown);
    public void Remove(Symbol s) => table.Remove(s);

    public void SubtituteRemove(Symbol s0, Symbol s1)
    {
        table.Remove(s1);
        foreach (var e in table)
            e.Value.Subtitute(s0, s1);
    }

    public void SubtituteWithKey(Symbol s0, Symbol s1)
    {
        if (!table.ContainsKey(s0))
        {
            var line = table[s1];
            table[s0] = new(s0, line);
            table.Remove(s1);
        }
        else
            throw new Exception("TO DO");

        foreach (var e in table)
            e.Value.Subtitute(s0, s1);
    }

    public (Symbol, Symbol) ApplyOp(SortedDictionary<OpKey, Symbol> opsTable, HashSet<Op> newOps)
    {
        var symbols = opsTable.Values.Distinct().Ascending();
        foreach (var s in symbols)
        {
            if (!table.ContainsKey(s))
                table[s] = NewLine(s);
        }

        foreach (var kv in table)
        {
            var err = kv.Value.ApplyOp(opsTable, newOps);
            if (err.Item1 != Symbol.Unknown)
                return err;
        }

        return new();
    }

    public void Display(int digits)
    {
        Console.WriteLine("# Relators table");
        header.DisplayHead(digits);
        int k = 0;
        foreach (var kv in table.OrderBy(p => p.Key))
        {
            if (k > 0 && k % 40 == 0) header.ReDisplayHead(digits);
            Console.WriteLine(kv.Value.Display(digits));
            ++k;
        }

        header.DisplayLineUp(digits);
        Console.WriteLine();
    }
}
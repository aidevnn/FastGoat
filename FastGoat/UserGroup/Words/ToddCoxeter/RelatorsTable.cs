using FastGoat.Commons;

namespace FastGoat.UserGroup.Words.ToddCoxeter;

public class RelatorsTable
{
    Dictionary<EqClass, Line> table { get; }
    Header header { get; }

    public RelatorsTable(Header head)
    {
        table = new();
        header = head;
        table[EqClass.One] = NewLine(EqClass.One);
    }

    public RelatorsTable(RelatorsTable rTable)
    {
        header = rTable.header;
        table = rTable.table.ToDictionary(s => s.Key, l => new Line(l.Value));
    }

    public bool ContainsKey(EqClass s) => table.ContainsKey(s);
    Line NewLine(EqClass s) => new(s, header);
    public IEnumerable<Op> GetOps() => table.SelectMany(kv => kv.Value.GetOps());
    public int CountUnknown => table.Sum(r => r.Value.CountUnknown);
    public void Remove(EqClass s) => table.Remove(s);

    public void SubtituteRemove(EqClass s0, EqClass s1)
    {
        table.Remove(s1);
        foreach (var e in table)
            e.Value.Subtitute(s0, s1);
    }

    public void SubtituteAddRemove(EqClass s0, EqClass s1)
    {
        if (!table.ContainsKey(s0))
        {
            var line = table[s1];
            table[s0] = new(s0, line);
            table.Remove(s1);
        }

        foreach (var e in table)
            e.Value.Subtitute(s0, s1);
    }

    public (EqClass, EqClass) ApplyOp(SortedDictionary<OpKey, EqClass> opsTable, HashSet<Op> newOps)
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
            if (err.Item1 != EqClass.Unknown)
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
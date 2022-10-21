namespace FastGoat.ToddCoxeter;

public class SubGroupTable
{
    Header header { get; }
    Line line { get; }
    public SubGroupTable(Header head)
    {
        header = head;
        line = new(Symbol.One, header);
    }
    public SubGroupTable(SubGroupTable sgtable)
    {
        header = sgtable.header;
        line = new(sgtable.line);
    }
    public int CountUnknown => line.CountUnknown;
    public void Subtitute(Symbol s0, Symbol s1) => line.Subtitute(s0, s1);
    public IEnumerable<Op> GetOps() => line.GetOps();
    public (Symbol, Symbol) ApplyOp(SortedDictionary<OpKey, Symbol> opsTable, HashSet<Op> newOps)
    {
        return line.ApplyOp(opsTable, newOps);
    }
    public void Display(int digits)
    {
        Console.WriteLine("# SubGroup table");
        header.DisplayHead(digits);
        Console.WriteLine(line.Display(digits));
        header.DisplayLineUp(digits);
        Console.WriteLine();
    }
}

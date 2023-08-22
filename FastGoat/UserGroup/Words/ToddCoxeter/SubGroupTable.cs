namespace FastGoat.UserGroup.Words.ToddCoxeter;

public class SubGroupTable
{
    Header header { get; }
    Line line { get; }

    public SubGroupTable(Header head)
    {
        header = head;
        line = new(EqClass.One, header);
    }

    public SubGroupTable(SubGroupTable sgtable)
    {
        header = sgtable.header;
        line = new(sgtable.line);
    }

    public int CountUnknown => line.CountUnknown;
    public void Subtitute(EqClass s0, EqClass s1) => line.Subtitute(s0, s1);
    public IEnumerable<Op> GetOps() => line.GetOps();

    public (EqClass, EqClass) ApplyOp(SortedDictionary<OpKey, EqClass> opsTable, HashSet<Op> newOps)
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
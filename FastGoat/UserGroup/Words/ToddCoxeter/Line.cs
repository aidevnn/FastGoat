using FastGoat.Commons;

namespace FastGoat.UserGroup.Words.ToddCoxeter;

public class Line
{
    EqClass[] row { get; }
    Header header { get; }

    public Line(EqClass key, Header head)
    {
        row = new EqClass[head.Count + 1];
        header = head;
        Key = key;

        foreach (var k in header.Separators)
            row[k] = key;
    }

    public Line(Line line)
    {
        header = line.header;
        row = line.row.Select(s => s).ToArray();
        Key = line.Key;
    }

    public Line(EqClass key, Line line)
    {
        header = line.header;
        row = line.row.Select(s => s).ToArray();
        Key = key;
    }

    public int CountUnknown => row.Count(s => s == EqClass.Unknown);

    public void Subtitute(EqClass s0, EqClass s1)
    {
        for (int k = 0; k < row.Length; ++k)
        {
            if (row[k] == s1)
                row[k] = s0;
        }
    }

    public (EqClass, EqClass) ApplyOp(SortedDictionary<OpKey, EqClass> opsTable, HashSet<Op> newOps)
    {
        for (int k = 0; k < header.Count; ++k)
        {
            var i = row[k];
            var g = header[k];
            var j = row[k + 1];
            if (i == EqClass.Unknown && j == EqClass.Unknown)
                continue;
            else if (i != EqClass.Unknown && j != EqClass.Unknown)
            {
                var opKey = new OpKey(i, g);
                var opiKey = new OpKey(j, g.Invert());
                var opCheck = opsTable.ContainsKey(opKey);
                var opiCheck = opsTable.ContainsKey(opiKey);

                if (opCheck && opsTable[opKey] != j)
                    return EqClass.MinMax(j, opsTable[opKey]);

                if (opiCheck && opsTable[opiKey] != i)
                    return EqClass.MinMax(i, opsTable[opiKey]);

                if (!opCheck && !opiCheck)
                    newOps.Add(Op.Create(i, g, j));
            }
            else if (j == EqClass.Unknown)
            {
                var opKey = new OpKey(i, g);
                if (opsTable.TryGetValue(opKey, out var s))
                    row[k + 1] = s;
            }
            else // if (i == Symbol.Unknown)
            {
                var opiKey = new OpKey(j, g.Invert());
                if (opsTable.TryGetValue(opiKey, out var s))
                    row[k] = s;
            }
        }

        return new();
    }

    public IEnumerable<Op> GetOps() => header.Select((g, k) => new Op(row[k], g, row[k + 1]));
    public EqClass Key { get; }

    public string Display(int digits)
    {
        var fmt = $"{{0,{digits + 1}}}";
        var s0 = row.Glue("", fmt);
        var s1 = (s0 + " ").ToArray();
        foreach (var k in header.Separators)
            s1[(digits + 1) * k] = s1[(digits + 1) * (k + 1)] = '│';

        var s2 = s1.Glue();
        return s2;
    }

    public override string ToString()
    {
        var digits = row.Max(s => s.ToString().Length);
        return Display(digits);
    }
}
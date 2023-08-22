using FastGoat.Commons;

namespace FastGoat.UserGroup.Words.ToddCoxeter;

public static class ToddCoxeterAlgo
{
    public static OpsTable Run(string rels, bool details = false)
    {
        var id = Generator.Id.Value;
        return Run($"{id}", $"{id}, {rels}", details);
    }

    public static OpsTable Run(string sg, string rels, bool details = false)
    {
        if (details)
            GlobalStopWatch.Restart();

        var gHeader = OpsTable.CreateHeader(sg.Split(',', StringSplitOptions.TrimEntries));
        var rHeader = OpsTable.CreateHeader(rels.Split(',', StringSplitOptions.TrimEntries));
        var tOps = new OpsTable(gHeader, rHeader);
        tOps.BuildTable();
        int k = 1;
        var fop = tOps.FirstOp();
        if (details)
        {
            Console.WriteLine($"#### Step {k++} Op : {fop} ####");
            tOps.Display();
        }

        while (true)
        {
            var op = tOps.NewOp();
            if (op.i == EqClass.Unknown)
                break;

            if (details)
                Console.WriteLine($"#### Step {k++} Op : {op} ####");

            tOps = new(tOps);
            tOps.ApplyOp(op);
            tOps.BuildTable();

            if (details)
                tOps.Display();
        }

        if (details)
        {
            Console.WriteLine($"####     End    ####");
            GlobalStopWatch.Show("TC");
        }

        tOps.GenerateWords();
        return tOps;
    }
}
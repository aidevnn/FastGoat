using System.Numerics;
using FastGoat.Commons;

namespace FastGoat.UserGroup.EllCurve;

public record EllDB(string name, BigInteger conductor, int rank, int[] torsType, BigInteger[] model)
{
    public int ordTors => torsType.Aggregate((ai, aj) => ai * aj);

    public override string ToString() => $"[{model.Glue(", ")}]; {name}; {conductor}; ({rank}) x [{torsType.Glue(", ")}]";
}
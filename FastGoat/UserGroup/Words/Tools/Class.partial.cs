namespace FastGoat.UserGroup.Words.Tools;

public partial class Class
{
    public static Class? Null => null;
    public Gen STGen { get; set; }
    public Class? STClass { get; set; }
    public List<Gen> Word { get; } = new();
    public List<Gen> WordInv  { get; } = new();
    public Dictionary<Gen, bool> Coloured { get; set; } = new();

    public void Color(Gen g)
    {
        Coloured[g] = Edges[g].Coloured[g.Invert()] = true;
    }
}
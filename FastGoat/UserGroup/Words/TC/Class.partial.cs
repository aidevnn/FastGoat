namespace FastGoat.UserGroup.Words.TC;

public partial class Class
{
    public static Class? Null => null;
    public bool Coloured { get; set; }
    public Gen STGen { get; set; }
    public Class? STClass { get; set; }
    public List<Gen> Word { get; } = new();
}
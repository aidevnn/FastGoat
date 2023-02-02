using FastGoat.UserGroup.Polynoms;

namespace FastGoat.UserGroup.Characters;

public class CnfCell : ACell
{
    public bool IsEmpty => !Content.HasValue;
    private Cnf? Content { get; }
    public Cnf E => Content ?? Cnf.CnfZero;

    public CnfCell()
    {
    }

    public CnfCell(Cnf c)
    {
        Content = c;
    }

    public override string Display() => IsEmpty ? "#" : $"{Content}";
}
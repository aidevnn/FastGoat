using FastGoat.Structures;

namespace FastGoat.UserGroup.Characters;

public abstract class ACell : IElt<ACell>
{
    public abstract string Display();
    public bool Equals(ACell? other) => false;
    public int CompareTo(ACell? other) => 0;
    public int Hash => GetHashCode();
    public override string ToString() => Display();
}
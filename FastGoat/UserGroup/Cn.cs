namespace FastGoat.UserGroup;

public class Cn : ConcreteGroup<ZnInt>
{
    public Cn(int n) : base($"C{n}", new Zn(n), new[] { new Zn(n)[1] })
    {
        Hash = (BaseGroup.Hash, "Cn").GetHashCode();
    }
}
namespace FastGoat.UserGroup;

public class Symm : ConcreteGroup<Perm>
{
    public Symm(int n) : base($"Symm{n}", new Sn(n))
    {
        
    }
}
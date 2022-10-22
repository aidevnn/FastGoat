using FastGoat.UserGroup;

namespace FastGoat;

public class ZentrumChain<T> where T : struct, IElt<T>
{
    public ZentrumChain(ConcreteGroup<T> g, int i = 0)
    {
        G = g;
        Z = Group.Zentrum(g);
        var zi = g.Over(Z, $"Z{i}");
        if (Z.Count() != 1)
        {
            Next = new ZentrumChain<Representative<T>>(zi, i + 1);
        }
    }
        
    public ConcreteGroup<T> G { get; }
    public ConcreteGroup<T> Z { get; }
    public ZentrumChain<Representative<T>>? Next { get; }
}
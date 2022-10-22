using FastGoat.UserGroup;

namespace FastGoat.Examples;

public static class NilpotentGroups
{
    static void ZentrumChain<T>(ConcreteGroup<T> g, int i = 1) where T : struct, IElt<T>
    {
        if(g.Count() == 1)
            return;
        
        var Z = Group.Zentrum(g);
        if(Z.Count() == 1)
            return;
        
        var zi = g.Over(Z, $"Z{i}");
        DisplayGroup.Head(zi);

        ZentrumChain(zi, i + 1);
    }

    public static void Nilpotent2()
    {
        var gr = new WordGroup("Q8", "a4, a2=b2, a-1=bab-1");
        DisplayGroup.Head(gr);
        ZentrumChain(gr);
    }

    public static void Nilpotent3()
    {
        var gr = new WordGroup("Q16", "a8, b2=a4, bab-1=a-1");
        DisplayGroup.Head(gr);
        ZentrumChain(gr);
    }

    public static void Nilpotent4()
    {
        var gr = new WordGroup("Q32", "a16, b2=a8, bab-1=a-1");
        DisplayGroup.Head(gr);
        ZentrumChain(gr);
    }
}
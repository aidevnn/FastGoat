using FastGoat.Commons;

namespace FastGoat.UserGroup.Words.Tools;

public readonly partial struct Circuit
{
    public Class Class { get; }
    public Relator Relator { get; }
    public Class?[] Content { get; }
    
    public Circuit(Class cl, Relator relator)
    {
        Class = cl;
        Relator = relator;
        Content = new Class?[relator.Length + 1];
        Content[0] = Content[relator.Length] = cl;
    }
    
    public int NbUnknowns => Content.Count(c => c is null);
    public bool IsComplete => NbUnknowns == 0;
    
    public (Class? cl1, Class? cl2) UpdateCircuit()
    {
        var leftCircuit = UpdateCircuit(true);
        if (leftCircuit.cl1 is not null)
            return leftCircuit;

        var rightCircuit = UpdateCircuit(false);
        if (rightCircuit.cl1 is not null)
            return rightCircuit;
        
        return (null, null);
    }
    
    public void Substitute(Class cl, Class err)
    {
        for (int i = 0; i < Content.Length; i++)
        {
            var cl0 = Content[i];
            if (cl0 is null)
                continue;

            if (cl0.V == err.V)
                Content[i] = cl;
        }
    }

    private (Class? cl1, Class? cl2) UpdateCircuit(bool leftDir)
    {
        var length = Relator.Length;
        var seq = leftDir ? length.Range() : length.Range(start: length, step: -1);
        var range = seq.Select(k => leftDir ? (k, k + 1, k) : (k, k - 1, k - 1)).ToArray();
        foreach (var (k0, k1, i) in range)
        {
            var cl0 = Content[k0]!;
            var g0 = leftDir ? Relator.Gens[i] : Relator.Gens[i].Invert();
            var cl1 = Content[k1];
            if (cl1 is null)
            {
                if (cl0.HasEdge(g0))
                    Content[k1] = cl0[g0];
                else
                    break;
            }
            else
            {
                if (cl0.HasEdge(g0))
                {
                    var cl2 = cl0[g0];
                    if (cl1.V != cl2.V)
                        return (cl1, cl2);
                }
                else
                    cl0[g0] = cl1;
            }
        }

        return (null, null);
    }

    public (Class? cl, Gen g) FindCandidate()
    {
        var (_, li) = Content.Select((c, k) => (c, k)).FirstOrDefault(e => e.c is null, (null, -1));
        var (_, ri) = Content.Select((c, k) => (c, k)).LastOrDefault(e => e.c is null, (null, -1));
        if (li == -1 && ri == -1)
            return (null, new());
        
        var lc = Content[li - 1]!;
        var rc = Content[ri + 1]!;
        var lg = Relator.Gens[li - 1];
        var rg = Relator.Gens[ri].Invert();
        if (lc.V == rc.V)
        {
            if (lg.CompareTo(rg) < 1)
                return (lc, lg);
            else
                return (rc, rg);
        }
        else if (lc.V < rc.V)
            return (lc, lg);
        else
            return (rc, rg);
    }
}
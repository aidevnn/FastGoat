namespace FastGoat.UserGroup.Words.Tools;

public readonly partial struct Circuit
{
    public List<(Class i, Gen s)> NotColoured
    {
        get
        {
            var gens = Relator.Gens;
            return Content.Select((c, k) => (c!, k)).SkipLast(1)
                .Select(e => (i: e.Item1, s: gens[e.k]))
                .Where(e => !e.i.Coloured[e.s])
                .ToList();
        }
    }
    
}
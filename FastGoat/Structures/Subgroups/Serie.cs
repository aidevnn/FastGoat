using FastGoat.Commons;

namespace FastGoat.Structures.Subgroups;

public struct Serie<T> : IEquatable<Serie<T>> where T : struct, IElt<T>
{
    public List<SubgroupConjugates<T>> Content { get; }
    public int Placeholder { get; }
    public SerieType SerieType { get; }
    public int Count => Content.Count;

    public Serie(List<SubgroupConjugates<T>> serie, SerieType serieType, int digits)
    {
        SerieType = serieType;
        Content = serie;
        Placeholder = digits;
    }

    public IEnumerable<string> ContentName => Content.Select(s => s.Representative.Name);

    public bool IsRefinementOf(Serie<T> serie) => ContentName.All(n => serie.ContentName.Contains(n));

    public bool AddTo(HashSet<Serie<T>> series)
    {
        var serie = this;
        if (series.Any(s => serie.IsRefinementOf(s)))
            return false;

        series.RemoveWhere(s => s.IsRefinementOf(serie));
        return series.Add(serie);
    }

    public bool Equals(Serie<T> other) => ContentName.SequenceEqual(other.ContentName);

    public override int GetHashCode() => Content.Count;

    public override string ToString()
    {
        var digits = Placeholder;
        var dash = Enumerable.Repeat('-', digits).Glue();
        var fmt = $"{{0,-{digits}}}";
        
        if (SerieType == SerieType.Serie)
            return Content.Select(s => s.Order == 1 ? "C1" : $"{s.FullName} {dash}".Substring(0, digits)).Glue("--> ", fmt).Trim();

        return ContentName.Reverse().Glue(" ⊲  ", fmt).Trim();
    }
}
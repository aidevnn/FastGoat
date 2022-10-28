using System.Collections;
using FastGoat;

namespace FastGoat.UserGroup;

public class Quotient<T> : IGroup<Coset<T>> where T : struct, IElt<T>
{
    public Quotient(ConcreteGroup<T> grG, ConcreteGroup<T> grH)
    {
        G = grG;
        H = grH;
        Hash = (G.Hash, G.Hash, "Quo").GetHashCode();
        var hName = H.Name.Contains(' ') ? $"({H.Name})" : H.Name;
        var gName = G.Name.Contains(' ') ? $"({G.Name})" : G.Name;
        Name = $"{gName}/{hName}";

        var cosets = Group.Cosets(G, H);
        Map = cosets.ToDictionary(a => a.Key, x => new Coset<T>(this, x.Value));
        Elements = Map.Values.ToHashSet();
        var lc = Group.LongestCycles(this, Elements);
        var lcKeys = lc.Select(kp => kp.Key).Ascending().ToArray();
        (_, PseudoGenerators) = Group.UniqueGenerators(this, lcKeys);
    }

    private IEnumerable<Coset<T>> PseudoGenerators { get; }
    private Dictionary<T, Coset<T>> Map { get; }
    private HashSet<Coset<T>> Elements { get; }
    public Coset<T> GetRepresentative(T x) => Map[x];
    public ConcreteGroup<T> G { get; }
    public ConcreteGroup<T> H { get; }
    public IEnumerator<Coset<T>> GetEnumerator() => GetElements().GetEnumerator();
    IEnumerator IEnumerable.GetEnumerator() => GetElements().GetEnumerator();
    public IEnumerable<Coset<T>> GetElements() => Elements;

    public bool Equals(IGroup<Coset<T>>? other) => other?.Hash == Hash;

    public int Hash { get; }
    public string Name { get; }

    public Coset<T> this[params ValueType[] us] => GetRepresentative(G[us]);

    public IEnumerable<Coset<T>> GetGenerators() => PseudoGenerators;

    public Coset<T> Neutral() => GetRepresentative(H.Neutral());

    public Coset<T> Invert(Coset<T> e) => GetRepresentative(H.Invert(e.X));

    public Coset<T> Op(Coset<T> e1, Coset<T> e2) => GetRepresentative(H.Op(e1.X, e2.X));

    public override int GetHashCode() => Hash;
    public override string ToString() => Name;
}
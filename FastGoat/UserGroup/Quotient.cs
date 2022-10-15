using System.Collections;
using FastGoat;

namespace FastGoat.UserGroup;

public class Quotient<T> : IGroup<Representative<T>> where T : struct, IElt<T>
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
        Map = cosets.ToDictionary(a => a.Key, x => new Representative<T>(this, x.Value));
        Elements = Map.Values.ToHashSet();
        var lc = Group.LongestCycles(this, Elements);
        var lcKeys = lc.Select(kp => kp.Key).Ascending().ToArray();
        (_, PseudoGenerators) = Group.UniqueGenerators(this, lcKeys);
    }

    private IEnumerable<Representative<T>> PseudoGenerators { get; }
    private Dictionary<T, Representative<T>> Map { get; }
    private HashSet<Representative<T>> Elements { get; }
    public Representative<T> GetRepresentative(T x) => Map[x];
    public ConcreteGroup<T> G { get; }
    public ConcreteGroup<T> H { get; }
    public IEnumerator<Representative<T>> GetEnumerator() => GetElements().GetEnumerator();
    IEnumerator IEnumerable.GetEnumerator() => GetElements().GetEnumerator();
    public IEnumerable<Representative<T>> GetElements() => Elements;

    public bool Equals(IGroup<Representative<T>>? other) => other?.Hash == Hash;

    public int Hash { get; }
    public string Name { get; }

    public Representative<T> this[params ValueType[] us] => GetRepresentative(G[us]);

    public IEnumerable<Representative<T>> GetGenerators() => PseudoGenerators;

    public Representative<T> Neutral() => GetRepresentative(H.Neutral());

    public Representative<T> Invert(Representative<T> e) => GetRepresentative(H.Invert(e.X));

    public Representative<T> Op(Representative<T> e1, Representative<T> e2) => GetRepresentative(H.Op(e1.X, e2.X));

    public override int GetHashCode() => Hash;
    public override string ToString() => Name;
}
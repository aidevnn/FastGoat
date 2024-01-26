using FastGoat.Commons;
using FastGoat.Structures;

namespace FastGoat.UserGroup.Words.Tools;

public struct Relator : IElt<Relator>
{
    public Gen[] Gens { get; }

    public Relator(int idx, Gen[] gens)
    {
        Hash = idx;
        Gens = gens;
    }

    public int Length => Gens.Length;

    public bool Equals(Relator other) => Hash == other.Hash;

    public int CompareTo(Relator other) => Hash.CompareTo(Hash);

    public int Hash { get; }
    public string Format(string fmt) => Gens.Glue("", fmt);

    public override int GetHashCode() => Hash;
    public override string ToString() => Gens.Glue(",");
}
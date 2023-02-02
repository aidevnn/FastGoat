using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.GenericGroup;
using FastGoat.UserGroup.Polynoms;

namespace FastGoat.UserGroup.Characters;

public readonly struct Character<T> : IElt<Character<T>>, IMap<T, Cnf> where T : struct, IElt<T>
{
    public ConjugacyClasses<T> Classes { get; }

    public Character(ConjugacyClasses<T> classes, Dictionary<T, Cnf> map)
    {
        Map = new(map);
        Classes = classes;
        Hash = Map.Aggregate(0, (hash, kp) => (hash, kp.Key.Hash, kp.Value.Hash).GetHashCode());
    }

    public Dictionary<T, Cnf> Map { get; }
    public bool Equals(IMap<T, Cnf>? other) => other?.Hash == Hash;

    public int CompareTo(IMap<T, Cnf>? other) => other is null ? 1 : IMap<T, Cnf>.CompareMap(this, other);

    public int Hash { get; }
    public ConcreteGroup<T> Domain => Classes.Gr;
    public IEnumerable<T> Kernel() => Map.Where(kp => kp.Value.IsZero()).Select(kp => kp.Key);

    public IEnumerable<Cnf> Image() => Map.Values.Distinct();

    public Cnf this[T index] => Map[index];
    public bool Equals(Character<T> other) => IMap<T, Cnf>.EqualiltyMap(this, other);

    public int CompareTo(Character<T> other) => IMap<T, Cnf>.CompareMap(this, other);
    public override int GetHashCode() => Hash;

    public override string ToString()
    {
        var str = new List<string>();
        foreach (var (t, c) in Map)
        {
            if (c.IsZero())
                continue;

            var cstr = $"{c}";
            if (Classes.Gr.Neutral().Equals(t))
            {
                if (cstr.Contains('+'))
                    str.Add($"({cstr})");
                else
                    str.Add($"{cstr}");
            }
            else
            {
                if (cstr.Contains('+'))
                    str.Add($"({cstr})({t})");
                else
                {
                    if (c.Equals(c.One))
                        str.Add($"({t})");
                    else if (c.Equals(-c.One))
                        str.Add($"-({t})");
                    else
                        str.Add($"{cstr}({t})");
                }
            }
        }

        return str.Glue("+");
    }
}
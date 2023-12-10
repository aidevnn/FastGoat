using System.Collections;
using FastGoat.Commons;
using FastGoat.Structures;

namespace FastGoat.UserGroup.Words;

public class WordGroupBase : IGroup<Word>
{
    public WordGroupBase(string relators)
    {
        Hash = Guid.NewGuid().GetHashCode();
        Relators = relators;
        Generators = relators.Where(char.IsLetter).Select(char.ToLower).Distinct().Ascending().ToArray();
        Name = $"WG[{Generators.Glue(",")}]";
        Elements = new() { Neutral() };
    }

    private HashSet<Word> Elements { get; }
    public string Relators { get; }

    private char[] Generators { get; }
    public string Definition => $"< ({Generators.Glue(",")}) | {Relators.Replace(" ", "").Replace(",", ", ").Replace("=", " = ")} >";

    public IEnumerator<Word> GetEnumerator() => GetElements().GetEnumerator();

    IEnumerator IEnumerable.GetEnumerator() => GetElements().GetEnumerator();

    public bool Equals(IGroup<Word>? other) => other?.Hash == Hash;

    public int Hash { get; }
    public string Name { get; }

    public Word this[string s]
    {
        get
        {
            var s0 = StringExt.ReducedWordForm2(s);
            return this[s0.Cast<ValueType>().ToArray()];
        }
    }

    public Word this[params ValueType[] us]
    {
        get
        {
            if (!us.All(c => c is char c1 && Generators.Contains(char.ToLower(c1))))
                throw new GroupException(GroupExceptionType.GroupDef);

            var s0 = StringExt.ReducedWordForm2(us.Glue());
            return new(this, s0);
        }
    }

    public IEnumerable<Word> GetElements() => Elements;

    public IEnumerable<Word> GetGenerators() => Generators.Select(c => new Word(this, new[] { c })).ToArray();

    public Word Neutral() => new(this);

    public Word Invert(Word e) => new(this, e.Get().Revert());

    public Word Op(Word e1, Word e2) => new(this, e1.Get().Add(e2.Get()));

    public override int GetHashCode() => Hash;
    public override string ToString() => Name;
}
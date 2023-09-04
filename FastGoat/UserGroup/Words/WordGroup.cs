using FastGoat.Structures;
using FastGoat.Structures.GenericGroup;
using FastGoat.UserGroup.Words.ToddCoxeter;

namespace FastGoat.UserGroup.Words;

public class WordGroup : ConcreteGroup<Word>
{
    public WordGroup(string name, WordGroupBase wg) : base(name, wg, true)
    {
        WGbase = wg;
        OpsTable = ToddCoxeterAlgo.Run(WGbase.Relators);
        OpsTable.BuildTable();
        Elements = OpsTable.Words().Select(s => new Word(wg, s)).ToHashSet();
        ElementsOrders = Group.ElementsOrders(this, Elements);
        PseudoGenerators = new(wg.GetGenerators().ToList());
        GroupType = (Group.IsCommutative(this, PseudoGenerators)
            ? GroupType.AbelianGroup
            : GroupType.NonAbelianGroup);
    }

    public WordGroup(WordGroupBase wg) : this(wg.Name, wg)
    {
    }

    public WordGroup(string relators) : this(new WordGroupBase(relators))
    {
    }

    public WordGroup(string name, string relators) : this(name, new WordGroupBase(relators))
    {
    }

    public WordGroupBase WGbase { get; }
    public string Definition => WGbase.Definition;
    private OpsTable OpsTable { get; }
    public IEnumerable<char> Rewrite(IEnumerable<char> s) => OpsTable.Rewrite(s);

    public Word this[string s]
    {
        get
        {
            var s0 = WGbase[s];
            return new(WGbase, Rewrite(s0.Get()));
        }
    }

    public new Word this[params ValueType[] us]
    {
        get
        {
            var s0 = WGbase[us];
            return new(WGbase, Rewrite(s0.Get()));
        }
    }

    public override Word Neutral() => new(WGbase);

    public override Word Invert(Word e) => new(WGbase, Rewrite(e.Get().Revert()));

    public override Word Op(Word e1, Word e2) => new(WGbase, Rewrite(e1.Get().Add(e2.Get())));
}
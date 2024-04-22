using FastGoat.Structures;
using FastGoat.Structures.GenericGroup;
using FastGoat.UserGroup.Words.Tools;

namespace FastGoat.UserGroup.Words;

public class WordGroup : ConcreteGroup<Word>
{
    public Dictionary<Word, Word> InvertTable { get; }
    public Dictionary<(Word, Word), Word> OpTable { get; }
    private int Ord { get; }
    public WordGroup(string name, WordGroupBase wg) : base(name, wg, true)
    {
        WGbase = wg;
        Graph = Graph.Run(WGbase.Relators);
        Elements = Graph.Words().Select(s => new Word(wg, s)).ToHashSet();
        Ord = Elements.Count;
        if (Ord < Group.GetStorageCapacity())
        {
            InvertTable = new(2 * Ord);
            OpTable = new(2 * Ord * Ord);
        }
        else
        {
            InvertTable = new();
            OpTable = new();
        }
        
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
    private Graph Graph { get; }
    public IEnumerable<char> Rewrite(IEnumerable<char> s) => Graph.Rewrite(s);

    public bool CheckHomomorphism<T>(ConcreteGroup<T> g, Dictionary<char, T> map) where T : struct, IElt<T>
    {
        return Graph.CheckHomomorphism(g, map);
    }

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

    public override Word Invert(Word e)
    {
        if (Ord >= Group.GetStorageCapacity())
            return new(WGbase, Rewrite(e.Get().Revert()));
        else
        {
            if (InvertTable.TryGetValue(e, out var r))
                return r;

            var e0 = new Word(WGbase, Rewrite(e.Get()));
            var ei = InvertTable[e0] = new(WGbase, Rewrite(e.Get().Revert()));
            return ei;
        }
    }

    public override Word Op(Word e1, Word e2)
    {
        if (Ord >= Group.GetStorageCapacity())
            return new(WGbase, Rewrite(e1.Get().Add(e2.Get())));
        else
        {
            var e12 = (e1, e2);
            if (OpTable.TryGetValue(e12, out var r))
                return r;

            var e120 = (new Word(WGbase, Rewrite(e1.Get())), new Word(WGbase, Rewrite(e2.Get())));
            var e3 = OpTable[e120] = new(WGbase, Rewrite(e1.Get().Add(e2.Get())));
            return e3;
        }
        
        return new(WGbase, Rewrite(e1.Get().Add(e2.Get())));
    }
}
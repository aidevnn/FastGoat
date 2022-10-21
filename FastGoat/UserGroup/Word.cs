namespace FastGoat.UserGroup;

public struct Word : IElt<Word>
{
    private string _word;
    public WordGroup WGroup { get; }
    public string Get() => _word;
    public Word(WordGroup wg)
    {
        WGroup = wg;
        _word = "";
        Hash = (_word, wg).GetHashCode();
    }

    public Word(WordGroup wg, IEnumerable<char> letters)
    {
        WGroup = wg;
        _word = letters.Glue();
        Hash = (_word, wg).GetHashCode();
    }

    public bool Equals(Word other) => Hash == other.Hash;

    public int CompareTo(Word other)
    {
        var compLength = _word.Length.CompareTo(other._word.Length);
        if (compLength != 0)
            return compLength;

        return _word.SequenceCompareTo(other._word);
    }

    public int Hash { get; }
    public IGroup<Word> BaseGroup => WGroup;
    public override int GetHashCode() => Hash;
    public override string ToString()
    {
        return _word.Length == 0 ? "()" : Group.ReducedWordForm1(_word.Glue());
    }
}
namespace FastGoat.Structures.VecSpace;

public readonly struct KAut<K> : IElt<KAut<K>> where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
{
    public KAutGroup<K> KAutGroup { get; }
    public EPoly<K> E { get; }
    
    public KAut(KAutGroup<K> kaut, EPoly<K> e)
    {
        if (!kaut.F.Equals(e.F))
            throw new GroupException(GroupExceptionType.GroupDef);
        
        E = e;
        KAutGroup = kaut;
        Hash = (e.Hash, KAutGroup.Hash).GetHashCode();
    }

    public KAut(EPoly<K> e)
    {
        E = e;
        KAutGroup = new KAutGroup<K>(e.F);
        Hash = (e.Hash, KAutGroup.Hash).GetHashCode();
    }

    public bool Equals(KAut<K> other) => E.Equals(other.E);

    public int CompareTo(KAut<K> other) => E.CompareTo(other.E);

    public int Hash { get; }
    public override int GetHashCode() => Hash;
    public override string ToString() => E.ToString();

    public static implicit operator KAut<K>(EPoly<K> e) => new(e);
    public static implicit operator EPoly<K>(KAut<K> e) => e.E;
}
namespace FastGoat.Structures.VecSpace;

public class TriVarFracSimplifier<K> where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
{
    public virtual TriVarFrac<K> Apply(TriVarPoly<K> num, TriVarPoly<K> denom) => new(this, num, denom);
    public virtual bool IsDivZero(TriVarFrac<K> P) => P.IsZero();
    public override string ToString() => $"TriVarFracSimplifier[{typeof(K)}][{GetHashCode()}]";
}
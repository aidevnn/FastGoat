namespace FastGoat.Structures;

public interface IRingElt<T> : IEquatable<T>, IComparable<T> where T : IElt<T>, IRingElt<T>
{
    bool IsZero();
    T Zero { get; }
    T One { get; }
    T Add(T e);
    T Sub(T e);
    T Opp();
    T Mul(T e);
    (T quo, T rem) Div(T e);
    T Mul(int k);
    T Pow(int k);
}

public interface IFieldElt<T> : IEquatable<T>, IComparable<T> where T : IElt<T>, IRingElt<T>
{
    int P { get; }
    T Inv();
}

public interface IVsElt<K, T> : IEquatable<T>, IComparable<T>
    where T : IElt<T>, IRingElt<T>, IVsElt<K, T> 
    where K : IElt<K>, IRingElt<K>, IFieldElt<K>
{
    int P { get; }
    T KMul(K k);
}
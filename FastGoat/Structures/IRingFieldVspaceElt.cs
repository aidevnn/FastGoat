using System.Numerics;

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
    static abstract T operator +(T a, T b);
    static abstract T operator +(int a, T b);
    static abstract T operator +(T a, int b);
    static abstract T operator -(T a, T b);
    static abstract T operator -(int a, T b);
    static abstract T operator -(T a, int b);
    static abstract T operator *(T a, T b);
    static abstract T operator *(int a, T b);
    static abstract T operator *(T a, int b);
    static abstract T operator /(T a, T b);
    static abstract T operator /(T a, int b);
    int P { get; }
}

public interface IFieldElt<T> : IEquatable<T>, IComparable<T> where T : IElt<T>, IRingElt<T>, IFieldElt<T>
{
    T Inv();
    static abstract T operator /(int a, T b);
}

public interface IVsElt<K, T> : IEquatable<T>, IComparable<T>
    where T : IElt<T>, IRingElt<T>, IVsElt<K, T>
    where K : IElt<K>, IRingElt<K>, IFieldElt<K>
{
    T KMul(K k);
    static abstract T operator +(T a, K b);
    static abstract T operator +(K a, T b);
    static abstract T operator -(T a, K b);
    static abstract T operator -(K a, T b);
    static abstract T operator *(T a, K b);
    static abstract T operator *(K a, T b);
    static abstract T operator /(T a, K b);
}

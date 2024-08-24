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
    static abstract T operator -(T a);
    static abstract T operator -(T a, T b);
    static abstract T operator -(int a, T b);
    static abstract T operator -(T a, int b);
    static abstract T operator *(T a, T b);
    static abstract T operator *(int a, T b);
    static abstract T operator *(T a, int b);
    static abstract T operator /(T a, T b);
    static abstract T operator /(T a, int b);
}

public interface IFieldElt<T> : IEquatable<T>, IComparable<T> where T : IElt<T>, IRingElt<T>, IFieldElt<T>
{
    int P { get; }
    T Inv();
    bool Invertible();
    static abstract T operator /(int a, T b);
    static abstract double Abs(T t);
    static abstract bool IsValuedField { get; }
}

public interface IVsElt<K, T> : IEquatable<T>, IComparable<T>
    where T : IElt<T>, IRingElt<T>, IVsElt<K, T>
    where K : IElt<K>, IRingElt<K>, IFieldElt<K>
{
    K KZero { get; }
    K KOne { get; }
    T KMul(K k);
    static abstract T operator +(T a, K b);
    static abstract T operator +(K a, T b);
    static abstract T operator -(T a, K b);
    static abstract T operator -(K a, T b);
    static abstract T operator *(T a, K b);
    static abstract T operator *(K a, T b);
    static abstract T operator /(T a, K b);
}

public interface IFloatElt<K> : IEquatable<K>, IComparable<K> where K : IElt<K>, IRingElt<K>, IFieldElt<K>, IFloatElt<K>
{
    public K RoundEven { get; }
    public K Absolute { get; }
    public int Sign { get; }
    public static abstract K Sqrt(K r);
    public static abstract K NthRoot(K r, int n);
}

public interface IFixedPrecisionElt<K> : IEquatable<K>, IComparable<K> 
    where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>, IFloatElt<K>, IFixedPrecisionElt<K>
{
    public static abstract K From<T>(T e) where T : IElt<T>, IRingElt<T>, IFieldElt<T>, IFloatElt<T>;
    public static abstract int Digits { get; }
    public static abstract double Eps { get; }
    static abstract bool operator ==(K a, K b);
    static abstract bool operator !=(K a, K b);
    static abstract bool operator <(K a, K b);
    static abstract bool operator >(K a, K b);
    static abstract bool operator <=(K a, K b);
    static abstract bool operator >=(K a, K b);
    static abstract K Min(K a, K b);
    static abstract K Max(K a, K b);
}
namespace FastGoat.Structures;

public interface IGmoduleNElt<N, G, T> : IEquatable<T>, IComparable<T>
    where N : struct, IElt<N>
    where G : struct, IElt<G>
    where T : struct, IElt<T>, IGmoduleNElt<N, G, T>
{
    T Add(T n);
    T Sub(T n);
    T Opp();
    T Act(int k);

    T Act(G a);

    static abstract T operator +(T a, T b);
    static abstract T operator -(T a, T b);
    static abstract T operator -(T a);
    static abstract T operator *(G a, T b);
    static abstract T operator *(int k, T b);
}
using System.Collections;
using FastGoat.Commons;
using FastGoat.Structures.CartesianProduct;

namespace FastGoat.Structures.VecSpace;

public readonly struct Vec<A> : IEnumerable<A>, IElt<Vec<A>>, IRingElt<Vec<A>>, IModuleElt<A, Vec<A>>
    where A : IElt<A>, IRingElt<A>
{
    public Ep<A> V { get; }
    public int Length => V.Ei.Length;

    public Vec(A[] elts)
    {
        V = new(elts);
        Hash = HashCode.Combine("vec", V);
        if (elts.Length == 0)
            throw new GroupException(GroupExceptionType.GroupDef);

        KZero = elts[0].Zero;
        KOne = elts[0].One;
    }

    public bool Equals(Vec<A> other) => V.Equals(other.V);

    public int CompareTo(Vec<A> other) => V.CompareTo(other.V);

    public A KZero { get; }
    public A KOne { get; }
    public Vec<A> KMul(A k) => new(V.Ei.Select(e => e * k).ToArray());
    public int Hash { get; }
    public bool IsZero() => V.Ei.All(e => e.IsZero());

    public Vec<A> Zero => new(V.Ei.Select(e => e.Zero).ToArray());
    public Vec<A> One => new(V.Ei.Select(e => e.One).ToArray());

    public Vec<A> Add(Vec<A> e)
    {
        if (e.Length != Length)
        {
            V.Ei.Println("first");
            e.Println("second");
            throw new GroupException(GroupExceptionType.GroupDef);
        }

        return new(V.Ei.Select((ei, i) => ei + e[i]).ToArray());
    }

    public Vec<A> Sub(Vec<A> e)
    {
        if (e.Length != Length)
            throw new GroupException(GroupExceptionType.GroupDef);

        return new(V.Ei.Select((ei, i) => ei - e[i]).ToArray());
    }

    public Vec<A> Opp() => new(V.Ei.Select(ei => -ei).ToArray());

    public Vec<A> Mul(Vec<A> e)
    {
        if (e.Length != Length)
            throw new GroupException(GroupExceptionType.GroupDef);

        return new(V.Ei.Select((ei, i) => ei * e[i]).ToArray());
    }

    public (Vec<A> quo, Vec<A> rem) Div(Vec<A> e)
    {
        throw new NotImplementedException();
    }

    public Vec<A> Mul(int k) => KMul(k * KOne);

    public A InnerProd(Vec<A> b) => Mul(b).Aggregate((ei, ej) => ei + ej);

    public Vec<A> Pow(int k) => new(V.Ei.Select(ei => ei.Pow(k)).ToArray());
    
    public Vec<U> AddA<U>(Vec<U> u) where U : IElt<U>, IRingElt<U>, IModuleElt<A, U>
    {
        if (u.Length != Length)
            throw new GroupException(GroupExceptionType.GroupDef);

        return new(V.Ei.Zip(u).Select(e => e.First + e.Second).ToArray());
    }

    public Vec<U> SubA<U>(Vec<U> u) where U : IElt<U>, IRingElt<U>, IModuleElt<A, U>
    {
        if (u.Length != Length)
            throw new GroupException(GroupExceptionType.GroupDef);

        return new(V.Ei.Zip(u).Select(e => e.First - e.Second).ToArray());
    }

    public Vec<U> MulA<U>(Vec<U> u) where U : IElt<U>, IRingElt<U>, IModuleElt<A, U>
    {
        if (u.Length != Length)
            throw new GroupException(GroupExceptionType.GroupDef);

        return new(V.Ei.Zip(u).Select(e => e.First * e.Second).ToArray());
    }

    public Vec<U> MulA<T, U>(Vec<U> u) 
        where U : IElt<U>, IRingElt<U>, IModuleElt<T, U> 
        where T : IElt<T>, IRingElt<T>, IModuleElt<A, T>
    {
        if (u.Length != Length)
        {
            V.Ei.Println("first");
            u.Println("second");
            throw new GroupException(GroupExceptionType.GroupDef);
        }

        var t = V;
        var o = u.KOne.KOne;
        var v = Length.SeqLazy().Select(i => (t[i] * o) * u[i]).ToArray();
        return new(v);
    }

    public Vec<U> MulA<T, U>(U u)
        where U : IElt<U>, IRingElt<U>, IModuleElt<T, U>
        where T : IElt<T>, IRingElt<T>, IModuleElt<A, T>
    {
        var t = V;
        var o = u.KOne;
        var v = Length.SeqLazy().Select(i => (t[i] * o) * u).ToArray();
        return new(v);
    }

    public Vec<U> Mul2A<U>(Vec<U> u, int size) where U : IElt<U>, IRingElt<U>, IModuleElt<A, U>
    {
        var idx = Length.SeqLazy().Grid2D(u.Length.SeqLazy()).Where(e => e.t1 + e.t2 <= size)
            .GroupBy(e => e.t1 + e.t2).OrderBy(e => e.Key).ToArray();
        var t = this;
        var v = idx.Select(k => k.Select(e => t[e.t1] * u[e.t2]).Aggregate((a0, a1) => a0 + a1)).ToArray();
        return new(v);
    }

    public Vec<U> Mul2A<U>(Vec<U> u) where U : IElt<U>, IRingElt<U>, IModuleElt<A, U>
        => Mul2A(u, Length + u.Length);

    public Vec<A> Mul2(Vec<A> u, int size)
    {
        var idx = Length.SeqLazy().Grid2D(u.Length.SeqLazy()).Where(e => e.t1 + e.t2 <= size)
            .GroupBy(e => e.t1 + e.t2).OrderBy(e => e.Key).ToArray();
        var t = this;
        var v = idx.Select(k => k.Select(e => t[e.t1] * u[e.t2]).Aggregate((a0, a1) => a0 + a1)).ToArray();
        return new(v);
    }

    public Vec<A> Mul2(Vec<A> u) => Mul2(u, Length + u.Length);

    public A this[int idx] => idx < 0 || idx >= Length ? KZero : V[idx];

    public IEnumerator<A> GetEnumerator() => V.Ei.Select(e => e).GetEnumerator();

    IEnumerator IEnumerable.GetEnumerator()
    {
        return GetEnumerator();
    }

    public override int GetHashCode() => Hash;

    public override string ToString() => $"{V}";
    public A Sum(Func<A, A> f) => V.Ei.Aggregate(KZero, (acc, e) => acc + f(e));

    public A Sum() => Sum(a => a);

    public static Vec<A> operator +(Vec<A> a, Vec<A> b) => a.Add(b);

    public static Vec<A> operator +(int a, Vec<A> b) => new(b.Select(e => e + a).ToArray());

    public static Vec<A> operator +(Vec<A> a, int b) => b + a;

    public static Vec<A> operator -(Vec<A> a) => a.Opp();

    public static Vec<A> operator -(Vec<A> a, Vec<A> b) => a.Sub(b);

    public static Vec<A> operator -(int a, Vec<A> b) => (-a) + b;

    public static Vec<A> operator -(Vec<A> a, int b) => a + (-b);

    public static Vec<A> operator *(Vec<A> a, Vec<A> b) => a.Mul(b);

    public static Vec<A> operator *(int a, Vec<A> b) => b.Mul(a);

    public static Vec<A> operator *(Vec<A> a, int b) => a.Mul(b);

    public static Vec<A> operator /(Vec<A> a, Vec<A> b) => a.Div(b).quo;

    public static Vec<A> operator /(Vec<A> a, int b) => new(a.Select(e => e / b).ToArray());

    public static Vec<A> operator +(Vec<A> a, A b) => new(a.Select(e => e * b).ToArray());

    public static Vec<A> operator +(A a, Vec<A> b) => b + a;

    public static Vec<A> operator -(Vec<A> a, A b) => a + b.Opp();

    public static Vec<A> operator -(A a, Vec<A> b) => b + a.Opp();

    public static Vec<A> operator *(Vec<A> a, A b) => a.KMul(b);

    public static Vec<A> operator *(A a, Vec<A> b) => b * a;
}
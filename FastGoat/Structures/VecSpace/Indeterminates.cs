using System.CodeDom;
using System.Collections;
using FastGoat.Commons;

namespace FastGoat.Structures.VecSpace;

public class Indeterminates<T> : IEnumerable<T>, IEquatable<Indeterminates<T>> where T : IElt<T>
{
    public T[] Content { get; private set; }

    public bool Graded { get; private set; }
    public bool Reverse { get; private set; }

    public Indeterminates()
    {
        throw new Exception();
    }

    public Indeterminates(T[] arr, bool g, bool r)
    {
        Content = arr.ToArray();
        if (Content.Length == 0)
            throw new ArgumentException();

        // if(Content.Distinct().Count() != arr.Count()) {} // TODO warning

        (Graded, Reverse) = (g, r);
    }

    public Indeterminates(params T[] arr) : this(arr, true, false)
    {
    }

    public Indeterminates(MonomOrder order, params T[] arr)
    {
        Content = arr.ToArray();
        if (Content.Length == 0)
            throw new ArgumentException();

        // if(Content.Distinct().Count() != arr.Count()) {} // TODO warning

        SetOrder(order);
    }

    public int Length => Content.Length;
    public MonomOrder Order =>
        (Graduate: Graded, Reverse) switch
        {
            (false, false) => MonomOrder.Lex,
            (true, false) => MonomOrder.GrLex,
            (false, true) => MonomOrder.RevLex,
            _ => MonomOrder.GrevLex
        };

    public void SetOrder(MonomOrder order)
    {
        (Graded, Reverse) = order switch
        {
            MonomOrder.Lex => (false, false),
            MonomOrder.GrLex => (true, false),
            MonomOrder.RevLex => (false, true),
            MonomOrder.GrevLex => (true, true),
            _ => (Graded, Reverse)
        };
    }

    public void Extend(params T[] xi)
    {
        Content = Content.Concat(xi).ToArray();
        if (Content.Length != Content.Distinct().Count())
            throw new();
    }

    public void Permute(int[] perm)
    {
        if (Enumerable.Range(0, Content.Length).SequenceEqual(perm.Order()))
        {
            var cont = perm.Select(i => Content[i]).ToArray();
            Content = cont;
        }
    }

    public bool Contains(T t) => Content.Contains(t);
    public bool Contains(string t) => Content.Any(xi => xi.ToString()!.Equals(t));

    public int Hash => Content.Aggregate(0, (acc, a) => (acc, a.Hash).GetHashCode());

    public IEnumerator<T> GetEnumerator()
    {
        if (Reverse)
            return Content.Reverse().GetEnumerator();
        else
            return Content.AsEnumerable().GetEnumerator();
    }

    IEnumerator IEnumerable.GetEnumerator()
    {
        return GetEnumerator();
    }

    public T this[int index] => Content[index];
    public T this[string name] => Content.First(xi => xi.ToString()!.Equals(name));
    public override int GetHashCode() => Hash;
    public override string ToString() => $"[{this.Glue(",")}]({Order})";

    public bool Equals(Indeterminates<T>? other)
    {
        if (ReferenceEquals(null, other)) return false;
        if (ReferenceEquals(this, other)) return true;
        return Content.Equals(other.Content) && Graded == other.Graded && Reverse == other.Reverse && Hash == other.Hash;
    }

    public override bool Equals(object? obj)
    {
        if (ReferenceEquals(null, obj)) return false;
        if (ReferenceEquals(this, obj)) return true;
        if (obj.GetType() != this.GetType()) return false;
        return Equals((Indeterminates<T>)obj);
    }
    
    public void Deconstruct(out T a, out T b)
    {
        (a, b) = (this[0], this[1]);
    }

    public void Deconstruct(out T a, out T b, out T c)
    {
        (a, b, c) = (this[0], this[1], this[2]);
    }

    public void Deconstruct(out T a, out T b, out T c, out T d)
    {
        (a, b, c, d) = (this[0], this[1], this[2], this[3]);
    }

    public void Deconstruct(out T a, out T b, out T c, out T d, out T e)
    {
        (a, b, c, d, e) = (this[0], this[1], this[2], this[3], this[4]);
    }

}
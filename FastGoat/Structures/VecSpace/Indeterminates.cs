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

        if (Content.Distinct().Count() != arr.Count())
            throw new();

        SetOrder(order);
    }

    public int Length => Content.Length;

    public MonomOrder Order =>
        (Graded, Reverse) switch
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

    public void Revert() => Reverse = !Reverse;

    public void ExtendPrepend(params T[] xi)
    {
        Content = xi.Concat(Content).ToArray();
        if (Content.Length != Content.Distinct().Count())
            throw new();
    }

    public void ExtendAppend(params T[] xi)
    {
        Content = Content.Concat(xi).ToArray();
        if (Content.Length != Content.Distinct().Count())
            throw new();
    }

    public void Remove(params T[] xi)
    {
        Content = Content.Except(xi).ToArray();
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

    public T this[int index] => Reverse ? Content[^(index + 1)] : Content[index];
    public T this[string name] => Content.First(xi => xi.ToString()!.Equals(name));

    public int this[T e]
    {
        get
        {
            if (Reverse)
                return Content.Reverse().ToList().FindIndex(d => d.Equals(e));
            
            return Content.ToList().FindIndex(d => d.Equals(e));
        }
    }
    public override int GetHashCode() => Hash;
    public override string ToString() => $"[{Content.Glue(",")}]({Order})";

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
}
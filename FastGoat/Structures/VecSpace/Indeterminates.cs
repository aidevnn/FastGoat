using System.Collections;
using FastGoat.Commons;

namespace FastGoat.Structures.VecSpace;

public class Indeterminates<T> : IEnumerable<T> where T : IElt<T>
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
        Hash = Content.Aggregate(0, (acc, a) => (acc, a.Hash).GetHashCode());
    }

    public Indeterminates(params T[] arr) : this(arr, true, false)
    {
    }

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

    public void Permute(int[] perm)
    {
        if (IntExt.CheckTable(Content.Length, perm))
        {
            var cont = perm.Select(i => Content[i]).ToArray();
            Content = cont;
        }
    }

    public bool Contains(T t) => Content.Contains(t);

    public int Hash { get; }

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
    public override int GetHashCode() => Hash;
    public override string ToString() => $"[{Content.Glue(",")}]({Order})";
}
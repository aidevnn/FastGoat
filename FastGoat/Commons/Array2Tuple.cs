namespace FastGoat.Commons;

public class Array2Tuple<T>
{
    private IEnumerable<T> array { get; }

    public Array2Tuple(IEnumerable<T> arr)
    {
        array = arr;
    }

    public T this[int index] => array.ElementAt(index);

    public void Deconstruct(out T a0, out T a1)
    {
        (a0, a1) = (this[0], this[1]);
    }

    public void Deconstruct(out T a0, out T a1, out T a2)
    {
        (a0, a1, a2) = (this[0], this[1], this[2]);
    }

    public void Deconstruct(out T a0, out T a1, out T a2, out T a3)
    {
        (a0, a1, a2, a3) = (this[0], this[1], this[2], this[3]);
    }

    public void Deconstruct(out T a0, out T a1, out T a2, out T a3, out T a4)
    {
        (a0, a1, a2, a3, a4) = (this[0], this[1], this[2], this[3], this[4]);
    }

    public void Deconstruct(out T a0, out T a1, out T a2, out T a3, out T a4, out T a5)
    {
        (a0, a1, a2, a3, a4, a5) = (this[0], this[1], this[2], this[3], this[4], this[5]);
    }

    public void Deconstruct(out T a0, out T a1, out T a2, out T a3, out T a4, out T a5, out T a6)
    {
        (a0, a1, a2, a3, a4, a5, a6) = (this[0], this[1], this[2], this[3], this[4], this[5], this[6]);
    }

    public void Deconstruct(out T a0, out T a1, out T a2, out T a3, out T a4, out T a5, out T a6, out T a7)
    {
        (a0, a1, a2, a3, a4, a5, a6, a7) = (this[0], this[1], this[2], this[3], this[4], this[5], this[6], this[7]);
    }

    public void Deconstruct(out T a0, out T a1, out T a2, out T a3, out T a4, out T a5, out T a6, out T a7, out T a8)
    {
        (a0, a1, a2, a3, a4, a5, a6, a7, a8) = (this[0], this[1], this[2], this[3], this[4], this[5], this[6], this[7], this[8]);
    }

    public override string ToString()
    {
        return $"({array.Glue(", ")})";
    }
}
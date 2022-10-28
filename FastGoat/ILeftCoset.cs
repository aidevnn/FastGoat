namespace FastGoat;

public interface ILeftCoset<out T, U> : IElt<U> where T : struct, IElt<T> where U : struct, IElt<U>
{
    T X { get; }
    IEnumerable<T> xH { get; }
}
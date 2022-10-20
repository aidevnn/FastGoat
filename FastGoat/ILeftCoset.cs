namespace FastGoat;

public interface ILeftCoset<out T> where T : struct, IElt<T>
{
    T X { get; }
    IEnumerable<T> xH { get; }
}
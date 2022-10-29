namespace FastGoat.Structures;

public interface ILeftCoset<T, U> : IEnumerable<T>, IElt<U> where T : struct, IElt<T> where U : struct, IElt<U>
{
    T X { get; }
    ConcreteGroup<T> G { get; }
    ConcreteGroup<T> H { get; }
}
using FastGoat.Structures.GenericGroup;

namespace FastGoat.Structures;

public interface ICoset<T, U> : IEnumerable<T>, IElt<U> where T : struct, IElt<T> where U : struct, IElt<U>
{
    CosetType CosetType { get; }
    T X { get; }
    ConcreteGroup<T> G { get; }
    ConcreteGroup<T> H { get; }
}
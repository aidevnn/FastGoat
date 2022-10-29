namespace FastGoat.Structures;

public interface IMap<T1, T2> : IElt<IMap<T1, T2>> where T1 : struct, IElt<T1> where T2 : struct, IElt<T2>
{
    ConcreteGroup<T1> Domain { get; }
    IEnumerable<T1> Kernel();
    IEnumerable<T2> Image();
    T2 this[T1 index] { get; }
}
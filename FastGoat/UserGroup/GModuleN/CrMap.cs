using System.Collections;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Integers;

namespace FastGoat.UserGroup.GModuleN;

public readonly struct CrMap<Tn, Tg> : IEnumerable<KeyValuePair<Ep<Tg>, ZNElt<Tn,Tg>>> 
    where Tg : struct, IElt<Tg> 
    where Tn : struct, IElt<Tn>
{
    private Dictionary<Ep<Tg>, ZNElt<Tn, Tg>> Map { get; }

    public CrMap()
    {
        Map = new();
    }
    
    public CrMap(Dictionary<Ep<Tg>, ZNElt<Tn, Tg>> map)
    {
        Map = map;
    }
    
    public int R => Map.Count == 0 ? 0 : Map.First().Key.Ei.Length;
    public bool IsZero() => Map.Values.All(c => c.IsZero());

    public ZNElt<Tn, Tg> this[Ep<Tg> index] => Map[index];

    public Dictionary<Ep<Tg>, ZNElt<Tn, Tg>>.KeyCollection Keys => Map.Keys;
    public Dictionary<Ep<Tg>, ZNElt<Tn, Tg>>.ValueCollection Values => Map.Values;
    public int Count => Map.Count;

    public IEnumerator<KeyValuePair<Ep<Tg>, ZNElt<Tn, Tg>>> GetEnumerator() => Map.GetEnumerator();

    IEnumerator IEnumerable.GetEnumerator()
    {
        return GetEnumerator();
    }
}
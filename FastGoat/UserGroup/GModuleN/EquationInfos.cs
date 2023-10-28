using System.Collections.ObjectModel;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Integers;

namespace FastGoat.UserGroup.GModuleN;

public record EquationInfos(Polynomial<ZnInt, Xi> equation, 
    int mod, 
    ReadOnlyDictionary<int, int> invertibles, 
    ReadOnlyDictionary<int, int> orders, 
    Monom<Xi> lm, 
    ZnInt lc);
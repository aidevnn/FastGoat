using System.Collections.ObjectModel;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Integers;

namespace FastGoat.UserGroup.GModuleN;

public record EquationInfos(Polynomial<ZnInt, Xi> equation, int mod, Monom<Xi> lma, ZnInt lca);